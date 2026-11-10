# test-call-scn.R: the pipeline (R/call_scn.R): caller parameters from caller.yaml, input normalisation, call_scn()
# end to end on the sim_norm_f simulation written in the extraction layout (features -> embedding ->
# caller_normal), the output writer / reader and the summary, the tumour-mode dispatch to the tumour mode, and the CLI.

# centromere table of sim_genome(): p = [0, 0.4 x length)
cs_centromeres = function(bins) {
  cs = bins[, .(len = max(end)), by = chrom]
  cs[, .(chrom, cstart = 1, cen_left = 0.4 * len, cen_right = 0.4 * len, cend = len)]
}

call_cache = new.env()
sim_dir = function(name) {
  key = paste0(name, "_dir")
  if (is.null(call_cache[[key]])) {
    d = file.path(tempdir(), paste0("callscn_", name))
    unlink(d, recursive = TRUE)
    invisible(capture.output(suppressMessages(write_sim_dataset(get(name), d))))
    call_cache[[key]] = d
  }
  call_cache[[key]]
}
norm_result = function() {
  if (is.null(call_cache$norm)) {
    t0 = Sys.time()
    d = sim_dir("sim_norm_f")
    call_cache$norm = call_scn(d, "normal", sex = "female", centromeres = cs_centromeres(sim_norm_f$scn$bins), verbose = FALSE)
    call_cache$norm_min = as.numeric(difftime(Sys.time(), t0, units = "mins"))
  }
  call_cache$norm
}

test_that("caller_params reads caller.yaml (= caller_defaults) and applies overrides", {
  f = system.file("config", "caller.yaml", package = "scphasedcn")
  expect_true(file.exists(f))
  pf = caller_params(f)
  expect_identical(pf[names(caller_defaults())], caller_defaults())
  expect_setequal(names(pf), c(names(caller_defaults()), "normal", "tumour"))
  expect_identical(caller_params(), pf)
  # the normal: / tumour: blocks are the full caller_normal() argument defaults / clone_defaults()
  expect_identical(pf$tumour, clone_defaults())
  expect_identical(pf$normal, normal_arg_defaults())
  expect_equal(pf$normal$n_draws, 20000)
  # the fallback parser (no yaml package) gives the same parameters
  y = parse_simple_yaml(readLines(f))
  p = caller_defaults()
  for (k in intersect(names(y), names(p))) { v = coerce_param(y[[k]], p[[k]]); if (!is.null(v)) p[[k]] = v }
  expect_equal(p, caller_defaults())
  expect_equal(block_params(y$tumour, clone_defaults(), "tumour"), clone_defaults())
  expect_setequal(names(y), c(names(caller_defaults()), "normal", "tumour"))
  p2 = caller_params(f, overrides = list(null_q = 0.99, shifts = c(loss = -1, gain = 0.6), extra_key = "x"))
  expect_equal(p2$null_q, 0.99)
  expect_equal(p2$shifts, c(loss = -1, gain = 0.6))
  expect_equal(p2$extra_key, "x")
  tf = tempfile(fileext = ".yaml")
  writeLines(c("null_q: 0.995", "n_pool: 10", "shifts:", "  loss: -0.9", "  gain: 0.5", "tabix: \"\"", "bogus: 1"), tf)
  expect_warning(p3 <- caller_params(tf), "bogus")
  expect_equal(p3$null_q, 0.995); expect_equal(p3$n_pool, 10); expect_equal(p3$shifts, c(loss = -0.9, gain = 0.5))
  expect_equal(p3$tabix, caller_defaults()$tabix)
  expect_null(p3$tumour); expect_null(p3$normal)
  # blocks: coerced to the defaults' types, a logical option takes a string (split_pass "disagree"), unknown keys warn
  tf3 = tempfile(fileext = ".yaml")
  writeLines(c("seed: 3", "tumour:", "  split_pass: \"disagree\"", "  em_rounds: 3", "  confirm_null: generic",
               "normal:", "  n_draws: 500", "  pool_loo: false"), tf3)
  p4 = caller_params(tf3)
  expect_equal(p4$seed, 3)
  expect_identical(p4$tumour$split_pass, "disagree")
  expect_identical(p4$tumour$em_rounds, 3)
  expect_identical(p4$tumour$confirm_null, "generic")
  expect_identical(p4$tumour$ref_null, clone_defaults()$ref_null)
  expect_equal(p4$normal$n_draws, 500); expect_false(p4$normal$pool_loo); expect_true(p4$normal$pool_standardize)
  writeLines(c("tumour:", "  nope: 1"), tf3)
  expect_warning(caller_params(tf3), "nope")
  # params.yaml written by write_scphasedcn reads back to the same parameters
  tf2 = tempfile(fileext = ".yaml")
  writeLines(yaml_lines(caller_defaults()), tf2)
  expect_equal(caller_params(tf2), caller_defaults())
})

test_that("multiplet_table normalises vectors, files and multiplet score tables", {
  expect_equal(multiplet_table(c("a", "b")), data.table(cell = c("a", "b"), prob = 1))
  expect_null(multiplet_table(NULL))
  e0 = data.table(cell = c("a", "b", "c"), clean_singlet = c(TRUE, FALSE, NA))
  expect_equal(multiplet_table(e0)$prob, c(0, 1, 0))
  expect_equal(multiplet_table(data.table(cell = "a", prob = 0.3))$prob, 0.3)
  tf = tempfile(); writeLines(c("x", "y"), tf)
  expect_equal(multiplet_table(tf)$cell, c("x", "y"))
})

test_that("call_scn runs end to end on sim_norm_f (normal mode) in under 2 minutes", {
  res = norm_result()
  expect_s3_class(res, "scphasedcn")
  expect_equal(res$mode, "normal")
  expect_lt(call_cache$norm_min, 2)
  expect_equal(nrow(res$cells), nrow(res$feature_info$cells))
  expect_true(all(c("cell", "qc", "stratum", "fold", "null_cell", "n_calls") %in% names(res$cells)))
  check_calls(res$calls); check_thresholds(res$thresholds)
  expect_true(all(c("cell", "pool") %in% res$calls$level))
  expect_true(any(res$pools$kind == "knn"))                         # chromatin embedding -> kNN pools
  expect_equal(nrow(res$embedding), nrow(res$cells))
  # truth: every reportable chr8 gain / chr1p loss / chr21 gain call is on a carrier, and most carriers are found
  tr = sim_norm_f$truth
  cells = sim_norm_f$scn$cells$cell
  car = lapply(tr$carriers, function(n) unique(cells[tr$nuclei$cell[n]]))
  cl = res$calls[level == "cell" & called & reportable]
  for (k in list(list("chr1p", "loss", 1), list("chr8", "gain", 2), list("chr21", "gain", 4))) {
    got = unique(cl[region == k[[1]] & state == k[[2]], cell])
    expect_true(all(got %in% car[[k[[3]]]]), info = k[[1]])
    expect_gte(length(intersect(got, car[[k[[3]]]])) / length(car[[k[[3]]]]), 0.7)
  }
  # null false calls: rates per class x stratum + "all", per-cell false events
  expect_true(all(c("class", "stratum", "n_null", "n_called", "per_1000") %in% names(res$false_calls)))
  expect_true("all" %in% res$false_calls$stratum)
  expect_equal(res$false_events$state, c("any", "loss", "gain", "cnloh"))
  expect_true(any(grepl("^call_scn normal", res$notes)))
})

test_that("write_scphasedcn writes every file, read_scphasedcn round-trips, scphasedcn_summary has the headline rows", {
  res = norm_result()
  out = file.path(tempdir(), "callscn_out")
  unlink(out, recursive = TRUE)
  files = write_scphasedcn(res, out, write_scores = TRUE)
  want = c("cells.tsv", "calls.tsv", "pooled_calls.tsv", "events.tsv", "thresholds.tsv", "null_fit.tsv", "false_calls.tsv",
           "false_events.tsv", "regions.tsv", "summary.tsv", "result.rds", "params.yaml", "NOTES.md", "scores.tsv.gz")
  expect_true(all(file.exists(file.path(out, want))))
  expect_setequal(names(files), want)
  ct = fread(file.path(out, "cells.tsv"))
  expect_equal(names(ct)[1:6], c("cell", "qc", "stratum", "fold", "null_cell", "n_calls"))
  expect_equal(nrow(ct), nrow(res$cells))
  calls = fread(file.path(out, "calls.tsv"))
  expect_true(all(calls$called) && all(calls$level == "cell"))
  expect_equal(nrow(calls), res$calls[level == "cell" & called == TRUE, .N])
  expect_equal(nrow(fread(file.path(out, "pooled_calls.tsv"))), res$calls[level == "pool" & called == TRUE, .N])
  nf = fread(file.path(out, "null_fit.tsv"))
  expect_true(all(c("depth_t", "allelic_rho") %in% nf$kind))
  expect_equal(nf[kind == "allelic_rho", rho], res$null_fit$rho)
  expect_equal(caller_params(file.path(out, "params.yaml")), res$params)
  # params records the caller_normal() arguments of the run (defaults filled in)
  expect_equal(res$params$normal$n_draws, 20000)
  expect_true(res$params$normal$pool_loo)
  expect_setequal(names(res$params$normal), names(normal_arg_defaults()))
  expect_true(any(grepl("## Summary", readLines(file.path(out, "NOTES.md")))))
  back = read_scphasedcn(out)
  expect_s3_class(back, "scphasedcn")
  expect_equal(back$calls, res$calls)
  expect_equal(back$cells, res$cells)
  expect_equal(back$thresholds, res$thresholds)
  expect_equal(back$events, res$events)
  expect_null(back$scores); expect_null(back$features)
  expect_equal(read_scphasedcn(file.path(out, "result.rds"))$mode, "normal")
  expect_error(read_scphasedcn(tempfile()), "no result.rds")
  # a result.rds written under the working name (class "scncall") reads back as scphasedcn
  old = readRDS(file.path(out, "result.rds")); class(old) = "scncall"
  f_old = tempfile(fileext = ".rds"); saveRDS(old, f_old)
  back_old = read_scphasedcn(f_old)
  expect_s3_class(back_old, "scphasedcn")
  expect_equal(back_old$calls, res$calls)
  expect_output(print(old), "<scphasedcn> mode: normal")
  out_old = file.path(tempdir(), "callscn_out_legacy")
  unlink(out_old, recursive = TRUE)
  write_scphasedcn(old, out_old)
  expect_s3_class(readRDS(file.path(out_old, "result.rds")), "scphasedcn")
  saveRDS(list(a = 1), f_old)
  expect_error(read_scphasedcn(f_old), "not a scphasedcn object")

  sm = scphasedcn_summary(res)
  expect_equal(names(sm), c("section", "item", "value", "n", "note"))
  expect_true(all(c("cells", "calls", "false_calls", "false_events", "events") %in% sm$section))
  expect_equal(sm[item == "cells", value], nrow(res$cells))
  expect_equal(sm[item == "null_cells", value], sum(res$cells$null_cell))
  expect_equal(sm[item == "pooled_calls", value], res$calls[level == "pool", sum(called)])
  expect_equal(sm[section == "calls", sum(value)], res$calls[called == TRUE, .N])
  expect_equal(sm[item == "arm_calls_per_null_cell|any", value], res$false_events[state == "any", arm_calls_per_cell])
  expect_true(all(c("pooled_calls|type", "pooled_calls|knn") %in% sm$item))
})

test_that("call_scn takes precomputed features, explicit null cells and multiplets", {
  res = norm_result()
  fe = list(long = copy(res$features), regions = res$regions, cells = res$feature_info$cells, use = res$feature_info$use,
            phasing = res$feature_info$phasing, notes = character())
  nulls = res$cells[null_cell == TRUE, cell][1:120]
  mult = res$cells$cell[1:5]
  r2 = call_scn(sim_dir("sim_norm_f"), "normal", features = fe, null_cells = nulls, multiplets = mult, embedding = FALSE,
                verbose = FALSE, keep_scores = FALSE)
  expect_setequal(r2$cells[null_cell == TRUE, cell], setdiff(nulls, mult))
  expect_true(all(r2$cells[cell %in% mult, multiplet]))
  expect_null(r2$scores)
  expect_false(any(r2$pools$kind == "knn"))
  expect_false(any(r2$pools$cell %in% mult))
  # normal mode: `groups` is ignored with a warning; caller_args override the params' normal: block and are recorded
  prm = caller_defaults(); prm$normal = utils::modifyList(normal_arg_defaults(), list(n_draws = 3000, pool_loo = FALSE))
  expect_warning(r3 <- call_scn(sim_dir("sim_norm_f"), "normal", features = fe, null_cells = nulls, embedding = FALSE, groups = "x",
                                params = prm, caller_args = list(n_draws = 2000), verbose = FALSE, keep_scores = FALSE),
                 "tumour mode only")
  expect_equal(r3$params$normal$n_draws, 2000)
  expect_false(r3$params$normal$pool_loo)
  expect_true(any(grepl("called on s_comb", r3$notes)))
  expect_null(r3$params$tumour)
})

test_that("scn_embedding is reproducible (seeded truncated SVD) and named by cell", {
  res = norm_result()
  d = sim_dir("sim_norm_f")
  e1 = scn_embedding(d, res$cells$cell, res$feature_info$use)
  expect_equal(rownames(e1), res$cells$cell)
  expect_identical(e1, scn_embedding(d, res$cells$cell, res$feature_info$use))
  expect_identical(e1, res$embedding)
})

tum_result = function(params = caller_defaults()) {
  s = sim_tum2
  cm = data.table(cell = s$scn$cells$cell, group = s$scn$cells$cell_type)
  call_scn(sim_dir("sim_tum2"), "tumour", cell_meta = cm, cells = cm$cell, sex = "male", groups = "group", verbose = FALSE,
           centromeres = cs_centromeres(s$scn$bins), params = params)
}

test_that("tumour mode dispatches to the tumour mode", {
  res = tum_result()
  call_cache$tum = res
  expect_s3_class(res, "scphasedcn")
  expect_equal(res$mode, "tumour")
  expect_true(all(c("clone", "posterior") %in% names(res$cells)))
  # the tumour options of the run are kept in params (call_scn no longer replaces the caller's record)
  expect_identical(res$params$tumour, clone_defaults())
  expect_null(res$params$normal)
  out = file.path(tempdir(), "callscn_tum")
  write_scphasedcn(res, out)
  expect_true(all(file.exists(file.path(out, c("clones.tsv", "assignments.tsv", "clone_segments.tsv", "clone_profiles.tsv.gz",
                                               "reference_pool.txt", "cells.tsv", "calls.tsv", "result.rds")))))
  expect_equal(caller_params(file.path(out, "params.yaml")), res$params)
})

test_that("inst/scripts/call_scn.R --help exits 0 without data", {
  script = system.file("scripts", "call_scn.R", package = "scphasedcn")
  expect_true(file.exists(script))
  rs = file.path(R.home("bin"), "Rscript")
  out = suppressWarnings(system2(rs, c(script, "--help"), stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"))
  expect_true(any(grepl("--sample", out)))
  expect_true(any(grepl("--null-cells", out)))
})

test_that("call_scn.R runs on --dir / --out alone, without a sample yaml", {
  # the script loads the installed package: run it under R CMD check (or SCPHASEDCN_TEST_CLI=true after an install)
  skip_if_not(testthat::is_checking() || isTRUE(as.logical(Sys.getenv("SCPHASEDCN_TEST_CLI", "FALSE"))))
  script = system.file("scripts", "call_scn.R", package = "scphasedcn")
  rs = file.path(R.home("bin"), "Rscript")
  d = sim_dir("sim_norm_f")
  cen = file.path(tempdir(), "cli_centromeres.tsv")
  fwrite(cs_centromeres(sim_norm_f$scn$bins), cen, sep = "\t")
  out = file.path(tempdir(), "callscn_cli")
  unlink(out, recursive = TRUE)
  log = suppressWarnings(system2(rs, c(script, "--dir", d, "--out", out, "--sex", "female", "--centromeres", cen,
                                       "--cores", "1"), stdout = TRUE, stderr = TRUE))
  expect_null(attr(log, "status"))
  expect_true(all(file.exists(file.path(out, c("cells.tsv", "calls.tsv", "params.yaml", "NOTES.md")))))
  # the same calls as the R API
  api = file.path(tempdir(), "callscn_cli_api")
  write_scphasedcn(norm_result(), api)
  for (f in c("calls.tsv", "cells.tsv", "thresholds.tsv"))
    expect_identical(readLines(file.path(out, f)), readLines(file.path(api, f)), info = f)
  # no --out and no sample yaml: a clear error
  log2 = suppressWarnings(system2(rs, c(script, "--dir", d, "--sex", "female"), stdout = TRUE, stderr = TRUE))
  expect_false(is.null(attr(log2, "status")))
  expect_true(any(grepl("give --out", log2)))
  # without --config the sample yaml comes from SCPHASEDCN_SAMPLE_CONFIG
  missing_yaml = file.path(tempdir(), "no_such_sample.yaml")
  log3 = suppressWarnings(system2(rs, c(script, "--sample", "s1"), stdout = TRUE, stderr = TRUE,
                                  env = paste0("SCPHASEDCN_SAMPLE_CONFIG=", missing_yaml)))
  expect_false(is.null(attr(log3, "status")))
  expect_true(any(grepl("sample yaml not found: .*no_such_sample.yaml", log3)))
})

test_that("make_cell_meta builds the cell metadata from cells_qc and user labels", {
  qc = data.table(cell = c("AAAC-1", "AAAG", "AAAT"), n_clean = c(40000, 25000, 900))
  types = data.table(cell = c("AAAC", "AAAG", "CCCC"), cell_type = c("epi", "immune", "epi"))
  expect_message(m <- make_cell_meta(qc, cell_types = types, flat = c("AAAG"), groups = c(AAAC = "g1")), "not in cells_qc")
  expect_equal(m$cell, c("AAAC", "AAAG", "AAAT"))
  expect_equal(m$qc_pass, c(TRUE, TRUE, FALSE))
  expect_equal(m$cell_type, c("epi", "immune", NA))
  expect_equal(m$flat, c(FALSE, TRUE, FALSE))
  expect_equal(m$group, c("g1", NA, NA))
  # no labels: NA columns; an extraction directory as input
  d = sim_dir("sim_norm_f")
  m2 = make_cell_meta(d)
  expect_equal(nrow(m2), nrow(fread(file.path(d, "cells_qc.tsv"))))
  expect_true(all(is.na(m2$cell_type)) && all(is.na(m2$flat)))
  expect_error(make_cell_meta(), "no cells")
})

# ---- one-copy model of the male sex chromosomes (params$haploid_model) ---------------------------------------------

male_result = function(haploid) {
  key = paste0("male_", haploid)
  if (is.null(call_cache[[key]])) {
    p = modifyList(caller_defaults(), list(haploid_model = haploid))
    call_cache[[key]] = call_scn(sim_dir("sim_norm_m"), "normal", sex = "male", centromeres = cs_centromeres(sim_norm_m$scn$bins),
                                 params = p, verbose = FALSE)
  }
  call_cache[[key]]
}
# singlet carriers of sim_norm_m's private events (1 = loss of Y, 2 = chrX gain, 3 = chr8 gain) and all carrier barcodes
male_carriers = function() {
  tr = sim_norm_m$truth
  cells = sim_norm_m$scn$cells$cell
  car = lapply(tr$carriers, function(n) unique(cells[tr$nuclei$cell[n]]))
  single = tr$cells[n_nuclei == 1, cell]
  list(all = car, single = lapply(car, intersect, single))
}

test_that("simulate_scn: a male has one chrX and one chrY copy, a female no chrY", {
  g = sim_genome(c(chr1 = 20e6, chrX = 10e6, chrY = 10e6))
  m = simulate_scn(n_cells = 20, genome = g, sex = "male", events = NULL, multiplet_rate = 0, seed = 3)
  y = m$scn$bins$chrom == "chrY"; x = m$scn$bins$chrom == "chrX"
  expect_true(all(m$truth$cn[y, ] == 1) && all(m$truth$cn[x, ] == 1) && all(m$truth$cn[!x & !y, ] == 2))
  expect_gt(sum(m$scn$counts[y, ]), 0)
  f = simulate_scn(n_cells = 20, genome = g, sex = "female", events = NULL, multiplet_rate = 0, seed = 3)
  expect_true(all(f$truth$cn[y, ] == 0) && all(f$truth$cn[!y, ] == 2))
  expect_equal(sum(f$scn$counts[y, ]), 0)
  # loss of Y in sim_norm_m: carriers have no chrY copy and no chrY fragments
  car = male_carriers()$single[[1]]
  expect_gte(length(car), 3)
  ym = sim_norm_m$scn$bins$chrom == "chrY"
  expect_true(all(sim_norm_m$truth$cn[ym, car] == 0))
  expect_equal(sum(sim_norm_m$scn$counts[ym, car]), 0)
})

test_that("haploid_model = FALSE: no chrY regions, no ploidy, chrY bins unused (the model without the option)", {
  r0 = male_result(FALSE)
  expect_false(any(r0$regions$chrom == "chrY"))
  expect_false("ploidy" %in% names(r0$regions))
  expect_false("ploidy" %in% names(r0$features))
  expect_true(all(is.na(r0$thresholds[level == "cell", sexgrp])))
  expect_false(any(grepl("one-copy model", r0$notes)))
})

test_that("haploid_model = TRUE: loss of Y and male X gains are called per cell; autosomes, pools, female and tumour unchanged", {
  r0 = male_result(FALSE); r1 = male_result(TRUE)
  car = male_carriers()
  # chrY regions (no centromere row is needed: the test table has one, so chrYp / chrYq exist here), all one copy
  yr = r1$regions[chrom == "chrY"]
  expect_true(all(c("chrY", "5Mb|chrY:0") %in% yr$region))
  expect_true(all(r1$regions[chrom %in% c("chrX", "chrY"), ploidy] == 1L))
  expect_true(all(r1$regions[!chrom %in% c("chrX", "chrY"), ploidy] == 2L))
  expect_false(any(yr$null_region) || any(yr$allelic_ok))
  expect_true(any(grepl("one-copy model", r1$notes)))
  # per-cell calls: loss of Y on its carriers, chrX gains on theirs; no chrX / chrY call on a non-carrier
  cl = r1$calls[level == "cell" & called & reportable]
  gy = unique(cl[region == "chrY" & state == "loss", cell])
  gx = unique(cl[region == "chrX" & state == "gain", cell])
  rec_y = length(intersect(gy, car$single[[1]])) / length(car$single[[1]])
  rec_x = length(intersect(gx, car$single[[2]])) / length(car$single[[2]])
  non_y = setdiff(r1$cells$cell, car$all[[1]]); non_x = setdiff(r1$cells$cell, car$all[[2]])
  fr_y = length(intersect(gy, non_y)) / length(non_y); fr_x = length(intersect(gx, non_x)) / length(non_x)
  message(sprintf("haploid_model on sim_norm_m: loss of Y recall %.2f (%d singlet carriers), false %.3f; chrX gain recall %.2f (%d), false %.3f",
                  rec_y, length(car$single[[1]]), fr_y, rec_x, length(car$single[[2]]), fr_x))
  expect_gte(rec_y, 0.8); expect_gte(rec_x, 0.8)
  expect_lte(fr_y, 0.01); expect_lte(fr_x, 0.01)
  expect_equal(cl[chrom %in% c("chrX", "chrY") & !cell %in% c(car$all[[1]], car$all[[2]]), .N], 0L)
  # without the option loss of Y is not callable (no chrY rows)
  expect_equal(r0$calls[chrom == "chrY", .N], 0L)
  # one-copy rows have thresholds of their own; their realised false-call rate on the null cells stays small
  th = r1$thresholds[level == "cell"]
  expect_setequal(unique(th$sexgrp), c("A", "X", "Y"))
  sc = r1$scores[ploidy == 1 & null_cell == TRUE & !(chrom == "chrY" & cell %in% car$all[[1]]) & !(chrom == "chrX" & cell %in% car$all[[2]])]
  expect_lte(mean(sc$call_s_comb), 0.01)
  # autosomes: features, scores, thresholds, cell and pooled calls identical to the run without the option
  auto = function(x) x[!chrom %in% c("chrX", "chrY")]
  key = c("cell", "region")
  f0 = auto(r0$features); f1 = auto(r1$features)[, names(f0), with = FALSE]
  setkeyv(f0, key); setkeyv(f1, key)
  expect_identical(f1, f0)
  s0 = auto(r0$scores); s1 = auto(r1$scores)[, names(s0), with = FALSE]
  setkeyv(s0, key); setkeyv(s1, key)
  expect_identical(s1, s0)
  t0 = r0$thresholds; t1 = r1$thresholds[sexgrp == "A" | level != "cell"]
  cn = setdiff(names(t0), "sexgrp")
  expect_identical(t1[level == "cell", cn, with = FALSE], t0[level == "cell", cn, with = FALSE])
  expect_identical(t1[level != "cell", names(t0), with = FALSE], t0[level != "cell"])
  c0 = r0$calls[!chrom %in% c("chrX", "chrY") | level == "pool"]; c1 = r1$calls[!chrom %in% c("chrX", "chrY") | level == "pool"]
  expect_identical(c1, c0)                                                 # incl. the chrX pooled calls (two-copy model)
  expect_identical(r1$embedding, r0$embedding)
  expect_identical(r1$pools, r0$pools)
  ev = function(x) x[!chrom %in% c("chrX", "chrY")]
  expect_identical(ev(r1$events), ev(r0$events))
  # a female run is unaffected by the option
  rf = call_scn(sim_dir("sim_norm_f"), "normal", sex = "female", centromeres = cs_centromeres(sim_norm_f$scn$bins), verbose = FALSE,
                params = modifyList(caller_defaults(), list(haploid_model = TRUE)))
  r_f = norm_result()
  for (k in c("calls", "thresholds", "events", "regions", "features", "scores", "pools", "embedding")) expect_identical(rf[[k]], r_f[[k]], info = k)
  expect_true(any(grepl("haploid_model: not used \\(sex female\\)", rf$feature_info$notes)))
  # tumour mode ignores the option with a warning and gives the same result
  r_t = call_cache$tum %||% tum_result()
  expect_warning(r_t1 <- tum_result(modifyList(caller_defaults(), list(haploid_model = TRUE))), "ignored in tumour mode")
  for (k in c("calls", "clones", "assignments", "segments", "regions")) expect_identical(r_t1[[k]], r_t[[k]], info = k)
  expect_false(r_t1$params$haploid_model)
})
