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

test_that("call_scn with mask_f0 runs the bias mask pass and records it; depth is untouched", {
  res = norm_result()
  expect_null(res$null_fit$allelic_mask)                                  # default: no mask pass
  r = call_scn(sim_dir("sim_norm_f"), "normal", sex = "female", centromeres = cs_centromeres(sim_norm_f$scn$bins), verbose = FALSE,
               params = modifyList(caller_defaults(), list(mask_f0 = TRUE)))
  am = r$null_fit$allelic_mask
  expect_true(all(c("f0", "d0", "z_f0", "masked_rho", "masked_f0") %in% names(am)))
  expect_identical(am$masked, am$masked_f0)                               # mask_allelic off: the bias rule alone
  expect_true(all(am$chrom %in% AUTOSOMES_CALLER))
  expect_true(any(grepl("allelic mask: .*bias rule .*: \\d+", r$notes)))
  # sim_norm_f has no shared allelic bias: nothing is masked at the default thresholds and the calls are those of
  # the default run
  expect_equal(sum(am$masked), 0)
  expect_identical(r$features$z, res$features$z)
  expect_identical(r$calls, res$calls)
  # both rules on, the bias rule forced to fire (thresholds below the simulation's noise): the union is masked,
  # the run notes give the count per rule, the allele votes of the masked tiles are gone, depth is untouched and
  # the truth checks of the default run still hold
  prm = modifyList(caller_defaults(), list(mask_allelic = TRUE, mask_f0 = TRUE, mask_f0_mad = 1.5, mask_f0_min = 0.005))
  r2 = call_scn(sim_dir("sim_norm_f"), "normal", sex = "female", centromeres = cs_centromeres(sim_norm_f$scn$bins), verbose = FALSE,
                params = prm)
  am2 = r2$null_fit$allelic_mask
  expect_gt(sum(am2$masked_f0), 0)
  expect_identical(am2$masked, am2$masked_rho | am2$masked_f0)
  nt = regmatches(r2$notes, regexpr("overdispersion rule [^;]*: \\d+; bias rule [^;]*: \\d+", r2$notes))
  expect_length(nt, 1)
  expect_equal(as.integer(regmatches(nt, gregexpr("\\d+(?=$|;)", nt, perl = TRUE))[[1]]), c(sum(am2$masked_rho), sum(am2$masked_f0)))
  mt = am2[masked == TRUE, tile]
  expect_true(all(r2$features[region %in% mt, n_hap] == 0))
  expect_gt(res$features[region %in% mt, sum(n_hap)], 0)
  expect_identical(r2$features$z, res$features$z)
  tr = sim_norm_f$truth
  cells = sim_norm_f$scn$cells$cell
  car = lapply(tr$carriers, function(n) unique(cells[tr$nuclei$cell[n]]))
  cl = r2$calls[level == "cell" & called & reportable]
  for (k in list(list("chr1p", "loss", 1), list("chr8", "gain", 2), list("chr21", "gain", 4))) {
    got = unique(cl[region == k[[1]] & state == k[[2]], cell])
    expect_true(all(got %in% car[[k[[3]]]]), info = k[[1]])
    expect_gte(length(intersect(got, car[[k[[3]]]])) / length(car[[k[[3]]]]), 0.7)
  }
  # the reportable cell calls on chromosomes without a masked tile are those of the default run (here 5 tiles of
  # chr1, chr2 and chr22 are masked; the 56 calls elsewhere are unchanged)
  key = function(x) x[level == "cell" & called & reportable & !chrom %in% am2[masked == TRUE, chrom], paste(cell, region, state)]
  expect_gt(length(key(res$calls)), 0)
  expect_setequal(key(r2$calls), key(res$calls))
})

test_that("call_scn does not apply the bias rule (mask_f0) under mirrored orientation: noted, warned, results unchanged", {
  cen = cs_centromeres(sim_norm_f$scn$bins)
  f = suppressWarnings(cell_features(sim_dir("sim_norm_f"), sex = "female", centromeres = cen, phasing = "mirrored"))
  expect_identical(features_orientation(f), "mirrored")
  # thresholds at zero: under mirrored orientation about half the tiles would be bias outliers if the rule ran
  prm = modifyList(caller_defaults(), list(mask_f0 = TRUE, mask_f0_mad = 0, mask_f0_min = 0))
  run = function(p) call_scn(sim_dir("sim_norm_f"), "normal", sex = "female", features = f, params = p, embedding = FALSE,
                             verbose = FALSE)
  expect_warning(r1 <- run(prm), "mask_f0.*mirrored")
  r0 = run(caller_defaults())
  expect_null(r1$null_fit$allelic_mask)                                   # no tile masked, no mask pass
  expect_equal(sum(grepl("bias rule \\(mask_f0\\) not applied: under mirrored orientation", r1$notes)), 1)
  expect_false(any(grepl("mask_f0", r0$notes)))
  for (k in c("features", "scores", "thresholds", "calls", "events", "pools", "null_fit", "false_calls", "false_events"))
    expect_identical(r1[[k]], r0[[k]], info = k)
  # the overdispersion rule still runs under mirrored orientation; the bias rule is not part of it
  w = character()
  r2 = withCallingHandlers(run(modifyList(prm, list(mask_allelic = TRUE))),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_equal(sum(grepl("mask_f0", w)), 1)                                # warned once
  expect_identical(r2$null_fit$allelic_mask$masked, r2$null_fit$allelic_mask$masked_rho)
  # the orientation of a features object without `phasing`: from its long table
  expect_identical(features_orientation(list(long = f$long)), "mirrored")
  l2 = copy(f$long)[allelic_ok %in% TRUE, orient := "pseudobulk"]
  expect_true(is.na(features_orientation(list(long = l2))))
})

test_that("scn_embedding is reproducible (seeded truncated SVD) and named by cell", {
  res = norm_result()
  d = sim_dir("sim_norm_f")
  e1 = scn_embedding(d, res$cells$cell, res$feature_info$use)
  expect_equal(rownames(e1), res$cells$cell)
  expect_identical(e1, scn_embedding(d, res$cells$cell, res$feature_info$use))
  expect_identical(e1, res$embedding)
})

test_that("tumour mode dispatches to the tumour mode", {
  d = sim_dir("sim_tum2")
  s = sim_tum2
  cm = data.table(cell = s$scn$cells$cell, group = s$scn$cells$cell_type)
  res = call_scn(d, "tumour", cell_meta = cm, cells = cm$cell, sex = "male", groups = "group", verbose = FALSE,
                 centromeres = cs_centromeres(s$scn$bins))
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
