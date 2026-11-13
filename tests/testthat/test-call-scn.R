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
# a simulation written with the phase sets of `chroms` replaced: `block` bp read-backed blocks (orientation kept), or
# one phase set per chromosome (block = NULL), e.g. a hybrid phase with one chromosome left in blocks
rephased_dir = function(name, chroms, block = NULL) {
  key = paste0(name, "_", paste(chroms, collapse = "_"), "_", if (is.null(block)) "chrom" else format(block), "_dir")
  if (is.null(call_cache[[key]])) {
    s = get(name)
    s$fragments = copy(s$fragments)
    if (is.null(block)) s$fragments[chrom %in% chroms & ps > 0, ps := min(ps), by = chrom]
    else s$fragments[chrom %in% chroms & ps > 0, ps := as.integer(floor(start / block) * block + 1)]
    d = file.path(tempdir(), paste0("callscn_", sub("_dir$", "", key)))
    unlink(d, recursive = TRUE)
    invisible(capture.output(suppressMessages(write_sim_dataset(s, d))))
    call_cache[[key]] = d
  }
  call_cache[[key]]
}
# call_scn with its warnings collected: list(r, w)
call_w = function(...) {
  w = character()
  r = withCallingHandlers(call_scn(...), warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  list(r = r, w = w)
}
norm_result = function() {
  if (is.null(call_cache$norm)) {
    t0 = Sys.time()
    d = sim_dir("sim_norm_f")
    w = character()
    call_cache$norm = withCallingHandlers(
      call_scn(d, "normal", sex = "female", centromeres = cs_centromeres(sim_norm_f$scn$bins), verbose = FALSE),
      warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
    call_cache$norm_warnings = w
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

test_that("the gain gate arguments of caller_normal() round-trip through caller.yaml and params.yaml", {
  pf = caller_params()
  expect_identical(pf$normal$gain_gate, "dl")
  expect_identical(pf$normal$gain_gate_q, 0.95)
  tf = tempfile(fileext = ".yaml")
  writeLines(c("normal:", "  gain_gate: \"z\"", "  gain_gate_q: 0.99"), tf)
  p = caller_params(tf)
  expect_identical(p$normal$gain_gate, "z")
  expect_identical(p$normal$gain_gate_q, 0.99)
  expect_identical(p$normal$gain_min_dl, 0)
  y = parse_simple_yaml(readLines(tf))                          # the fallback parser reads the block the same way
  expect_equal(block_params(y$normal, normal_arg_defaults(), "normal"), p$normal)
  tf2 = tempfile(fileext = ".yaml")                             # params.yaml of a run reads back
  writeLines(yaml_lines(p), tf2)
  expect_equal(caller_params(tf2)$normal, p$normal)
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
  # every normal-mode run on mirrored features also warns that this orientation is a diagnostic; collect the warnings
  run = function(p) call_w(sim_dir("sim_norm_f"), "normal", features = f, params = p, embedding = FALSE, verbose = FALSE)
  x1 = run(prm); r1 = x1$r
  expect_length(grep("mask_f0.*mirrored", x1$w), 1L)
  x0 = run(caller_defaults()); r0 = x0$r
  expect_false(any(grepl("mask_f0", x0$w)))
  expect_null(r1$null_fit$allelic_mask)                                   # no tile masked, no mask pass
  expect_equal(sum(grepl("bias rule \\(mask_f0\\) not applied: under mirrored orientation", r1$notes)), 1)
  expect_false(any(grepl("mask_f0", r0$notes)))
  for (k in c("features", "scores", "thresholds", "calls", "events", "pools", "null_fit", "false_calls", "false_events"))
    expect_identical(r1[[k]], r0[[k]], info = k)
  # the overdispersion rule still runs under mirrored orientation; the bias rule is not part of it
  x2 = run(modifyList(prm, list(mask_allelic = TRUE))); r2 = x2$r
  expect_equal(sum(grepl("mask_f0", x2$w)), 1)                             # warned once
  expect_identical(r2$null_fit$allelic_mask$masked, r2$null_fit$allelic_mask$masked_rho)
  # the orientation of a features object without `phasing`: from its long table
  expect_identical(features_orientation(list(long = f$long)), "mirrored")
  l2 = copy(f$long)[allelic_ok %in% TRUE, orient := "pseudobulk"]
  expect_true(is.na(features_orientation(list(long = l2))))
})

test_that("call_scn with both allelic masks on a run without any allele votes: nothing masked, depth-only with one warning", {
  # an unphased extraction: no molecule vote and nA = nB = n_hap = 0 on every row; caller_normal warns and runs
  # depth-only, and the mask pass (mask_allelic, mask_f0) has no tile to assess
  f = suppressWarnings(cell_features(sim_dir("sim_norm_f"), sex = "female", centromeres = cs_centromeres(sim_norm_f$scn$bins)))
  f$long[, `:=`(nA = 0, nB = 0, n_hap = 0)]
  if ("f" %in% names(f$long)) f$long[, f := NA_real_]
  f$votes = f$votes[0]
  run = function(p) {
    ws = character()
    r = withCallingHandlers(call_scn(sim_dir("sim_norm_f"), "normal", features = f, params = p, embedding = FALSE, verbose = FALSE),
                            warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(r = r, ws = ws)
  }
  x = run(modifyList(caller_defaults(), list(mask_allelic = TRUE, mask_f0 = TRUE)))
  x0 = run(caller_defaults())
  expect_length(grep("no allele votes", x$ws), 1)
  expect_false(any(grepl("mask", x$ws)))
  am = x$r$null_fit$allelic_mask
  expect_gt(nrow(am), 0)
  expect_false(any(am$assessed))
  expect_false(any(am$masked))
  expect_true(any(grepl("allelic mask: 0 of 0 assessed 5 Mb tiles masked", x$r$notes)))
  expect_true(any(grepl("no allele votes", x$r$notes)))
  expect_gt(nrow(x$r$calls), 0)
  expect_identical(x$r$calls, x0$r$calls)
})

test_that("normal mode warns on block-only phasing (loo orientation), not on chromosome-wide phasing", {
  # sim_norm_f: one phase set per chromosome -> chromosome orientation, no warning, no run note
  res = norm_result()
  expect_equal(unname(res$feature_info$phasing), c("chromosome", "chromosome"))
  expect_false(any(grepl("block phasing", call_cache$norm_warnings)))
  expect_false(any(grepl("block phasing", res$notes)))
  pc0 = res$feature_info$phasing_chrom
  expect_s3_class(pc0, "data.table")
  expect_setequal(pc0$chrom, unique(sim_norm_f$scn$bins$chrom))
  expect_equal(unique(pc0$regime), "chromosome")
  expect_false(any(grepl("while the run's regime", res$notes, fixed = TRUE)))
  # sim_norm_blk: unoriented 200 kb blocks -> loo orientation -> one warning and a run note that reaches NOTES.md
  s = sim_norm_blk
  d = sim_dir("sim_norm_blk")
  w = character()
  r = withCallingHandlers(
    call_scn(d, "normal", sex = "female", centromeres = cs_centromeres(s$scn$bins), embedding = FALSE, verbose = FALSE,
             keep_scores = FALSE),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  hit = grepl("normal mode on block phasing", w)
  expect_equal(sum(hit), 1L)
  expect_match(w[hit], "scphasedcn-inputs", fixed = TRUE)
  expect_match(w[hit], "scphasedcn hybrid-phase", fixed = TRUE)
  expect_equal(unname(r$feature_info$phasing), c("block", "loo"))
  pc = r$feature_info$phasing_chrom
  expect_s3_class(pc, "data.table")
  expect_equal(anyDuplicated(pc$chrom), 0L)
  expect_setequal(pc$chrom, unique(s$scn$bins$chrom))
  expect_true(all(pc$n_votes > 0))
  expect_equal(unique(pc$regime), "block")
  expect_true(any(grepl("phasing per chromosome", r$feature_info$notes)))
  expect_equal(sum(grepl("normal mode on block phasing", r$notes)), 1L)
  out = file.path(tempdir(), "callscn_blk")
  unlink(out, recursive = TRUE)
  write_scphasedcn(r, out)
  expect_true(any(grepl("normal mode on block phasing", readLines(file.path(out, "NOTES.md")))))
  expect_equal(read_scphasedcn(out)$feature_info$phasing_chrom, pc)
  # the note follows the features' phasing record: precomputed block features warn too, tumour mode never
  expect_null(phasing_run_note(c(detected = "block", orientation = "loo"), "tumour"))
  expect_null(phasing_run_note(c(detected = "chromosome", orientation = "chromosome"), "normal"))
  expect_null(phasing_run_note(c(detected = "block", orientation = "mirrored"), "normal"))
  expect_null(phasing_run_note(NULL, "normal"))
  expect_type(phasing_run_note(c(detected = "block", orientation = "block"), "normal"), "character")
  # `phasing` is checked up front
  expect_error(call_scn(d, "normal", phasing = "bogus"), "should be one of")
})

test_that("chromosomes left in blocks under chromosome-scale phasing reach the run notes and NOTES.md", {
  s = sim_norm_f
  x = call_w(rephased_dir("sim_norm_f", "chr22", block = 2e5), "normal", sex = "female", centromeres = cs_centromeres(s$scn$bins),
             embedding = FALSE, verbose = FALSE, keep_scores = FALSE)
  call_cache$mixed = x
  r = x$r
  expect_equal(unname(r$feature_info$phasing), c("chromosome", "chromosome"))
  pc = r$feature_info$phasing_chrom
  expect_equal(pc[regime == "block", chrom], "chr22")
  expect_setequal(pc[regime == "chromosome", chrom], setdiff(unique(s$scn$bins$chrom), "chr22"))
  msg = "chr22 phased in blocks while the run's regime is chromosome"
  expect_equal(sum(grepl(msg, r$notes, fixed = TRUE)), 1L)
  # the note says what chromosome orientation does to the blocks (summed unoriented); not under loo
  expect_match(r$notes[grepl(msg, r$notes, fixed = TRUE)], "attenuates or cancels their allele signal", fixed = TRUE)
  pcb = data.table(chrom = c("chr1", "chr2"), n_votes = c(100, 80), n_ps = c(1L, 40L), k = c(1L, 30L),
                   regime = c("chromosome", "block"))
  expect_false(grepl("attenuates", phasing_mixed_note(pcb, "block", "loo"), fixed = TRUE))
  expect_match(phasing_mixed_note(pcb, "block", "loo"), "chr1 phased chromosome-wide", fixed = TRUE)
  out = file.path(tempdir(), "callscn_mixed")
  unlink(out, recursive = TRUE)
  write_scphasedcn(r, out)
  expect_equal(sum(grepl(msg, readLines(file.path(out, "NOTES.md")), fixed = TRUE)), 1L)
})

test_that("chromosome orientation on block phasing and mirrored orientation in normal mode warn", {
  s = sim_norm_blk
  d = sim_dir("sim_norm_blk")
  run = function(ph) call_w(d, "normal", sex = "female", centromeres = cs_centromeres(s$scn$bins), embedding = FALSE,
                            verbose = FALSE, keep_scores = FALSE, phasing = ph)
  forced = 'phasing = "chromosome" on block phasing'
  mirr = 'normal mode with phasing = "mirrored"'
  # `chromosome` forced on unoriented blocks: one warning, one run note in NOTES.md
  x = run("chromosome")
  expect_equal(unname(x$r$feature_info$phasing), c("block", "chromosome"))
  hit = grepl(forced, x$w, fixed = TRUE)
  expect_equal(sum(hit), 1L)
  expect_match(x$w[hit], "carries no signal", fixed = TRUE)
  expect_false(any(grepl(mirr, x$w, fixed = TRUE)))
  out = file.path(tempdir(), "callscn_forced")
  unlink(out, recursive = TRUE)
  write_scphasedcn(x$r, out)
  expect_equal(sum(grepl(forced, readLines(file.path(out, "NOTES.md")), fixed = TRUE)), 1L)
  # `mirrored`: one warning, one run note
  y = run("mirrored")
  hit = grepl(mirr, y$w, fixed = TRUE)
  expect_equal(sum(hit), 1L)
  expect_match(y$w[hit], "not calibrated for calling", fixed = TRUE)
  expect_false(any(grepl(forced, y$w, fixed = TRUE)))
  expect_equal(sum(grepl(mirr, y$r$notes, fixed = TRUE)), 1L)
  # valid cases are silent: auto on chromosome-wide phasing (sim_norm_f) and on blocks (loo note only), chromosome
  # requested on chromosome-wide phasing
  res = norm_result()
  expect_false(any(grepl(paste0(forced, "|", mirr), call_cache$norm_warnings)))
  z = run("auto")
  expect_false(any(grepl(forced, z$w, fixed = TRUE) | grepl(mirr, z$w, fixed = TRUE)))
  expect_length(phasing_orientation_notes(res$feature_info$phasing, res$feature_info$phasing_chrom, "normal", "chromosome"), 0L)
  # one chromosome in blocks: `chromosome` requested names it, `auto` does not warn
  pc = data.table(chrom = c("chr1", "chr2", "chrX"), n_votes = c(100, 80, 40), n_ps = c(1L, 2L, 30L), k = c(1L, 1L, 20L),
                  regime = c("chromosome", "chromosome", "block"))
  ph = c(detected = "chromosome", orientation = "chromosome")
  n1 = phasing_orientation_notes(ph, pc, "normal", "chromosome")
  expect_length(n1, 1L)
  expect_match(n1, paste0(forced, " (chrX)"), fixed = TRUE)
  expect_length(phasing_orientation_notes(ph, pc, "normal", "auto"), 0L)
  expect_length(phasing_orientation_notes(ph, pc[regime == "chromosome"], "normal", "chromosome"), 0L)
  # tumour mode: the forced warning applies, mirrored is not flagged
  expect_length(phasing_orientation_notes(c(detected = "block", orientation = "chromosome"), NULL, "tumour", "chromosome"), 1L)
  expect_length(phasing_orientation_notes(c(detected = "block", orientation = "mirrored"), NULL, "tumour", "mirrored"), 0L)
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
  # multiplet scores between 0.5 and the tumour mode's multiplet_prob (0.9) are not flagged: the cells are used as singlets
  mult = data.table(cell = cm$cell[1:4], prob = c(0.6, 0.85, 0.95, 1))
  x = call_w(d, "tumour", cell_meta = cm, cells = cm$cell, sex = "male", groups = "group", verbose = FALSE,
             centromeres = cs_centromeres(s$scn$bins), multiplets = mult)
  res = x$r
  expect_false(any(grepl("chromosome-scale phasing", x$w, fixed = TRUE)))
  expect_false(any(grepl("chromosome-scale phasing", res$notes, fixed = TRUE)))
  expect_s3_class(res, "scphasedcn")
  expect_equal(res$mode, "tumour")
  expect_true(all(c("clone", "posterior") %in% names(res$cells)))
  expect_identical(res$cells[match(mult$cell, cell), multiplet], c(FALSE, FALSE, TRUE, TRUE))
  asg = as.data.table(res$assignments)
  expect_identical(res$cells[match(asg$cell, cell), multiplet], asg$multiplet)
  # the tumour options of the run are kept in params (call_scn no longer replaces the caller's record)
  expect_identical(res$params$tumour, clone_defaults())
  expect_null(res$params$normal)
  out = file.path(tempdir(), "callscn_tum")
  write_scphasedcn(res, out)
  expect_true(all(file.exists(file.path(out, c("clones.tsv", "assignments.tsv", "clone_segments.tsv", "clone_profiles.tsv.gz",
                                               "reference_pool.txt", "cells.tsv", "calls.tsv", "result.rds")))))
  expect_equal(caller_params(file.path(out, "params.yaml")), res$params)
})

test_that("tumour mode on chromosome-scale phasing warns once and notes it", {
  s = sim_tum2
  d = rephased_dir("sim_tum2", unique(s$scn$bins$chrom))
  cm = data.table(cell = s$scn$cells$cell, group = s$scn$cells$cell_type)
  x = call_w(d, "tumour", cell_meta = cm, cells = cm$cell, sex = "male", groups = "group", verbose = FALSE,
             centromeres = cs_centromeres(s$scn$bins))
  expect_equal(x$r$feature_info$phasing[["detected"]], "chromosome")
  msg = "tumour mode on chromosome-scale phasing"
  hit = grepl(msg, x$w, fixed = TRUE)
  expect_equal(sum(hit), 1L)
  expect_match(x$w[hit], "cannot be re-oriented per clone", fixed = TRUE)
  out = file.path(tempdir(), "callscn_tum_chrom")
  unlink(out, recursive = TRUE)
  write_scphasedcn(x$r, out)
  expect_equal(sum(grepl(msg, readLines(file.path(out, "NOTES.md")), fixed = TRUE)), 1L)
  expect_type(phasing_run_note(c(detected = "chromosome", orientation = "chromosome"), "tumour"), "character")
  expect_null(phasing_run_note(c(detected = "chromosome", orientation = "chromosome"), "tumour", phasing_chrom = data.table(chrom = character(), regime = character())))
})

# a simulation written without phase sets (ps = 0: an unphased extraction), so the run has no allele votes at all
unphased_dir = function(name) {
  key = paste0(name, "_unphased_dir")
  if (is.null(call_cache[[key]])) {
    s = get(name)
    s$fragments = copy(s$fragments)[, ps := 0L]
    d = file.path(tempdir(), paste0("callscn_", name, "_unphased"))
    unlink(d, recursive = TRUE)
    invisible(capture.output(suppressMessages(write_sim_dataset(s, d))))
    call_cache[[key]] = d
  }
  call_cache[[key]]
}

test_that("a run without allele votes: the phasing checks stay quiet (normal: depth-only calls; tumour: the caller's stop)", {
  s = sim_norm_f
  # phasing = "chromosome" requested: the check of block-phased chromosomes sees an empty per-chromosome table
  x = call_w(unphased_dir("sim_norm_f"), "normal", sex = "female", centromeres = cs_centromeres(s$scn$bins),
             phasing = "chromosome", verbose = FALSE)
  expect_s3_class(x$r, "scphasedcn")
  expect_equal(nrow(x$r$feature_info$phasing_chrom), 0L)
  expect_equal(unname(x$r$feature_info$phasing[["orientation"]]), "chromosome")
  expect_length(grep("no allele votes", x$w), 1L)
  expect_false(any(grepl("block phasing|chromosome-scale phasing|phased in blocks|chromosome-wide", x$w)))
  expect_false(any(grepl("block phasing|chromosome-scale phasing|phased in blocks|chromosome-wide", x$r$notes)))
  out = file.path(tempdir(), "callscn_norm_unphased")
  unlink(out, recursive = TRUE)
  expect_no_error(write_scphasedcn(x$r, out))
  # the other requested regimes on the same (empty) phasing record give no phasing note either
  pc = x$r$feature_info$phasing_chrom
  for (ori in c("loo", "chromosome")) {
    ph = c(detected = "chromosome", orientation = ori)
    expect_null(phasing_run_note(ph, "normal", pc))
    expect_length(phasing_orientation_notes(ph, pc, "normal", if (ori == "loo") "block" else "chromosome"), 0L)
    expect_null(phasing_mixed_note(pc, "chromosome", ori))
  }
  # tumour mode stops in the caller with its own message, not in the phasing code
  t2 = sim_tum2
  cm = data.table(cell = t2$scn$cells$cell, group = t2$scn$cells$cell_type)
  expect_error(suppressWarnings(call_scn(unphased_dir("sim_tum2"), "tumour", cell_meta = cm, cells = cm$cell, sex = "male",
                                         groups = "group", verbose = FALSE, centromeres = cs_centromeres(t2$scn$bins))),
               "caller_tumour: no allele votes")
})

test_that("inst/scripts/call_scn.R --help exits 0 without data", {
  script = system.file("scripts", "call_scn.R", package = "scphasedcn")
  expect_true(file.exists(script))
  rs = file.path(R.home("bin"), "Rscript")
  out = suppressWarnings(system2(rs, c(script, "--help"), stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"))
  expect_true(any(grepl("--sample", out)))
  expect_true(any(grepl("--null-cells", out)))
  expect_true(any(grepl("--phasing", out)))
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

test_that("cell selection: the cell QC applies, a cell_meta qc_pass column restricts it further; `cells` bypasses QC", {
  all = paste0("c", 1:6)
  qct = data.table(cell = all, qc_pass = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))
  expect_identical(feature_cells(all, qct)$keep, c("c1", "c2", "c3", "c6"))
  # a cell_meta without qc_pass: its cells that pass the QC (the CLI no longer passes them as `cells`)
  m = data.table(cell = c("c6", "c4", "c2", "x9"))
  s = feature_cells(all, qct, cell_meta = m)
  expect_identical(s$keep, c("c2", "c6"))
  expect_identical(s$how, "QC pass")
  # with qc_pass: a cell must pass both the QC and the table's qc_pass (c4 fails the QC, c2 / c1 the table)
  m2 = data.table(cell = c("c6", "c4", "c2", "c1"), qc_pass = c(TRUE, TRUE, FALSE, NA))
  s2 = feature_cells(all, qct, cell_meta = m2)
  expect_identical(s2$keep, "c6")
  expect_identical(s2$how, "QC pass and cell_meta qc_pass")
  # explicit cells: no QC, in the sample's order
  expect_identical(feature_cells(all, qct, cells = c("c5", "c4", "zz"), cell_meta = m2)$keep, c("c4", "c5"))
})

test_that("cell_features: a cell_meta qc_pass column does not override the `qc` argument", {
  # make_cell_meta() writes qc_pass with its own QC arguments; a stricter `qc` of the run still applies
  d = sim_dir("sim_norm_f")
  m = make_cell_meta(d)
  expect_true("qc_pass" %in% names(m))
  qc = list(min_frag = 20000)
  qct = do.call(cell_qc, c(list(read_sample_depth(d)), qc))
  want = intersect(qct[qc_pass == TRUE, cell], m[qc_pass %in% TRUE, cell])
  expect_lt(length(want), sum(m$qc_pass))                            # the stricter QC removes cells
  f = suppressWarnings(cell_features(d, cell_meta = m, sex = "female", centromeres = cs_centromeres(sim_norm_f$scn$bins),
                                     qc = qc, classes = c("arm", "chrom")))
  expect_setequal(f$cells$cell, want)
})

test_that("cell_features: the sample's sex from a cells table must be one sex", {
  s = sim_norm_f
  s$scn$cells = copy(s$scn$cells)[, sex := rep(c("F", "male"), length.out = .N)]
  d = file.path(tempdir(), "mixed_sex")
  unlink(d, recursive = TRUE)
  invisible(capture.output(suppressMessages(write_sim_dataset(s, d))))
  expect_error(cell_features(d, centromeres = cs_centromeres(s$scn$bins), classes = c("arm", "chrom")),
               "cell_features: the cells' `sex` column mixes both sexes")
  # spellings normalised, missing / empty values ignored
  expect_identical(cells_sex(c("F", NA, "female", " ")), "female")
  expect_identical(cells_sex(c(NA, "M")), "male")
  expect_identical(cells_sex(c(NA, NA)), NA_character_)
  expect_error(cells_sex(c("m", "Female")), "mixes both sexes")
  expect_error(cells_sex(c("m", "x")), "`sex` must be")
})

test_that("call_scn.R --cell-meta without qc_pass applies the cell QC (as the R API does)", {
  skip_if_not(testthat::is_checking() || isTRUE(as.logical(Sys.getenv("SCPHASEDCN_TEST_CLI", "FALSE"))))
  script = system.file("scripts", "call_scn.R", package = "scphasedcn")
  rs = file.path(R.home("bin"), "Rscript")
  d = sim_dir("sim_norm_f")
  all_cells = fread(file.path(d, "cells_qc.tsv"))$cell
  api = norm_result()$cells$cell                                     # QC-pass cells
  expect_gt(length(setdiff(all_cells, api)), 0)                     # some cells of the simulation fail the QC
  cm = file.path(tempdir(), "cli_meta_noqc.tsv")
  fwrite(data.table(cell = all_cells), cm, sep = "\t")
  cen = file.path(tempdir(), "cli_centromeres.tsv")
  fwrite(cs_centromeres(sim_norm_f$scn$bins), cen, sep = "\t")
  out = file.path(tempdir(), "callscn_cli_qc_out")
  unlink(out, recursive = TRUE)
  log = suppressWarnings(system2(rs, c(script, "--dir", d, "--out", out, "--sex", "female", "--centromeres", cen,
                                       "--cell-meta", cm, "--no-embedding", "--cores", "1"), stdout = TRUE, stderr = TRUE))
  expect_null(attr(log, "status"))
  expect_setequal(fread(file.path(out, "cells.tsv"))$cell, api)
})

test_that("params.yaml round-trips empty vectors and strings with backslashes / quotes", {
  p = caller_defaults()
  p$tabix = "C:\\tools\\tabix \"v2\""
  p$tumour = clone_defaults()
  p$tumour$allelic_exclude = character()
  p$tumour$checkpoint = "ck\\run\t1.rds"
  p$normal = normal_arg_defaults()
  tf = tempfile(fileext = ".yaml")
  writeLines(c("# scphasedcn caller parameters", yaml_lines(p)), tf)
  r = caller_params(tf)
  expect_identical(r$tumour$allelic_exclude, character())             # not the default chrX, chrY
  expect_identical(r$tabix, p$tabix)
  expect_identical(r$tumour$checkpoint, p$tumour$checkpoint)
  expect_equal(r, p)
  # a null keeps the default, an empty list is empty
  writeLines(c("tumour:", "  allelic_exclude: null", "  mix_w: []"), tf)
  r2 = caller_params(tf)
  expect_identical(r2$tumour$allelic_exclude, c("chrX", "chrY"))
  expect_identical(r2$tumour$mix_w, numeric())
})

test_that("headerless one-column barcode files keep their first barcode (make_cell_meta, cli_cells, call_scn)", {
  f = tempfile(); writeLines(c("AAAC-1", "AAAG", "AAAT"), f)
  fh = tempfile(); writeLines(c("cell", "AAAC", "AAAT"), fh)
  expect_identical(cli_cells(f, c("AAAC-1", "AAAG", "AAAT", "CCCC")), c("AAAC-1", "AAAG", "AAAT"))
  expect_identical(cli_cells(fh, c("AAAC", "AAAG", "AAAT")), c("AAAC", "AAAT"))
  m = make_cell_meta(data.table(cell = c("AAAC", "AAAG", "AAAT"), n_clean = 1e4), flat = f)
  expect_identical(m$flat, c(TRUE, TRUE, TRUE))
  expect_identical(make_cell_meta(data.table(cell = c("AAAC", "AAAG", "AAAT"), n_clean = 1e4), flat = fh)$flat, c(TRUE, FALSE, TRUE))
  expect_identical(as_cell_vector(f), c("AAAC-1", "AAAG", "AAAT"))
  expect_identical(multiplet_table(f)$cell, c("AAAC-1", "AAAG", "AAAT"))
})

test_that("barcodes of every input lose a trailing -1; matched counts in the notes; a warning when none match", {
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  ff = fake_features(120, rg, seed = 9)
  ids = unique(ff$cell)
  fe = list(long = ff, regions = rg, cells = unique(ff[, .(cell, stratum, n_frag)]))
  run = function(...) call_scn(tempdir(), "normal", features = fe, embedding = FALSE, verbose = FALSE, keep_scores = FALSE,
                               caller_args = list(n_draws = 500), ...)
  r = run(null_cells = paste0(ids[1:100], "-1"), multiplets = c("cell_0101-1", "nope"))
  expect_setequal(r$cells[null_cell == TRUE, cell], ids[1:100])
  expect_true(r$cells[cell == "cell_0101", multiplet])
  expect_true(any(grepl("null_cells: 100 given, 100 matched", r$notes)))
  expect_true(any(grepl("multiplets: 2 given, 1 matched", r$notes)))
  expect_warning(run(multiplets = "AAAC-1"), "none of the 1 multiplets barcodes")
  expect_identical(cli_cells({ f = tempfile(); writeLines(c("cell_0002-1", "cell_0003"), f); f }, ids[1:4]), ids[2:3])
})

test_that("write_scphasedcn clears the files of a previous result and nothing else", {
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  ff = fake_features(100, rg, seed = 10)
  fe = list(long = ff, regions = rg, cells = unique(ff[, .(cell, stratum, n_frag)]))
  res = call_scn(tempdir(), "normal", features = fe, embedding = FALSE, verbose = FALSE, caller_args = list(n_draws = 500))
  out = file.path(tempdir(), "write_clear")
  unlink(out, recursive = TRUE); dir.create(out)
  stale = c("clones.tsv", "assignments.tsv", "scores.tsv.gz", "reference_pool.txt")
  for (f in c(stale, "my_notes.txt")) writeLines("old", file.path(out, f))
  dir.create(file.path(out, "plots"))
  files = write_scphasedcn(res, out)
  expect_false(any(file.exists(file.path(out, stale))))
  expect_true(file.exists(file.path(out, "my_notes.txt")))                 # not ours: kept
  expect_true(dir.exists(file.path(out, "plots")))
  expect_setequal(list.files(out), c(basename(files), "my_notes.txt", "plots"))
  expect_true(all(basename(files) %in% scphasedcn_output_files()))
})

test_that("write_scphasedcn: rewriting a read-back result into its own directory keeps its scores and features", {
  # result.rds is stored without scores / features, so they cannot be written again from a read-back object
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  ff = fake_features(100, rg, seed = 10)
  fe = list(long = ff, regions = rg, cells = unique(ff[, .(cell, stratum, n_frag)]))
  res = call_scn(tempdir(), "normal", features = fe, embedding = FALSE, verbose = FALSE, caller_args = list(n_draws = 500))
  out = file.path(tempdir(), "write_roundtrip")
  unlink(out, recursive = TRUE)
  write_scphasedcn(res, out, write_scores = TRUE, write_features = TRUE)
  big = file.path(out, c("scores.tsv.gz", "features.tsv.gz"))
  expect_true(all(file.exists(big)))
  md5 = unname(tools::md5sum(big))
  back = read_scphasedcn(out)
  expect_null(back$scores)
  expect_null(back$features)
  write_scphasedcn(back, out, write_scores = TRUE, write_features = TRUE)
  expect_true(all(file.exists(big)))
  expect_identical(unname(tools::md5sum(big)), md5)
  # an object that carries them: written again when asked, cleared as stale when not
  write_scphasedcn(res, out, write_scores = TRUE)
  expect_true(file.exists(big[1]))
  expect_false(file.exists(big[2]))
  write_scphasedcn(res, out)
  expect_false(file.exists(big[1]))
})

test_that("precomputed features: unused feature-step inputs warn; tumour-mode null_cells = 'select' warns", {
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  ff = fake_features(100, rg, seed = 12)
  fe = list(long = ff, regions = rg, cells = unique(ff[, .(cell, stratum, n_frag)]))
  cm = data.table(cell = unique(ff$cell), cell_type = "a", flat = TRUE)
  expect_warning(call_scn(tempdir(), "normal", features = fe, sex = "female", events = c(x = "chr8"), cell_meta = cm,
                          centromeres = "cen.tsv", embedding = FALSE, verbose = FALSE, caller_args = list(n_draws = 500)),
                 "not used: sex; events; centromeres \\(\\.\\.\\.\\); cell_meta column\\(s\\) flat")
  expect_silent(call_scn(tempdir(), "normal", features = fe, cell_meta = cm[, .(cell, cell_type)], centromeres = NULL,
                         embedding = FALSE, verbose = FALSE, caller_args = list(n_draws = 500)))
  # tumour mode: "select" is a normal-mode selection (the warning comes before any work; the call then fails on the toy)
  expect_warning(try(call_scn(tempdir(), "tumour", features = fe, null_cells = "select", verbose = FALSE), silent = TRUE),
                 "tumour mode fits its input null")
})
