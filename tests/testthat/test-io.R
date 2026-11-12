# test-io.R: h5ad / fragments / config round trips and tolerance to missing layers.
test_that("bins h5ad round trip is exact and tolerates missing layers", {
  dir = file.path(tempdir(), "sim_io")
  paths = write_sim_dataset(sim_small, dir)
  scn2 = read_bins_h5ad(paths$bins)
  expect_equal(as.matrix(scn2$counts), as.matrix(sim_small$scn$counts), ignore_attr = TRUE)
  expect_equal(as.matrix(scn2$layers$hapA), as.matrix(sim_small$scn$layers$hapA), ignore_attr = TRUE)
  expect_equal(as.matrix(scn2$layers$span), as.matrix(sim_small$scn$layers$span), ignore_attr = TRUE)
  expect_equal(scn2$obsm$fraglen_hist, sim_small$scn$obsm$fraglen_hist, ignore_attr = TRUE)
  expect_identical(scn2$cells$cell, sim_small$scn$cells$cell)
  expect_identical(scn2$bins$chrom, sim_small$scn$bins$chrom)
  expect_equal(scn2$bins$start, sim_small$scn$bins$start)
  expect_true(is.logical(scn2$bins$is_centromere))
  # the caller's sample reader: same counts and cells, the fragments path and the cells_qc.tsv columns
  scn3 = scphasedcn:::read_sample_depth(dir, layers = "span")
  expect_equal(as.matrix(scn3$counts), as.matrix(sim_small$scn$counts), ignore_attr = TRUE)
  expect_identical(scn3$cells$cell, sim_small$scn$cells$cell)
  expect_setequal(names(scn3$layers), "span")
  expect_identical(scn3$fragments, paths$fragments)
  expect_true(all(setdiff(names(sim_small$scn$cells), "cell") %in% names(scn3$cells)))
  # subset of layers written; asking for an absent layer is tolerated
  f2 = file.path(dir, "partial.h5ad")
  write_bins_h5ad(sim_small$scn, f2, layers = c("span", "frag_long"))
  expect_message(s3 <- read_bins_h5ad(f2, layers = c("span", "hapA")), "not in file")
  expect_setequal(names(s3$layers), "span")
  expect_null(s3$layers$hapA)
})

test_that("fragments are written with tabix and read by region", {
  skip_if_not_installed("Rsamtools")
  dir = file.path(tempdir(), "sim_io")
  if (!file.exists(file.path(dir, "fragments.tsv.gz"))) write_sim_dataset(sim_small, dir)
  fr = read_fragments(file.path(dir, "fragments.tsv.gz"), regions = parse_region("chr21:0-2000000"))
  expect_true(all(fr$chrom == "chr21"))
  expect_true(all(fr$start < 2e6))
  expect_true(all(fr$hap %in% c("U", "A", "B", "N")))
  all_fr = read_fragments(file.path(dir, "fragments.tsv.gz"))
  expect_equal(nrow(all_fr), nrow(sim_small$fragments))
  expect_true(all(c("len", "nsnp", "ps") %in% names(all_fr)))
})

test_that("read_fragments with tabix: each fragment once, as without the index; absent contigs give no rows", {
  skip_if_not_installed("Rsamtools")
  dir = file.path(tempdir(), "sim_io")
  if (!file.exists(file.path(dir, "fragments.tsv.gz"))) write_sim_dataset(sim_small, dir)
  f = file.path(dir, "fragments.tsv.gz")
  f2 = file.path(tempdir(), "frag_noindex.tsv.gz")             # the same file without its index: the fallback reader
  file.copy(f, f2, overwrite = TRUE)
  key = function(x) x[order(chrom, start, end, cell)]
  # adjacent and disjoint ranges: fragments spanning a boundary once
  rg = parse_region(c("chr21:0-1000000", "chr21:1000000-2000000", "chr21:2000500-3000000", "chr21:5000000-6000000"))
  a = read_fragments(f, regions = rg)
  b = read_fragments(f2, regions = rg)
  expect_gt(nrow(a), 0)
  expect_equal(key(a), key(b), ignore_attr = TRUE)
  # a contig absent from the index: no rows (also next to a present one)
  expect_equal(nrow(read_fragments(f, regions = parse_region("chrNOPE"))), 0)
  one = read_fragments(f, regions = parse_region("chr21:0-1000000"))
  expect_equal(read_fragments(f, regions = rbind(parse_region("chrNOPE"), parse_region("chr21:0-1000000"))), one)
  expect_equal(key(one), key(read_fragments(f2, regions = parse_region("chr21:0-1000000"))), ignore_attr = TRUE)
})

test_that("the fallback yaml parser handles the config subset", {
  txt = c("# comment", "reference:", "  fasta: /a/b.fa", "outputs:", "  root: /x", "samples:", "  S1:",
          "    sex: male", "    chroms: [chr1, chrX]", "    events: {g20: \"chr20:28000000-\", c5: chr5}", "    n: 3")
  f = tempfile(fileext = ".yaml"); writeLines(txt, f)
  cfg = scphasedcn:::parse_simple_yaml(readLines(f))
  expect_equal(cfg$reference$fasta, "/a/b.fa")
  expect_equal(cfg$samples$S1$sex, "male")
  expect_equal(unlist(cfg$samples$S1$chroms), c("chr1", "chrX"))
  expect_equal(cfg$samples$S1$events$g20, "chr20:28000000-")
  expect_equal(cfg$samples$S1$n, 3)
  if (requireNamespace("yaml", quietly = TRUE)) {
    ref = yaml::read_yaml(f)
    expect_equal(ref$samples$S1$events$c21, cfg$samples$S1$events$c21)
  }
})
