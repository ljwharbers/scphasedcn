#!/usr/bin/env Rscript
# check_extraction_io.R
# What: reads a real extraction directory written by the Python extractor with the R package's readers and checks shapes,
#       columns and a few known values.
# Why: the two halves of scphasedcn meet at the extraction directory (bins h5ad, fragment table + tabix index, cells QC);
#      the R tests use simulated data written by R, so only this check reads what the Python side writes.
# Usage: Rscript tools/check_extraction_io.R <extraction dir>   (or the directory in SCPHASEDCN_EXTRACT_DIR)
#        The directory is the Python test suite's synthetic sample (python/tests/conftest.py, fixture `extracted`),
#        written by tools/tiny_extraction.py with the same sample at 50 kb bins in <extraction dir>_50kb; the expected
#        values below are documented there.
# Inputs: <dir>/bins_100kb.h5ad (found with bins_file()), fragments.tsv.gz (+ .tbi), cells_qc.tsv;
#         <dir>_50kb/bins_50kb.h5ad
# Outputs: none (stops at the first failed check; prints "extraction io ok" at the end)

suppressPackageStartupMessages({
  if (!"scphasedcn" %in% loadedNamespaces()) library(scphasedcn)
  library(data.table)
})
args = commandArgs(trailingOnly = TRUE)
dir = if (length(args)) args[1] else Sys.getenv("SCPHASEDCN_EXTRACT_DIR")
if (!nzchar(dir) || !dir.exists(dir)) stop("usage: Rscript check_extraction_io.R <extraction dir>")
ns = asNamespace("scphasedcn")
check = function(ok, what) {
  if (!isTRUE(ok)) stop("check failed: ", what, call. = FALSE)
  cat("ok  ", what, "\n")
}
CELL_A = "AAAACCCCGGGGTTTT"
CELL_B = "CCCCGGGGTTTTAAAA"

# ---- bins h5ad --------------------------------------------------------------------------------------------------
check(identical(bins_file(dir), file.path(dir, "bins_100kb.h5ad")), "bins_file: bins_100kb.h5ad")
scn = read_bins_h5ad(bins_file(dir))
check(inherits(scn, "scn"), "read_bins_h5ad returns an scn object")
check(identical(colnames(scn$counts), sort(c(CELL_A, CELL_B))), "two cells, in the file's order")
check(identical(scn$cells$cell, colnames(scn$counts)), "cells table matches the matrix columns")
check(nrow(scn$bins) == 19 && nrow(scn$counts) == 19, "19 bins (chr1 10, chr2 5, chrM 1, chrX 3)")
check(all(c("bin", "chrom", "start", "end", "bin_id", "gc", "frac_blacklist", "is_centromere") %in% names(scn$bins)),
      "bins columns")
check(is.numeric(scn$bins$start) && is.numeric(scn$bins$end) && is.logical(scn$bins$is_centromere), "bins column types")
check(identical(as.integer(table(scn$bins$chrom)[c("chr1", "chr2", "chrM", "chrX")]), c(10L, 5L, 1L, 3L)),
      "bins per chromosome")
b0 = which(scn$bins$chrom == "chr1" & scn$bins$start == 0)
check(scn$counts[b0, CELL_A] == 12, "depth: 12 fragments of cell A in chr1:0-100000")
check(all(c("span", "frag_long", "frag_short", "hapA", "hapB") %in% names(scn$layers)), "depth layers present")
check(all(vapply(scn$layers, function(m) identical(dim(m), dim(scn$counts)), logical(1))), "layer shapes")
check(sum(scn$layers$hapA) == 1 && sum(scn$layers$hapB) == 1, "one haplotype-A and one haplotype-B fragment")
check(all(c("n_pass", "n_records") %in% names(scn$cells)), "per-cell QC in obs")
check(scn$cells[cell == CELL_A, n_pass] == 20 && scn$cells[cell == CELL_B, n_pass] == 7, "n_pass per cell")

# ---- cells_qc.tsv and the caller's sample reader ----------------------------------------------------------------
qc = read_cells_qc(file.path(dir, "cells_qc.tsv"))
check(names(qc)[1] == "cell" && setequal(qc$cell, c(CELL_A, CELL_B)), "read_cells_qc: cell column first, two cells")
sd = ns$read_sample_depth(dir, layers = "span")
check(identical(sd$fragments, file.path(dir, "fragments.tsv.gz")), "read_sample_depth: fragment table found")
check(identical(sd$cells$cell, colnames(sd$counts)) && all(names(qc) %in% names(sd$cells)),
      "read_sample_depth: QC columns merged in matrix order")
check(identical(names(sd$layers), "span"), "read_sample_depth: requested layer only")
meta = make_cell_meta(cells_qc = dir)
check(setequal(meta$cell, c(CELL_A, CELL_B)) && "qc_pass" %in% names(meta), "make_cell_meta on the directory")

# ---- fragment table ---------------------------------------------------------------------------------------------
ff = file.path(dir, "fragments.tsv.gz")
fr = read_fragments(ff)
cols = c("chrom", "start", "end", "cell", "strand", "mapq", "flags", "len", "piece", "npieces", "hap", "nsnp", "ps", "read_id")
check(identical(names(fr), cols), "read_fragments: the extractor's columns, in order")
check(setequal(fr$cell, c(CELL_A, CELL_B)) && sum(fr$cell == CELL_A) == 20, "read_fragments: 20 fragments of cell A")
check(is.character(fr$hap) && all(fr$hap %in% c("U", "A", "B", "N")), "read_fragments: haplotype labels")
r1 = fr[chrom == "chr2" & start == 900]
check(nrow(r1) == 1 && r1$cell == CELL_A && r1$hap == "A" && r1$nsnp == 3 && r1$ps == 1000, "read_fragments: read R1")
reg = data.table(chrom = "chr2", start = 0, end = NA_real_)
fr2 = read_fragments(ff, regions = reg)
check(nrow(fr2) == 5 && all(fr2$chrom == "chr2"), "read_fragments with a region (tabix index): 5 chr2 fragments")
check(nrow(read_fragments(ff, regions = data.table(chrom = "chr1", start = 10000, end = 10001))) == 2,
      "read_fragments: 2 fragments overlap chr1:10000-10001")

# ---- allele votes: tabix command line and the R fallback agree --------------------------------------------------
cells = colnames(scn$counts)
atoms = atoms_from_regions(scn$bins, data.table(chrom = "chr2", start = 0, end = NA_real_))
check_atoms(atoms)
want = data.table(cell = match(c(CELL_A, CELL_B), cells), nA = c(1L, 0L), nB = c(0L, 1L))
v_r = read_votes(ff, cells, atoms, chroms = "chr2", tabix = NULL)
check_votes(v_r)
check(identical(v_r[order(cell), .(cell, nA, nB)], want[order(cell)]) && all(v_r$ps == 1000),
      "read_votes (R fallback): one A vote (cell A) and one B vote (cell B) in phase set 1000")
tbx = unname(Sys.which("tabix"))
if (nzchar(tbx)) {
  v_t = read_votes(ff, cells, atoms, chroms = "chr2", tabix = tbx)
  check(identical(v_t, v_r), "read_votes (tabix command line) equals the R fallback")
  hf = ns$read_hap_fragments(ff, "chr2", tabix = tbx)
  check(identical(names(hf), c("start", "end", "cell", "hap", "ps")) && nrow(hf) == 2, "read_hap_fragments (tabix)")
} else {
  cat("skip tabix not on the PATH: the tabix path of read_votes was not checked\n")
}
# ---- the same sample at 50 kb bins: bins_50kb.h5ad, found by the sample reader -------------------------------------
dir50 = paste0(sub("/+$", "", dir), "_50kb")
if (!dir.exists(dir50)) stop("no ", dir50, " (tools/tiny_extraction.py writes it next to the extraction directory)")
check(identical(bins_file(dir50), file.path(dir50, "bins_50kb.h5ad")), "bins_file: the one bins_50kb.h5ad")
sd50 = ns$read_sample_depth(dir50)
check(identical(as.integer(table(sd50$bins$chrom)[c("chr1", "chr2", "chrM", "chrX")]), c(20L, 10L, 1L, 6L)),
      "50 kb: 37 bins (chr1 20, chr2 10, chrM 1, chrX 6)")
check(all(sd50$bins[chrom != "chrM", end - start] == 50000), "50 kb: bin width")
check(sum(sd50$counts) == sum(scn$counts) && identical(sd50$fragments, file.path(dir50, "fragments.tsv.gz")),
      "50 kb: the same fragments as at 100 kb")
b50 = which(sd50$bins$chrom == "chr1" & sd50$bins$start %in% c(0, 50000))
check(identical(as.numeric(sd50$counts[b50, CELL_A]), c(8, 4)),
      "50 kb: cell A's 12 fragments of chr1:0-100000 split 8 + 4 at 50 kb")
cat("extraction io ok\n")
