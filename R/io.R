# io.R
# What: the `scn` container (bins x cells count matrix + layers + obs/var +
#       obsm + alleles) and readers / writers for the extraction-layer data
#       contract: bins_100kb.h5ad (anndata HDF5, gzip only), alleles_dna.h5ad (writer only: simulated data sets),
#       fragments.tsv.gz (+ tabix), cells_qc.tsv, and the sample YAML.
# Why: one object feeds every module; readers tolerate missing layers so the
#      R side runs on partial extractions and on simulated data written by
#      write_bins_h5ad() (round-trip tested).
# Inputs: file paths. Outputs: scn objects, data.tables, files on disk.
# Layout assumptions (anndata >= 0.8): X and layers are CSR groups
# (data/indices/indptr, attr shape), obs/var are dataframe groups (attr
# `_index`), obsm are 2-D arrays (n_obs x k, C order), uns scalars/arrays.
# rhdf5 reverses array dimensions relative to h5py, so 2-D obsm arrays are
# transposed on read and write.

#' Construct an scn object
#'
#' @param counts bins x cells matrix (coerced to dgCMatrix)
#' @param bins data.table of bins (chrom, start, end, ...)
#' @param cells data.table of cells (column `cell` first)
#' @param layers named list of bins x cells matrices (same dims as counts)
#' @param obsm named list of cells x k matrices
#' @param uns list of unstructured metadata
#' @param alleles NULL or list(AD, DP (snps x cells dgCMatrix), snps data.table)
#' @param fragments NULL, a path to fragments.tsv.gz, or a data.table
#' @return object of class `scn`
#' @export
new_scn = function(counts, bins, cells, layers = list(), obsm = list(), uns = list(), alleles = NULL, fragments = NULL) {
  counts = as_dgc(counts)
  bins = as.data.table(bins); cells = as.data.table(cells)
  stopifnot(nrow(bins) == nrow(counts), nrow(cells) == ncol(counts))
  if (!"cell" %in% names(cells)) stop("cells needs a `cell` column")
  if (!"bin_id" %in% names(bins)) bins[, bin_id := seq_len(.N)]
  layers = lapply(layers, function(m) {
    m = as_dgc(m)
    stopifnot(all(dim(m) == dim(counts)))
    m
  })
  obj = list(counts = counts, bins = bins, cells = cells, layers = layers, obsm = obsm, uns = uns,
             alleles = alleles, fragments = fragments)
  class(obj) = "scn"
  obj
}

#' @export
print.scn = function(x, ...) {
  cat(sprintf("scn object: %d bins x %d cells\n", nrow(x$counts), ncol(x$counts)))
  cat(sprintf("  chromosomes: %s\n", paste(unique(x$bins$chrom), collapse = " ")))
  cat(sprintf("  layers: %s\n", if (length(x$layers)) paste(names(x$layers), collapse = " ") else "(none)"))
  cat(sprintf("  obsm: %s\n", if (length(x$obsm)) paste(names(x$obsm), collapse = " ") else "(none)"))
  cat(sprintf("  cells columns: %s\n", paste(names(x$cells), collapse = " ")))
  if (!is.null(x$alleles)) cat(sprintf("  alleles: %d SNPs\n", nrow(x$alleles$snps)))
  if (!is.null(x$fragments)) cat(sprintf("  fragments: %s\n", if (is.character(x$fragments)) x$fragments else sprintf("%d rows", nrow(x$fragments))))
  invisible(x)
}

#' Subset an scn object by cells and/or bins
#' @param scn scn object
#' @param cells logical/integer/character selector of cells (NULL = all)
#' @param bins logical/integer selector of bins (NULL = all)
#' @return scn
#' @export
subset_scn = function(scn, cells = NULL, bins = NULL) {
  ci = if (is.null(cells)) seq_len(ncol(scn$counts)) else if (is.character(cells)) match(cells, scn$cells$cell) else which(seq_len(ncol(scn$counts)) %in% seq_len(ncol(scn$counts))[cells])
  bi = if (is.null(bins)) seq_len(nrow(scn$counts)) else seq_len(nrow(scn$counts))[bins]
  out = scn
  out$counts = scn$counts[bi, ci, drop = FALSE]
  out$layers = lapply(scn$layers, function(m) m[bi, ci, drop = FALSE])
  out$bins = scn$bins[bi]
  out$cells = scn$cells[ci]
  out$obsm = lapply(scn$obsm, function(m) m[ci, , drop = FALSE])
  if (!is.null(scn$alleles)) {
    out$alleles$AD = scn$alleles$AD[, ci, drop = FALSE]
    out$alleles$DP = scn$alleles$DP[, ci, drop = FALSE]
  }
  out
}

# ---- HDF5 helpers ---------------------------------------------------------------

h5_attr_scalar = function(value, file, name, attr) {
  fid = rhdf5::H5Fopen(file)
  on.exit(rhdf5::H5Fclose(fid), add = TRUE)
  root = name == "/" || name == ""
  if (!root && !rhdf5::H5Lexists(fid, name)) stop("no such object: ", name)
  obj = if (root) fid else rhdf5::H5Oopen(fid, name)
  if (!root) on.exit(rhdf5::H5Oclose(obj), add = TRUE)
  if (is.character(value) && length(value) == 1) {
    rhdf5::h5writeAttribute(value, obj, attr, asScalar = TRUE, variableLengthString = TRUE, encoding = "UTF-8")
  } else if (is.character(value)) {
    rhdf5::h5writeAttribute(value, obj, attr, variableLengthString = TRUE, encoding = "UTF-8")
  } else rhdf5::h5writeAttribute(value, obj, attr)
  invisible(TRUE)
}

h5_write_encoding = function(file, name, type, version = "0.1.0") {
  h5_attr_scalar(type, file, name, "encoding-type")
  h5_attr_scalar(version, file, name, "encoding-version")
}

# write a bins x cells matrix as an anndata CSR group (cells x bins): the CSC
# representation of bins x cells is exactly the CSR of its transpose
h5_write_csr = function(m, file, name) {
  m = as_dgc(m)
  rhdf5::h5createGroup(file, name)
  rhdf5::h5write(as.numeric(m@x), file, paste0(name, "/data"))
  rhdf5::h5write(as.integer(m@i), file, paste0(name, "/indices"))
  rhdf5::h5write(as.integer(m@p), file, paste0(name, "/indptr"))
  h5_write_encoding(file, name, "csr_matrix")
  h5_attr_scalar(as.integer(c(ncol(m), nrow(m))), file, name, "shape")
  invisible(TRUE)
}

h5_read_sparse = function(file, name, n_obs = NULL, n_var = NULL) {
  at = tryCatch(rhdf5::h5readAttributes(file, name), error = function(e) list())
  if (!is.null(at$shape)) {
    shape = as.integer(at$shape)
  } else shape = c(n_obs, n_var)
  enc = as.character(at[["encoding-type"]] %||% "csr_matrix")
  data = as.numeric(rhdf5::h5read(file, paste0(name, "/data"), bit64conversion = "double"))
  indices = as.integer(rhdf5::h5read(file, paste0(name, "/indices"), bit64conversion = "double"))
  indptr = as.integer(rhdf5::h5read(file, paste0(name, "/indptr"), bit64conversion = "double"))
  if (enc == "csr_matrix") {
    # rows = obs (cells), columns = var (bins): CSC of bins x cells
    m = new("dgCMatrix", i = indices, p = indptr, x = data, Dim = c(shape[2], shape[1]))
  } else {
    m = new("dgCMatrix", i = indices, p = indptr, x = data, Dim = c(shape[1], shape[2]))
    m = t(m)   # csc of obs x var -> transpose to var x obs
  }
  m
}

# full paths and object types of every entry in an HDF5 file
h5_entries = function(file) {
  ls = rhdf5::h5ls(file)
  data.table(path = sub("^//", "/", paste0(ls$group, "/", ls$name)), otype = ls$otype, group = ls$group, name = ls$name)
}

h5_read_matrix_any = function(file, name, n_obs, n_var) {
  en = h5_entries(file)
  entry = en[path == paste0("/", sub("^/", "", name))]
  if (nrow(entry) && entry$otype[1] == "H5I_GROUP") return(h5_read_sparse(file, name, n_obs, n_var))
  d = rhdf5::h5read(file, name, bit64conversion = "double")
  as_dgc(matrix(as.numeric(d), nrow = n_var, ncol = n_obs))  # dense obs x var in C order -> var x obs in R
}

h5_write_df = function(df, file, name, index_col) {
  rhdf5::h5createGroup(file, name)
  cols = setdiff(names(df), index_col)
  rhdf5::h5write(as.character(df[[index_col]]), file, paste0(name, "/", index_col))
  for (cn in cols) {
    v = df[[cn]]
    path = paste0(name, "/", cn)
    if (is.factor(v)) v = as.character(v)
    if (is.logical(v)) { rhdf5::h5write(as.integer(v), file, path); h5_write_encoding(file, path, "array", "0.2.0"); next }
    if (is.character(v)) {
      v[is.na(v)] = "NA"
      rhdf5::h5write(v, file, path)
      h5_write_encoding(file, path, "string-array", "0.2.0")
      next
    }
    rhdf5::h5write(as.numeric(v), file, path)
    h5_write_encoding(file, path, "array", "0.2.0")
  }
  h5_write_encoding(file, name, "dataframe", "0.2.0")
  h5_attr_scalar(index_col, file, name, "_index")
  h5_attr_scalar(if (length(cols)) cols else numeric(0), file, name, "column-order")
  invisible(TRUE)
}

# an HDF5 vector as rhdf5 reads it, made a plain R vector: no dim, raw -> integer, and an enum {FALSE, TRUE} (how h5py /
# anndata write a boolean array; read as a factor) -> logical. Other enums keep their labels (character).
h5_plain = function(v) {
  if (is.factor(v) && length(levels(v)) && all(levels(v) %in% c("FALSE", "TRUE"))) return(as.logical(as.character(v)))
  if (is.array(v)) v = as.vector(v)
  if (is.raw(v)) v = as.integer(v)
  v
}

h5_read_df = function(file, name) {
  at = tryCatch(rhdf5::h5readAttributes(file, name), error = function(e) list())
  index_col = as.character(at[["_index"]] %||% "_index")
  ls = rhdf5::h5ls(file)
  here = ls[ls$group == paste0("/", name), ]
  out = list()
  for (k in seq_len(nrow(here))) {
    nm = here$name[k]
    path = paste0(name, "/", nm)
    if (here$otype[k] == "H5I_GROUP") {
      sub = ls[ls$group == paste0("/", path), "name"]
      if (all(c("codes", "categories") %in% sub)) {
        codes = as.integer(rhdf5::h5read(file, paste0(path, "/codes"), bit64conversion = "double"))
        cats = as.character(rhdf5::h5read(file, paste0(path, "/categories")))
        v = cats[codes + 1L]
        v[codes < 0] = NA
      } else if (all(c("values", "mask") %in% sub)) {
        # nullable integer / boolean / string array: values with a mask of the missing entries
        v = h5_plain(rhdf5::h5read(file, paste0(path, "/values"), bit64conversion = "double"))
        mask = as.logical(h5_plain(rhdf5::h5read(file, paste0(path, "/mask"))))
        v[mask %in% TRUE] = NA
      } else next
    } else {
      v = h5_plain(rhdf5::h5read(file, path, bit64conversion = "double"))
    }
    out[[nm]] = v
  }
  dt = as.data.table(out)
  order_attr = as.character(at[["column-order"]] %||% character())
  ord = c(index_col, intersect(order_attr, names(dt)))
  setcolorder(dt, c(intersect(ord, names(dt)), setdiff(names(dt), ord)))
  dt
}

h5_read_obsm = function(file, name, n_obs) {
  ls = rhdf5::h5ls(file)
  here = ls[ls$group == paste0("/", name), ]
  out = list()
  for (k in seq_len(nrow(here))) {
    if (here$otype[k] != "H5I_DATASET") next
    m = rhdf5::h5read(file, paste0(name, "/", here$name[k]), bit64conversion = "double")
    m = as.matrix(m)
    if (nrow(m) != n_obs && ncol(m) == n_obs) m = t(m)
    out[[here$name[k]]] = m
  }
  out
}

# ---- bins h5ad ---------------------------------------------------------------------

#' Write an scn object as `bins_100kb.h5ad` (anndata layout, gzip only)
#' @param scn scn object
#' @param file output path
#' @param layers which layers to write (default all)
#' @return file path (invisibly)
#' @export
write_bins_h5ad = function(scn, file, layers = names(scn$layers)) {
  if (file.exists(file)) file.remove(file)
  rhdf5::h5createFile(file)
  h5_write_encoding(file, "/", "anndata", "0.1.0")
  h5_write_csr(scn$counts, file, "X")
  h5_write_df(scn$cells, file, "obs", "cell")
  bins = copy(scn$bins)
  bins[, bin_id := NULL]
  bins[, bin := sprintf("%s:%.0f-%.0f", chrom, start, end)]
  setcolorder(bins, "bin")
  h5_write_df(bins, file, "var", "bin")
  rhdf5::h5createGroup(file, "layers")
  for (nm in layers) h5_write_csr(scn$layers[[nm]], file, paste0("layers/", nm))
  rhdf5::h5createGroup(file, "obsm")
  for (nm in names(scn$obsm)) {
    rhdf5::h5write(t(as.matrix(scn$obsm[[nm]])), file, paste0("obsm/", nm))   # R (k x n_obs) -> h5py (n_obs x k)
    h5_write_encoding(file, paste0("obsm/", nm), "array", "0.2.0")
  }
  rhdf5::h5createGroup(file, "uns")
  if (length(scn$uns$params)) {
    rhdf5::h5createGroup(file, "uns/params")
    for (nm in names(scn$uns$params)) {
      v = scn$uns$params[[nm]]
      if (is.logical(v)) v = as.integer(v)
      if (is.list(v)) next
      rhdf5::h5write(v, file, paste0("uns/params/", nm))
    }
  }
  rhdf5::h5closeAll()
  invisible(file)
}

#' Read `bins_100kb.h5ad` into an scn object (missing layers tolerated)
#' @param file path
#' @param layers layers to read (default: all present)
#' @param cells optional character vector of cells to keep
#' @return scn object
#' @export
read_bins_h5ad = function(file, layers = NULL, cells = NULL) {
  obs = h5_read_df(file, "obs")
  var = h5_read_df(file, "var")
  setnames(obs, 1, "cell")
  n_obs = nrow(obs); n_var = nrow(var)
  X = h5_read_matrix_any(file, "X", n_obs, n_var)
  ls = rhdf5::h5ls(file)
  present = ls[ls$group == "/layers", "name"]
  want = if (is.null(layers)) present else intersect(layers, present)
  missing = setdiff(layers %||% character(), present)
  if (length(missing)) message("layers not in file (skipped): ", paste(missing, collapse = ", "))
  lay = list()
  for (nm in want) lay[[nm]] = h5_read_matrix_any(file, paste0("layers/", nm), n_obs, n_var)
  obsm = if ("obsm" %in% ls[ls$group == "/", "name"]) h5_read_obsm(file, "obsm", n_obs) else list()
  uns = list()
  if (any(ls$group == "/uns/params")) {
    for (nm in ls[ls$group == "/uns/params", "name"]) uns$params[[nm]] = tryCatch(as.vector(rhdf5::h5read(file, paste0("uns/params/", nm))), error = function(e) NULL)
  }
  rhdf5::h5closeAll()
  bins = copy(var)
  setnames(bins, 1, "bin")
  if ("width" %in% names(bins)) bins[, width := NULL]
  if (!all(c("chrom", "start", "end") %in% names(bins))) {
    p = tstrsplit(bins$bin, "[:-]")
    bins[, chrom := p[[1]]][, start := as.numeric(p[[2]])][, end := as.numeric(p[[3]])]
  }
  bins[, start := as.numeric(start)][, end := as.numeric(end)]
  bins[, bin_id := seq_len(.N)]
  if ("is_centromere" %in% names(bins)) bins[, is_centromere := as.logical(is_centromere)]
  colnames(X) = obs$cell
  for (nm in names(lay)) colnames(lay[[nm]]) = obs$cell
  scn = new_scn(counts = X, bins = bins, cells = obs, layers = lay, obsm = obsm, uns = uns)
  if (!is.null(cells)) scn = subset_scn(scn, cells = cells)
  scn
}

# ---- alleles h5ad -----------------------------------------------------------------

#' Write the allele counts as `alleles_dna.h5ad` (layers AD, DP; var = SNPs)
#' @param alleles list(AD, DP, snps) (snps x cells matrices)
#' @param cells data.table with a `cell` column (same order as the matrix columns)
#' @param file output path
#' @export
write_alleles_h5ad = function(alleles, cells, file) {
  if (file.exists(file)) file.remove(file)
  rhdf5::h5createFile(file)
  h5_write_encoding(file, "/", "anndata", "0.1.0")
  h5_write_csr(alleles$DP, file, "X")
  h5_write_df(cells[, .(cell)], file, "obs", "cell")
  snps = copy(alleles$snps)
  if (!"snp" %in% names(snps)) snps[, snp := paste0(chrom, ":", pos)]
  setcolorder(snps, "snp")
  h5_write_df(snps, file, "var", "snp")
  rhdf5::h5createGroup(file, "layers")
  h5_write_csr(alleles$AD, file, "layers/AD")
  h5_write_csr(alleles$DP, file, "layers/DP")
  rhdf5::h5closeAll()
  invisible(file)
}

# ---- fragments ----------------------------------------------------------------------

#' Write a fragments table as bgzipped TSV with a tabix index
#'
#' Columns written: chrom start end cell strand mapq flags (+ any extra columns
#' present, e.g. hap, n_snps). The header line starts with `#` so tabix skips it.
#' @param fr data.table
#' @param file output path ending in `.tsv.gz`
#' @param extra extra columns to keep (default: any beyond the contract)
#' @return path
#' @export
write_fragments = function(fr, file, extra = NULL) {
  core = c("chrom", "start", "end", "cell", "strand", "mapq", "flags")
  extra = extra %||% setdiff(intersect(names(fr), c("len", "piece", "npieces", "hap", "nsnp", "ps", "read_id")), core)
  d = copy(fr)[, c(core, extra), with = FALSE]
  d[, start := as.integer(start)][, end := as.integer(end)]
  setorderv(d, c("chrom", "start"))
  tmp = sub("\\.gz$", "", file)
  hdr = paste0("#", paste(names(d), collapse = "\t"))
  writeLines(hdr, tmp)
  fwrite(d, tmp, sep = "\t", col.names = FALSE, append = TRUE)
  if (requireNamespace("Rsamtools", quietly = TRUE)) {
    if (file.exists(file)) file.remove(file)
    Rsamtools::bgzip(tmp, dest = file, overwrite = TRUE)
    file.remove(tmp)
    Rsamtools::indexTabix(file, seq = 1, start = 2, end = 3, skip = 0, comment = "#", zeroBased = TRUE)
  } else {
    system2("gzip", c("-f", tmp))
    warning("Rsamtools not available: wrote plain gzip without a tabix index")
  }
  invisible(file)
}

#' Read fragments (whole file or tabix regions)
#' @param file fragments.tsv.gz
#' @param regions optional data.table (chrom, start, end) from [parse_region()]; uses the tabix index
#' @param cells optional character vector of cells to keep
#' @return data.table chrom, start, end, cell, strand, mapq, flags (+ extras), hap decoded from flags when absent
#' @export
read_fragments = function(file, regions = NULL, cells = NULL) {
  con = gzfile(file)
  hdr = readLines(con, n = 1)
  close(con)
  cols = if (startsWith(hdr, "#")) strsplit(sub("^#", "", hdr), "\t")[[1]] else c("chrom", "start", "end", "cell", "strand", "mapq", "flags")
  if (!is.null(regions) && file.exists(paste0(file, ".tbi")) && requireNamespace("Rsamtools", quietly = TRUE)) {
    tb = Rsamtools::TabixFile(file)
    rr = copy(as.data.table(regions))
    cs = tryCatch(Rsamtools::seqnamesTabix(tb), error = function(e) NULL)
    if (!is.null(cs)) rr = rr[chrom %in% cs]     # a contig absent from the index has no fragments (as without tabix)
    rr[is.na(end), end := 536870911]            # tabix limit (2^29 - 1), beyond any chromosome
    lines = character()
    if (nrow(rr)) {
      # overlapping / adjacent query ranges merged; a fragment spanning two disjoint ranges is returned by both
      # queries and kept from the first only (start before the end of the previous range on the same contig)
      gr = GenomicRanges::reduce(GenomicRanges::GRanges(rr$chrom, IRanges::IRanges(rr$start + 1, rr$end)))
      res = Rsamtools::scanTabix(tb, param = gr)
      gch = as.character(GenomicRanges::seqnames(gr)); gend = GenomicRanges::end(gr)
      js = match("start", cols); if (is.na(js)) js = 2L
      for (i in seq_along(res)) {
        li = res[[i]]
        if (i > 1 && length(li) && gch[i] == gch[i - 1]) {
          st = as.numeric(vapply(strsplit(li, "\t", fixed = TRUE), `[`, "", js))
          li = li[!(st < gend[i - 1])]
        }
        lines = c(lines, li)
      }
    }
    d = if (length(lines)) fread(text = lines, header = FALSE, col.names = cols) else as.data.table(setNames(replicate(length(cols), character()), cols))
  } else {
    d = fread(cmd = paste("gzip -dc", shQuote(file)), header = FALSE, skip = if (startsWith(hdr, "#")) 1 else 0, col.names = cols)
    if (!is.null(regions)) {
      keep = rep(FALSE, nrow(d))
      for (r in seq_len(nrow(regions))) {
        e = regions$end[r]; if (is.na(e)) e = Inf
        keep = keep | (d$chrom == regions$chrom[r] & d$end > regions$start[r] & d$start < e)
      }
      d = d[keep]
    }
  }
  if (!is.null(cells)) d = d[cell %in% cells]
  if (nrow(d)) {
    d[, start := as.numeric(start)][, end := as.numeric(end)]
    if ("hap" %in% names(d)) d[, hap := hap_label(hap)] else d[, hap := "U"]
    if (!"len" %in% names(d)) d[, len := end - start]
  }
  d[]
}

#' Read `cells_qc.tsv`
#' @param file path
#' @return data.table with column `cell` first
#' @export
read_cells_qc = function(file) {
  d = fread(file)
  if (!"cell" %in% names(d)) setnames(d, 1, "cell")
  setcolorder(d, "cell")
  d[]
}

#' Write a simulated data set in the extraction-layer layout
#' @param sim result of [simulate_scn()]
#' @param dir output directory
#' @return named list of paths
#' @export
write_sim_dataset = function(sim, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  paths = list(bins = file.path(dir, "bins_100kb.h5ad"), alleles = file.path(dir, "alleles_dna.h5ad"),
               fragments = file.path(dir, "fragments.tsv.gz"), cells_qc = file.path(dir, "cells_qc.tsv"),
               truth = file.path(dir, "truth_cells.tsv"))
  write_bins_h5ad(sim$scn, paths$bins)
  write_alleles_h5ad(sim$scn$alleles, sim$scn$cells, paths$alleles)
  fr = copy(sim$fragments)
  fr[, cell := sim$scn$cells$cell[cell]]
  write_fragments(fr, paths$fragments)
  fwrite(sim$scn$cells, paths$cells_qc, sep = "\t")
  fwrite(sim$truth$cells, paths$truth, sep = "\t")
  paths
}

# ---- config -----------------------------------------------------------------------------

#' Read a YAML config (sample yaml or caller yaml)
#'
#' Uses the yaml package when installed, otherwise a minimal parser for the
#' subset used in the sample and caller yaml files (nested maps, scalars,
#' inline `[a, b]` lists and `{k: v}` maps, `#` comments). The sample yaml layout is
#' `system.file("extdata", "sample_template.yaml", package = "scphasedcn")`.
#' @param file path
#' @return nested list
#' @export
#' @examples
#' str(read_config(system.file("extdata", "sample_template.yaml", package = "scphasedcn"))$samples)
read_config = function(file) {
  if (length(file) != 1 || !nzchar(file) || !file.exists(file)) stop("read_config: no such file: '", file, "'")
  if (requireNamespace("yaml", quietly = TRUE)) return(yaml::read_yaml(file))
  parse_simple_yaml(readLines(file))
}

parse_scalar = function(s) {
  s = trimws(s)
  s = sub("\\s+#.*$", "", s)
  if (grepl("^\\[.*\\]$", s)) {
    inner = trimws(substr(s, 2, nchar(s) - 1))
    if (!nzchar(inner)) return(list())
    return(lapply(trimws(strsplit(inner, ",")[[1]]), parse_scalar))
  }
  if (grepl("^\\{.*\\}$", s)) {
    inner = trimws(substr(s, 2, nchar(s) - 1))
    items = trimws(strsplit(inner, ",")[[1]])
    out = list()
    for (it in items) {
      kv = strsplit(it, ":", fixed = TRUE)[[1]]
      out[[trimws(kv[1])]] = parse_scalar(paste(kv[-1], collapse = ":"))
    }
    return(out)
  }
  s = gsub('^"|"$', "", s)
  s = gsub("^'|'$", "", s)
  if (s %in% c("true", "True", "TRUE")) return(TRUE)
  if (s %in% c("false", "False", "FALSE")) return(FALSE)
  if (grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", s)) return(as.numeric(s))
  s
}

parse_simple_yaml = function(lines) {
  lines = lines[!grepl("^\\s*(#|$)", lines)]
  indent = nchar(sub("^(\\s*).*$", "\\1", lines))
  parse_block = function(from, to, level) {
    out = list()
    k = from
    while (k <= to) {
      if (indent[k] != level) { k = k + 1; next }
      line = trimws(lines[k])
      key = sub(":.*$", "", line)
      val = sub("^[^:]*:\\s*", "", line)
      # children: following lines with deeper indent
      nxt = k + 1
      while (nxt <= to && indent[nxt] > level) nxt = nxt + 1
      if (nxt > k + 1) {
        out[[key]] = parse_block(k + 1, nxt - 1, indent[k + 1])
      } else out[[key]] = parse_scalar(val)
      k = nxt
    }
    out
  }
  parse_block(1, length(lines), 0)
}
