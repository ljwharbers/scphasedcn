# Per-cell per-region features of a sample

Reads an extraction directory, applies the cell QC
([`cell_qc()`](https://ljwharbers.github.io/scphasedcn/reference/cell_qc.md);
with a `qc_pass` column of `cell_meta` a cell must also have `qc_pass`
TRUE; `cells` restricts the run to a set of cells instead, without QC),
annotates arms, filters bins (default:
[`bin_filter()`](https://ljwharbers.github.io/scphasedcn/reference/bin_filter.md)
on the kept cells, chr1-22 and chrX), normalises depth
([`normalise_depth()`](https://ljwharbers.github.io/scphasedcn/reference/normalise_depth.md):
`celltype` for normal tissue, `model` for tumours), builds the regions
([`region_set()`](https://ljwharbers.github.io/scphasedcn/reference/region_set.md)),
the depth features
([`depth_features()`](https://ljwharbers.github.io/scphasedcn/reference/depth_features.md)),
the molecule votes
([`read_votes()`](https://ljwharbers.github.io/scphasedcn/reference/read_votes.md))
and the oriented allele features
([`allelic_features()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_features.md)),
and joins them into the long features table
([`vignette("scphasedcn-outputs", package = "scphasedcn")`](https://ljwharbers.github.io/scphasedcn/articles/scphasedcn-outputs.md)).

## Usage

``` r
cell_features(
  dir,
  cells = NULL,
  cell_meta = NULL,
  method = c("celltype", "model"),
  ref_cells = NULL,
  min_cells_type = 20,
  use = NULL,
  sex = NULL,
  centromeres = NULL,
  classes = c("arm", "chrom", "5Mb", "10Mb", "20Mb", "50Mb"),
  events = NULL,
  null_exclude = character(),
  phasing = "auto",
  groups = NULL,
  mol_gap = params$mol_gap,
  keep_votes = TRUE,
  layers = character(),
  qc = list(),
  params = caller_defaults(),
  verbose = FALSE,
  exclude_atoms = NULL
)
```

## Arguments

- dir:

  extraction directory (bins_100kb.h5ad, fragments.tsv.gz + .tbi,
  cells_qc.tsv)

- cells:

  NULL (QC-pass cells) or the cells to use (no QC)

- cell_meta:

  NULL or a data.table with `cell` and optional
  `qc_pass, cell_type, flat, fold, group, stratum, n_frag` (used as-is;
  NA cell_type = no type-matched reference; NA or missing `stratum` /
  `n_frag` keep the QC values). With `cells` NULL only cells of
  `cell_meta` are used, and with a `qc_pass` column (e.g.
  [`make_cell_meta()`](https://ljwharbers.github.io/scphasedcn/reference/make_cell_meta.md))
  only those with `qc_pass` TRUE that also pass
  [`cell_qc()`](https://ljwharbers.github.io/scphasedcn/reference/cell_qc.md)
  with the `qc` arguments

- method:

  depth normalisation, "celltype" or "model" (any
  [`normalise_depth()`](https://ljwharbers.github.io/scphasedcn/reference/normalise_depth.md)
  method is accepted)

- ref_cells:

  reference cells of the normalisation (selector; NULL = all)

- min_cells_type:

  celltype: minimal reference cells of a type

- use:

  NULL or a logical bin filter

- sex:

  "female" / "male" (case-insensitive, or "F" / "M"; any other value
  stops; NULL: the cells' `sex` column, whose non-missing values must
  all give the same sex)

- centromeres:

  centromere table or path (Battenberg layout, see
  [`arm_table()`](https://ljwharbers.github.io/scphasedcn/reference/arm_table.md);
  NULL: the hg38 table shipped with the package,
  [`hg38_centromeres()`](https://ljwharbers.github.io/scphasedcn/reference/hg38_centromeres.md))

- classes:

  region classes

- events:

  named events (see
  [`region_set()`](https://ljwharbers.github.io/scphasedcn/reference/region_set.md))

- null_exclude:

  chromosomes / regions without null rows (known events)

- phasing:

  "auto", "chromosome", "block" (= "loo"), "loo" or "mirrored"

- groups:

  NULL, the name of a `cell_meta` column, a vector named by cell, or one
  label per kept cell: the pseudobulk groups of the block orientation
  (NA = all-cell pseudobulk)

- mol_gap:

  molecule gap (bp; default `params$mol_gap`, 2000)

- keep_votes:

  keep the votes table in the result

- layers:

  bin layers to read (none are needed by the depth normalisation)

- qc:

  named list of
  [`cell_qc()`](https://ljwharbers.github.io/scphasedcn/reference/cell_qc.md)
  arguments

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- verbose:

  log each step (the run notes) as it happens

- exclude_atoms:

  NULL or atom ids whose votes are left out of the allele features (see
  [`allelic_features()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_features.md))

## Value

`scn_features` list: `long` (features, keyed by cell and region),
`regions`, `cells` (row order = the votes' cell index), `votes`,
`atoms`, `win`
([`window_matrices()`](https://ljwharbers.github.io/scphasedcn/reference/window_matrices.md)),
`params`, `notes`, `use`, `phasing`
