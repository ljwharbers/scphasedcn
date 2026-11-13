# Cell metadata for the caller

Builds the `cell_meta` table of
[`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
(the `--cell-meta` TSV of `call_scn.R`). Every column but `cell` is
optional for the caller:

- `qc_pass`: the cells to call.
  [`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
  (and `call_scn.R --cell-meta`) without `cells` calls the cells of the
  table that pass
  [`cell_qc()`](https://ljwharbers.github.io/scphasedcn/reference/cell_qc.md)
  (with the run's `qc` arguments) and, with this column, also have
  `qc_pass` TRUE.

- `cell_type`: the depth reference of normal-tissue mode is matched by
  cell type (types with fewer than 20 cells, and cells with `NA`, use
  the reference of all cells). Without cell types the reference is the
  all-cell one.

- `flat`: cells known to be copy-number neutral (e.g. normal cells from
  an orthogonal assay). Normal mode uses them as the null cells with
  `null_cells = "features"`; without `flat`, use `null_cells = "select"`
  (`--null-cells select`): the caller selects its own null cells. Tumour
  mode fits its input depth / allelic null on them (when the features
  carry no null fit); its reference pool is not taken from them but is
  the normal clone the caller finds and cleans itself.

- `group`: tumour mode only, pseudobulk groups of the phase-block
  orientation (`groups = "group"`, `--groups group`), e.g. spatial
  regions or RNA clusters; NA = all-cell pseudobulk.

## Usage

``` r
make_cell_meta(
  cells_qc = NULL,
  cell_types = NULL,
  flat = NULL,
  groups = NULL,
  qc_args = list()
)
```

## Arguments

- cells_qc:

  `cells_qc.tsv` written by the extractor (`scphasedcn run`; path), its
  table, or the extraction directory; NULL = the cells of the other
  inputs

- cell_types:

  NULL, a path or table with columns `cell` and `cell_type` (else the
  first two columns), or a character vector named by cell

- flat:

  NULL, cell barcodes of known copy-number-flat cells (vector, or a file
  with one cell per line, without a header unless the first line is
  `cell`, or with a `cell` column)

- groups:

  NULL, a path or table with `cell` and `group` (else the first two
  columns), or a vector named by cell

- qc_args:

  named list of
  [`cell_qc()`](https://ljwharbers.github.io/scphasedcn/reference/cell_qc.md)
  arguments for `qc_pass` (e.g. `list(min_frag = 10000)`)

## Value

data.table `cell`, `cell_type`, `flat`, `group`, plus `qc_pass`
([`cell_qc()`](https://ljwharbers.github.io/scphasedcn/reference/cell_qc.md)
on `cells_qc`) when `cells_qc` is given

## Details

Barcodes are matched as given (a trailing `-1` is stripped from every
input, as the extractor does).

## Examples

``` r
qc = data.table::data.table(cell = c("AAAC", "AAAG", "AAAT"), n_clean = c(40000, 25000, 900))
make_cell_meta(qc, cell_types = c(AAAC = "epithelial", AAAG = "immune"), flat = "AAAG")
#>      cell qc_pass  cell_type   flat  group
#>    <char>  <lgcl>     <char> <lgcl> <char>
#> 1:   AAAC    TRUE epithelial  FALSE   <NA>
#> 2:   AAAG    TRUE     immune   TRUE   <NA>
#> 3:   AAAT   FALSE       <NA>  FALSE   <NA>
```
