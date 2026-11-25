# Clean the normal reference pool with the allele channel

Reference pool = cells of the normal clone(s) (no non-neutral segment
\>= `min_event` bp) with posterior at least `min_post`, not doublet-like
with a tumour clone in the mixture (a mixture of two normal clones is
not a reason to exclude a cell), and without a per-cell allelic call on
the clone events: for every tumour clone k and each of its imbalanced
non-neutral segments of \>= `min_event` bp, the cell's votes are
oriented by k's pseudobulk (so the tumour's major haplotype is "A") and
scored with the one-sided
`ll_allelic(nA, n, f_grid, rho, two_sided = FALSE)` (best excess of the
tumour haplotype over the grid); the cell's event score is the maximum
over these clone x event regions, and a call is event score \> `thr`.
Without `thr` (`null = "max"`), it is the type-7 `q` quantile over the
candidate cells of the same maximum over the control regions (neutral in
every clone; orientation by the same clones), so the multiplicity of the
maximum is in the null as well; `null = "region"` is the quantile of the
single-region control scores (the first build's rule, which flags too
many cells when there are many events).

## Usage

``` r
clean_reference(
  feat,
  assign,
  segments,
  thr = NULL,
  par = NULL,
  members = NULL,
  min_post = 0.9,
  min_event = 1e+07,
  q = 0.999,
  control = NULL,
  f_grid = c(0.6, 2/3, 0.75, 0.85, 0.95),
  null = c("max", "region")
)
```

## Arguments

- feat:

  scn_features list (or `clone_input`)

- assign:

  assignments
  ([`assign_cells()`](https://ljwharbers.github.io/scphasedcn/reference/assign_cells.md))

- segments:

  segments table

- thr:

  allelic call threshold (NULL = from the control-region null)

- par:

  data.table clone, c, a (attenuation; default a = 1)

- members:

  data.table cell, clone used for orientation (default the assigned
  cells)

- min_post:

  minimal posterior

- min_event:

  minimal length of a non-neutral segment that makes a clone non-normal

- q:

  null quantile

- control:

  optional control regions (default
  [`control_regions()`](https://ljwharbers.github.io/scphasedcn/reference/control_regions.md)
  with min_len 5 Mb)

- f_grid:

  major-haplotype fractions of the one-sided allelic score

- null:

  "max" or "region" (see above)

## Value

list(pool = cell names, table = per-cell event / control max score and
flag, thr, normal = normal clones)
