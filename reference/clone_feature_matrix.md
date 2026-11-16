# Cells x features matrix for clone discovery

One depth column (`d|<region>`: the window z, capped at +/- `cap`, NA
-\> 0) and one signed allelic column (`a|<region>`) per eligible region
of the given classes. Allelic z =
`(nA - n / 2) / sqrt(n / 4 * (1 + (n - 1) * rho))` on the oriented
votes, capped at +/- `cap`, 0 below `min_n` molecules; the sign keeps
which haplotype is in excess, so subclones that are imbalanced on
opposite haplotypes separate (e.g. two subclones with a loss of
different chr21 homologues). Input whose `orient` is not in
`keep_orient` (unoriented block-phased `none`, or `clone` = oriented by
known groups such as spatial regions, which would leak the grouping into
discovery) is first re-oriented by the all-cell leave-one-out pseudobulk
(`orient_blocks`, NULL groups) when `votes`, `atoms`, `regions` and
`cells` are given.

## Usage

``` r
clone_feature_matrix(
  feat_long,
  classes = c("arm", "20Mb"),
  rho,
  cap = 6,
  min_n = 5,
  votes = NULL,
  atoms = NULL,
  regions = NULL,
  cells = NULL,
  keep_orient = c("anchored", "pseudobulk"),
  allele = c("z", "fraction")
)
```

## Arguments

- feat_long:

  features long table
  ([`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)
  `long`)

- classes:

  region classes used

- rho:

  beta-binomial overdispersion of the molecule votes

- cap:

  cap of \|z\|

- min_n:

  minimal molecules for an allelic value

- votes, atoms, regions, cells:

  optional: orient unoriented input (cells = names in vote index order)

- keep_orient:

  `orient` labels used as given (chromosome-wide anchored phase,
  all-cell pseudobulk)

- allele:

  scale of the allelic columns: `"z"` (default) = the standardised count
  above, whose expectation in an imbalanced region grows like
  `(f - 0.5) sqrt(n)`, so deep and shallow cells of one clone separate;
  `"fraction"` = coverage-free: the cell's oriented fraction `nA / n`
  minus 0.5, in units of the null SD at the region's median molecule
  count `n0` (`sqrt(0.25 (1 + (n0 - 1) rho) / n0)`; same scale as `"z"`
  for a typical cell, expectation independent of the cell's depth), and
  cells below `min_n` molecules get the region's mean (no information)
  instead of 0 (balanced)

## Value

numeric matrix cells x features (rownames = cells)
