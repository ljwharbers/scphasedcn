# Random pools of cells within groups

Shuffles the cells of each group (cell type) and splits them into
`round(n / n_pool)` pools of (nearly) equal size; groups with fewer than
`ceiling(n_pool / 2)` cells give no pool. Groups are drawn (and listed)
in C-locale order of their names, so a seed gives the same pools under
any collation.

## Usage

``` r
make_pools(cells, group = NULL, n_pool = 20, seed = 1)
```

## Arguments

- cells:

  cell ids

- group:

  group per cell (NULL = one group)

- n_pool:

  target pool size

- seed:

  RNG seed

## Value

data.table pool, cell, group
