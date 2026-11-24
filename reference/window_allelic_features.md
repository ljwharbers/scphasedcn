# Windowed, orientation-free allele features per cell and region

For each unit u of
[`window_units()`](https://ljwharbers.github.io/scphasedcn/reference/window_units.md)
(one cell x subwindow of at most `params$win_size` bp x phase set) and
each fraction f of
[`win_f_grid()`](https://ljwharbers.github.io/scphasedcn/reference/win_f_grid.md)
(default 0.98: near-complete LOH), the mirrored beta-binomial LLR
[`mirrored_llr()`](https://ljwharbers.github.io/scphasedcn/reference/mirrored_llr.md)
`L_u(f) = log(0.5 BB(a_u; n_u, f) + 0.5 BB(a_u; n_u, 1 - f)) - log BB(a_u; n_u, 1/2)`
(overdispersion `rho`, allele error floor `eps`), summed over the units
of each region: `W(f) = sum_u L_u(f)`. The phase is trusted only inside
a unit, so neither the orientation of the cell's haplotypes relative to
other cells nor phase switches between windows matter: an allelic
imbalance private to one cell is visible on chromosome-scale (hybrid)
phasing, where the oriented chromosome sums of
[`allelic_features()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_features.md)
cancel it at population-phase switches. A unit with one vote contributes
exactly 0, and under the null model `E[exp W(f)] = 1`, so
`P(W(f) > t) <= exp(-t)` (the null-cell guard of
[`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md);
the bound assumes independent votes inside a window, which holds well at
f near 1 but not at gain fractions, see
[`win_f_grid()`](https://ljwharbers.github.io/scphasedcn/reference/win_f_grid.md)).

## Usage

``` r
window_allelic_features(
  votes,
  atoms,
  regions,
  params = caller_defaults(),
  rho = NULL,
  exclude_atoms = NULL
)
```

## Arguments

- votes:

  votes table

- atoms:

  atoms table

- regions:

  regions table (`allelic_ok` regions are used)

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  (`eps`, `win_f`, `win_size`, `win_bias_max`, `win_bias_min_n`,
  `win_chrx`, `seed`)

- rho:

  NULL (fitted:
  [`fit_bb_rho()`](https://ljwharbers.github.io/scphasedcn/reference/fit_bb_rho.md)
  on the units with n \>= 2 of the null regions, at most 2e6 seeded
  units) or a value

- exclude_atoms:

  NULL or atom ids whose votes are dropped (allelic mask)

## Value

data.table `cell` (index), `region`, `n_win` (units with n \>= 2), one
`wal_f<fff>` column per fraction of
[`win_f_grid()`](https://ljwharbers.github.io/scphasedcn/reference/win_f_grid.md)
([`win_col()`](https://ljwharbers.github.io/scphasedcn/reference/win_col.md));
regions without a unit of the cell are absent. Attributes: `rho`,
`windows` (the subwindows with their pooled `nA`, `n`, `f` and
`masked`), `f_grid`, `regions` (the regions in scope)

## Details

Needs phase sets that are reliable over a window: on block-only phasing
a window holds many short unoriented blocks and the units are blocks of
one or two molecules (the per-block mirrored statistic, which has no
power). Subwindows whose votes over all cells have a shared allelic bias
(`|f - 0.5| > params$win_bias_max` on at least `params$win_bias_min_n`
votes; germline CNV, mapping bias, imprinting) are left out. The same
mask hides a clonal LOH that more than about `2 x win_bias_max` (20 %)
of the cells carry in one orientation of the population phase:
acceptable in normal mode, whose window term is for events of single
cells (shared events are the pools' job). Only autosomes are scored
unless `params$win_chrx` (female chrX carries X-inactivation driven
imbalance).
