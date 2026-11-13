# Headline numbers of a caller result

Headline numbers of a caller result

## Usage

``` r
scphasedcn_summary(obj)
```

## Arguments

- obj:

  `scphasedcn`

## Value

data.table section, item, value, n (denominator or count), note.
Sections: `cells` (cells, null cells, multiplets, cells with a
reportable call / event), `calls` (called cell / pool / clone rows per
level, class and state; `n` = reportable among them), `false_calls`
(null false calls per 1,000 per class over all strata; `n` = null rows),
`false_events` (false autosomal arm calls per null cell, any state and
per state), `events` (reportable events per cell and per null cell,
pooled calls; normal mode with
[`event_fdr()`](https://ljwharbers.github.io/scphasedcn/reference/event_fdr.md)
output: `cell_events_in_null_cells`, `expected_false_cell_events` (per
chromosome group when there is more than one) and
`expected_false_cell_event_fraction` (lower bounds, see
[`event_fdr()`](https://ljwharbers.github.io/scphasedcn/reference/event_fdr.md)),
`cell_events_q_below_0.05` (a descriptive count at a row-level q: not a
false discovery rate and not a filter, see
[`event_fdr()`](https://ljwharbers.github.io/scphasedcn/reference/event_fdr.md)),
`empirical_p_floor`; `calls` then also has
`expected_false_reportable_cell_calls`), `null_check` (selection-free
exceedances per 1,000 per class, an upper bound of the null rate: real
events count too,
[`null_selection_check()`](https://ljwharbers.github.io/scphasedcn/reference/null_selection_check.md)),
`clones` (clones, assigned cells, reference pool).
