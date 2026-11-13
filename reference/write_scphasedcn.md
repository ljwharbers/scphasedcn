# Write a caller result to a directory

Flat tables: `cells.tsv` (cell, qc, stratum, fold, null_cell, n_calls,
n_events, clone / posterior / doublet_like / ref_pool when present),
`calls.tsv` (cell- and clone-level calls: called rows in normal mode,
all inherited / confirmed rows in tumour mode), `pooled_calls.tsv`
(called pool rows), `events.tsv`, `thresholds.tsv`, `null_fit.tsv`
(depth t rows + rho), `false_calls.tsv` + `false_events.tsv` (normal
mode), tumour mode: `clones.tsv`, `assignments.tsv`,
`clone_segments.tsv`, `clone_profiles.tsv.gz`, `reference_pool.txt`;
`regions.tsv`; `result.rds` (the object without scores and features),
`params.yaml`, `summary.tsv`
([`scphasedcn_summary()`](https://ljwharbers.github.io/scphasedcn/reference/scphasedcn_summary.md))
and `NOTES.md`. Optional: `scores.tsv.gz`, `features.tsv.gz`. The files
of a previous result in `out` (the names above) are removed first, so a
re-run never leaves stale tables behind (e.g. `clones.tsv` of a tumour
run, `scores.tsv.gz` of a run with `write_scores`); other files in `out`
are left alone. An object without scores (or features), e.g. one read
back with
[`read_scphasedcn()`](https://ljwharbers.github.io/scphasedcn/reference/read_scphasedcn.md),
leaves an existing `scores.tsv.gz` (`features.tsv.gz`) in place: it
cannot write them again.

## Usage

``` r
write_scphasedcn(obj, out, write_scores = FALSE, write_features = FALSE)
```

## Arguments

- obj:

  `scphasedcn` from
  [`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)

- out:

  output directory (created)

- write_scores:

  also write the per cell x region scores

- write_features:

  also write the features long table

## Value

named character vector of the written paths, invisibly
