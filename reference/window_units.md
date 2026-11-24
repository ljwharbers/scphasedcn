# Window units of the votes

Groups the atoms into subwindows (the `size` grid of each chromosome cut
at every region start / end, so that each region is a union of
subwindows) and sums the votes per cell x subwindow x phase set: one
unit is a stretch of at most `size` bp in one phase set, inside which
the phase is trusted.

## Usage

``` r
window_units(votes, atoms, regions, size = 1e+06)
```

## Arguments

- votes:

  votes table
  ([`read_votes()`](https://ljwharbers.github.io/scphasedcn/reference/read_votes.md))

- atoms:

  atoms table

- regions:

  regions table (the subwindows are cut at their boundaries; only their
  chromosomes are used)

- size:

  window size (bp)

## Value

list `units` (cell, sw, nA, n; units with n \>= 1) and `windows` (sw,
chrom, start, end)
