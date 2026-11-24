# Haplotype fractions of the windowed allele term

The fractions at which
[`window_allelic_features()`](https://ljwharbers.github.io/scphasedcn/reference/window_allelic_features.md)
evaluates its windowed LLR: `1 - eps` (the loss / CN-LOH rows of the
normal state set, always included, so that the scorer finds their column
whatever `eps`) and `params$win_f` (default 0.98). Only near-complete
LOH is specific at window level: cells of normal tissue show a mild
chromosome-scale allelic skew (best window fraction 0.67-0.80 on many
chromosomes in the deepest cells of a normal thyroid sample), which a
window term at gain fractions would call; at 0.98 no cell of that sample
has a chromosome with W \> 10, while 99 % of the COLO829BL rows with W
\> 20 lie on ASCAT.sc-HP LOH segments.

## Usage

``` r
win_f_grid(params = caller_defaults())
```

## Arguments

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  (`win_f`, `eps`)

## Value

sorted numeric vector
