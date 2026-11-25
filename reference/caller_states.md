# Copy-number state sets of the scorer

One row per (state, allelic configuration). A state may have several
rows (e.g. gain = 2:1 at f 2/3 and 3:1 at f 3/4); its allelic LLR is the
best row (two-sided), its depth LLR the row's log2 shift.

- `normal`: neutral 1:1; loss 1:0 at f = 1 - eps (shift loss); gain 2:1
  (2/3) and 3:1 (3/4), both at the gain shift; cnloh 2:0 at 1 - eps
  (shift 0). The normal-tissue state set.

- `normal_att`: as `normal` with the loss and cnloh rows over
  `f_att_normal` (the attenuation-aware opt-in of normal mode,
  `params$states_normal`; `f_att` when the parameters have no
  `f_att_normal`), always together with `1 - eps`, the fraction of the
  `normal` rows (so the set contains the `normal` loss / cnloh rows, and
  the windowed allele term of
  [`score_cells()`](https://ljwharbers.github.io/scphasedcn/reference/score_cells.md)
  always finds its row). Loss and cnloh keep one shared row set, so
  their allelic class stays `loss_or_cnloh`.

- `tumour` (per-cell detection and assignment): loss 1:0 at f in c(0.80,
  f_att); cnloh 2:0 over f_att; gain 2:1 (2/3); gain_as 3:0 over f_att
  (gain shift); amp 3:1 (0.75), 2:2 (0.5), 4:1 (0.80) with shift
  `k_amp * log2(cn / 2)`.

## Usage

``` r
caller_states(
  set = c("normal", "normal_att", "tumour"),
  params = caller_defaults()
)
```

## Arguments

- set:

  "normal", "normal_att" or "tumour"

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

## Value

data.table state, row_id, nA, nB, cn, shift (log2), f (major-haplotype
fraction); the first row is `neutral`
