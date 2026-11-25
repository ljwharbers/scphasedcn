# State set of normal mode

Resolves `params$states_normal`
([`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md))
to a
[`caller_states()`](https://ljwharbers.github.io/scphasedcn/reference/caller_states.md)
set: `"normal"`, `"normal_att"` (opt-in: loss and CN-LOH over
`f_att_normal`) or `"auto"`: `normal_att` when the features' phasing
regime (`phasing["detected"]`, the run-level rule of
[`detect_phasing()`](https://ljwharbers.github.io/scphasedcn/reference/detect_phasing.md)
that
[`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)
records in `features$phasing`) is `chromosome`, i.e. merged
chromosome-scale phase sets, whose switch errors attenuate the oriented
fraction of a real event under any orientation; `normal` on block
phasing, without a phasing record and without allele votes (the allelic
terms are 0, so both sets score alike). Mixed regimes (some chromosomes
in blocks,
[`phasing_by_chrom()`](https://ljwharbers.github.io/scphasedcn/reference/phasing_by_chrom.md))
follow the run-level regime: the state set is one per run, not per
chromosome.

## Usage

``` r
normal_state_set(params = caller_defaults(), phasing = NULL, votes = TRUE)
```

## Arguments

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  (`states_normal`, `f_att_normal`)

- phasing:

  the features' phasing record (`features$phasing`: named `detected`,
  `orientation`) or NULL

- votes:

  FALSE when the run has no allele votes

## Value

list(set = "normal" or "normal_att", requested, reason): `reason` is the
text of the run note
