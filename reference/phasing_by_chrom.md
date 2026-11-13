# Phasing regime per chromosome

Per chromosome of the votes, the number of phase sets needed to hold
`frac` of the chromosome's votes (largest first) and the regime by the
rule of
[`detect_phasing()`](https://ljwharbers.github.io/scphasedcn/reference/detect_phasing.md):
`"chromosome"` when at most 2 phase sets hold `frac` of the votes, else
`"block"`.
[`detect_phasing()`](https://ljwharbers.github.io/scphasedcn/reference/detect_phasing.md)
is the same rule on the median of `k` over chromosomes; this table shows
chromosomes whose phasing differs from the global regime (e.g. a hybrid
phase whose chrX is still in read-backed blocks, so that plain sums over
its unoriented blocks attenuate or cancel the chrX allele signal).

## Usage

``` r
phasing_by_chrom(votes, atoms, frac = 0.9)
```

## Arguments

- votes:

  votes table
  ([`read_votes()`](https://ljwharbers.github.io/scphasedcn/reference/read_votes.md))

- atoms:

  atoms table

- frac:

  vote fraction the counted phase sets must hold

## Value

data.table `chrom, n_votes, n_ps` (distinct phase sets), `k` (phase sets
holding `frac` of the votes), `regime` ("chromosome" / "block"), in the
chromosome order of `atoms`
