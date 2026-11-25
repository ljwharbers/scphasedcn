# Mini genome for simulations

Mini genome for simulations

## Usage

``` r
sim_genome(
  chroms = c(chr1 = 6e+07, chr2 = 4e+07, chr8 = 3e+07, chr21 = 2e+07, chr22 = 2e+07, chrX
    = 4e+07),
  bin_size = 1e+05,
  centromere_frac = 0.4,
  seed = 1
)
```

## Arguments

- chroms:

  named numeric vector of chromosome lengths (bp)

- bin_size:

  bin width in bp

- centromere_frac:

  position of the centromere as a fraction of length

- seed:

  RNG seed for the bin covariates: the genome depends on `seed` only
  (not on the caller's random stream), and the caller's random stream is
  left as it was (the RNG state is saved and restored), also when
  `sim_genome()` is evaluated as the `genome` argument of
  [`simulate_scn()`](https://ljwharbers.github.io/scphasedcn/reference/simulate_scn.md)

## Value

data.table of bins: bin_id, chrom, start, end, gc, frac_blacklist,
frac_trf, is_centromere, peak_bp, arm, n_peaks
