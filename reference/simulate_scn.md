# Simulate a SPLONGGET-like single-nucleus DNA data set

Simulate a SPLONGGET-like single-nucleus DNA data set

## Usage

``` r
simulate_scn(
  n_cells = 200,
  genome = sim_genome(),
  sex = c("male", "female"),
  cell_types = c(normal = 0.6, tumour = 0.4),
  depth = c(meanlog = log(3000), sdlog = 0.35),
  events = sim_events_tumour(),
  private_frac = 0.05,
  multiplet_rate = 0.15,
  multiplet_sizes = c(`2` = 0.8, `3` = 0.2),
  peak_frac = 0.5,
  peak_sd = 200,
  access_sd = 0.5,
  wave_sd = 0.15,
  noise_phi = 0.25,
  frag_len = c(meanlog = log(700), sdlog = 0.6),
  snp_density = 1/2000,
  phase_block_mean = 1500000,
  switch_rate = 1e-04,
  block_flip_prob = 0,
  allele_error = 0.04,
  seed = 1
)
```

## Arguments

- n_cells:

  number of barcodes (cells) in the output

- genome:

  bins table from
  [`sim_genome()`](https://ljwharbers.github.io/scphasedcn/reference/sim_genome.md)

- sex:

  "male" (chrX one copy, hap A) or "female"

- cell_types:

  named numeric vector of type proportions

- depth:

  lognormal parameters (meanlog, sdlog) of fragments per nucleus

- events:

  data.table (chrom, start, end, cnA, cnB, type); `type` may be a cell
  type name, "private" (then `private_frac` of the nuclei carry it) or
  "nucleus:*,,..." (exactly the given nuclei carry it, e.g. several
  events in one nucleus: a highly aneuploid cell)*

- private_frac:

  fraction of nuclei carrying a "private" event

- multiplet_rate:

  fraction of barcodes that are multi-nucleus pools

- multiplet_sizes:

  sizes (number of nuclei) and their probabilities

- peak_frac:

  fraction of fragments placed at accessibility peaks (0 = uniform
  placement within a bin)

- peak_sd:

  spread (bp) of fragment starts around a peak

- access_sd:

  log-sd of the shared per-bin accessibility

- wave_sd:

  log-sd of the cell-type chromatin wave (2-5 Mb scale)

- noise_phi:

  cell-private multiplicative noise (log-sd), autocorrelated over ~5
  bins

- frag_len:

  lognormal parameters of fragment length (bp)

- snp_density:

  het SNPs per bp (uniform positions, one SNP per position: a draw that
  falls on an occupied position is moved to the next free one)

- phase_block_mean:

  mean phase-block length (bp)

- switch_rate:

  per-SNP switch-error probability inside a block (1e-4 keeps a 750-SNP
  block mostly in one parity, attenuation ~0.93 by
  [`switch_attenuation()`](https://ljwharbers.github.io/scphasedcn/reference/switch_attenuation.md);
  0.02 scrambles it)

- block_flip_prob:

  probability that a block's orientation is flipped relative to the true
  haplotypes (0 = anchored blocks)

- allele_error:

  per-observation allele error rate (ONT)

- seed:

  RNG seed of the simulation (cells, nuclei, copy numbers, fragments,
  SNPs, alleles). It does not change the genome: `genome` is evaluated
  before seeding, and
  [`sim_genome()`](https://ljwharbers.github.io/scphasedcn/reference/sim_genome.md)
  uses its own `seed`. After `set.seed(seed)` the simulator first makes
  (and discards) the random draws of
  [`sim_genome()`](https://ljwharbers.github.io/scphasedcn/reference/sim_genome.md)
  for `nrow(genome)` bins, then draws the data; the result depends only
  on the genome table, `seed` and the other parameters, not on how
  `genome` was given. The caller's random stream is saved and restored.

## Value

list with `scn` (an `scn` object, see
[`new_scn()`](https://ljwharbers.github.io/scphasedcn/reference/new_scn.md)),
`truth` (cells, nuclei, cn/cnA per nucleus, events), `fragments`
(data.table in the Python `fragments.tsv.gz` layout, `ps` = phase set of
the first covered SNP or -1, plus nucleus, copy, hap_true, hap, n_snps),
`snps` (with hap_of_alt_true), `params`
