# Collapse fragments into molecules

Fragments with haplotype A or B, sorted by cell, chromosome and start; a
new molecule starts when `start > cummax(end) + gap` of the previous
fragments of the same cell and chromosome; one vote per molecule for the
majority haplotype (ties dropped); `ps` is the first fragment's phase
set. The running maximum is computed per group with an integer offset
(exact for coordinates below 2^53 / number of groups).

## Usage

``` r
collapse_molecules(fr, gap = 2000, split_ps = FALSE)
```

## Arguments

- fr:

  fragments with `cell, chrom, start, end, hap, ps`; `hap` "A"/"B"
  (other labels dropped) or integer codes 1 / 2 (= A / B)

- gap:

  maximal gap (bp) between fragments of one molecule

- split_ps:

  one vote per molecule x phase set instead of one per molecule (default
  FALSE; the caller parameter `mol_split_ps`)

## Value

data.table `cell, chrom, start, end, ps, n_frag, hap` (hap "A"/"B"), one
row per molecule (per molecule x phase set with `split_ps`)

## Details

Phase sets are not oriented relative to each other, so a molecule whose
fragments lie in different phase sets (a molecule across a block
junction) mixes unrelated haplotype labels. With `split_ps = TRUE` the
molecules are formed as above and then split by phase set: the fragments
of each phase set of a molecule give their own vote (majority haplotype
of those fragments, ties dropped) in that phase set, with `start` /
`end` / `n_frag` of those fragments. Molecules within one phase set are
the same under both settings.
