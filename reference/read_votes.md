# Molecule votes per cell, atom and phase set

Per chromosome: `tabix <file> <chrom> | cut -f2,3,4,7,11,13` read with
`fread(cmd = )` (an awk filter keeps hap 1 / 2 with `ps > 0` before R
sees the rows), or
[`read_fragments()`](https://ljwharbers.github.io/scphasedcn/reference/read_fragments.md)
when the tabix binary or index is missing; fragments of `cells` with
haplotype A / B and a phase set are collapsed into molecules
([`collapse_molecules()`](https://ljwharbers.github.io/scphasedcn/reference/collapse_molecules.md)),
each molecule is assigned to the atom containing its start and the votes
are aggregated per cell x atom x phase set. Memory stays per chromosome
(integer cell indices, [`gc()`](https://rdrr.io/r/base/gc.html) between
chromosomes). A failing tabix command stops with its error; a chromosome
with atoms but no haplotype-assigned fragment of `cells` gives a
warning.

## Usage

``` r
read_votes(
  frag_file,
  cells,
  atoms,
  chroms = NULL,
  gap = 2000,
  tabix = caller_defaults()$tabix
)
```

## Arguments

- frag_file:

  fragments.tsv.gz (+ .tbi)

- cells:

  cell barcodes; the votes' `cell` is the index into this vector

- atoms:

  [`atoms_from_regions()`](https://ljwharbers.github.io/scphasedcn/reference/atoms_from_regions.md)
  output

- chroms:

  chromosomes to read (default: the atoms' chromosomes)

- gap:

  molecule gap (bp)

- tabix:

  path of the tabix binary ("" = the `tabix` on the PATH; NULL / not
  found = the R fallback through Rsamtools)

## Value

votes data.table `cell, atom, ps, nA, nB` (passes
[`check_votes()`](https://ljwharbers.github.io/scphasedcn/reference/check_tables.md))
