# scphasedcn (Python extractor)

Extraction layer of the [scphasedcn](https://github.com/ljwharbers/scphasedcn) single-cell copy-number caller (the R package of the same name). One pass
per chromosome over a SPLONGGET DNA BAM (10x Multiome barcodes in the `CB` tag, Oxford Nanopore reads) writes the
directory the scphasedcn R package reads:

- `bins_100kb.h5ad`: fragment counts per cell x 100 kb bin, with GC, blacklist / tandem-repeat fractions and
  centromere flags per bin;
- `fragments.tsv.gz` + `.tbi`: one row per fragment with cell, flags, haplotype (from phased heterozygous SNPs) and
  phase set;
- `cells_qc.tsv`: per-cell QC; `alleles_dna.h5ad`: allele counts per cell x SNP; `summary.json`.

```bash
pip install ./python            # or: the scphasedcn container / conda environment
scphasedcn run --bam sample.bam --barcodes barcodes.tsv.gz --fasta GRCh38.fa \
    --snps het_phased.vcf.gz --phased --out results/sample/extract --jobs 8
```

`scphasedcn run` extracts every chromosome in parallel and merges; `scphasedcn extract --chrom <c>` and
`scphasedcn merge` are the same steps for one cluster job per chromosome. `scphasedcn select-hets` selects the
heterozygous SNVs of a phased VCF (Clair3 + longphase) for `--snps`. `scphasedcn hybrid-phase` orients those
read-backed phase blocks with a population-phased VCF of the same sample and writes one phase set per chromosome
(recommended for the normal-tissue mode of the caller; see below). `scphasedcn index-check` checks the BAM index
and can write a fresh one to a separate path (`--reindex`, then `--index` in `extract` / `run`).
See `scphasedcn <command> --help` and the scphasedcn documentation for the input requirements and the phasing recipe.

Options worth knowing:

- `--barcode-map`: a TSV with columns `atac_barcode` and `rna_barcode`. A whitelist in RNA-barcode space (10x
  Multiome GEX barcodes) is translated to the ATAC barcodes of the BAM's `CB` tag.
- `--peaks`: an optional peak BED file (a file, not a directory) for FRiP per cell and `peak_bp` per bin (QC only).
- `--layers`: a comma list of `depth`, `alleles`, `fragments`, `qc` (default: all four, which is what the R package
  reads). `--layers alleles,fragments` recomputes the allele layer alone, e.g. with a new VCF.

`scphasedcn hybrid-phase`:

```bash
scphasedcn hybrid-phase --blocks hets.vcf.gz --population popphased_{chrom}.vcf.gz --out hybrid.vcf.gz \
    --block-tsv hybrid_blocks.tsv --chrom-tsv hybrid_chroms.tsv
```

- `--blocks`: the read-backed phased VCF (GT `|` + PS, e.g. longphase or the `select-hets` output).
- `--population`: the same sample's genotypes phased with a reference panel (e.g. Beagle 5 or SHAPEIT), one VCF or a
  pattern with `{chrom}` for one file per chromosome. Contig names may differ by the `chr` prefix; plain, gzip or bgzip,
  with or without an index.
- A site votes when it is a heterozygous biallelic SNV with the same REF / ALT in both files; population markers with
  the INFO flag `--imputed-flag` (default `IMP`) never vote. Each block takes the majority orientation and joins the
  chromosome's phase set; blocks with a tie or fewer than `--min-votes` votes keep their own PS. `--fill-unphased`
  gives hets outside the blocks the population phase. INFO `PHASE_SRC` / `POP_VOTE` record the source and the vote of
  each site; the TSV reports give per block and per chromosome the votes, blocks oriented / flipped, the
  switch-disagreement rate and the minority fraction.

Genome: GRCh38 / hg38 with `chr` names (bundled chromosome sizes, ENCODE blacklist v2 and centromeres; see
`src/scphasedcn/resources/README.md`).
