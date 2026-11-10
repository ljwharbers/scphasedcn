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
heterozygous SNVs of a phased VCF (Clair3 + longphase) for `--snps`. `scphasedcn index-check` checks the BAM index
and can write a fresh one to a separate path (`--reindex`, then `--index` in `extract` / `run`).
See `scphasedcn <command> --help` and the scphasedcn documentation for the input requirements and the phasing recipe.

Options worth knowing:

- `--barcode-map`: a TSV with columns `atac_barcode` and `rna_barcode`. A whitelist in RNA-barcode space (10x
  Multiome GEX barcodes) is translated to the ATAC barcodes of the BAM's `CB` tag.
- `--peaks`: an optional peak BED file (a file, not a directory) for FRiP per cell and `peak_bp` per bin (QC only).
- `--layers`: a comma list of `depth`, `alleles`, `fragments` (default: all three, which is what the R package
  reads; the per-cell QC tables are written with `depth` and `alleles`, and the old name `qc` is accepted and
  ignored). `--layers alleles,fragments` recomputes the allele layer alone, e.g. with a new VCF.
- `--snps`: a bgzipped, tabix-indexed VCF whose contig names match the BAM (`chr2`, not `2`); both are checked before
  the BAM pass. Phased genotypes without a `PS` field (e.g. population phasing) form one phase set per contig.
  A VCF with several samples needs `--vcf-sample` (yaml: `vcf_sample`).
- Cell barcodes: a trailing `-1` is optional on both sides (whitelist and BAM `CB`). A chromosome where no record
  is in a whitelisted cell while many are in other barcodes stops the extraction (wrong whitelist).
- `merge` fails when a chromosome has no finished partial (`--allow-missing` merges without it) and when the partials
  were extracted from different inputs (cells, bin size, BAM, SNP file, masks: recorded in each partial's `DONE`).
  `--region` writes `partials/region_<chrom>_<start>-<end>`, which `merge` does not merge.

Genome: GRCh38 / hg38 with `chr` names (bundled chromosome sizes, ENCODE blacklist v2 and centromeres; see
`src/scphasedcn/resources/README.md`).
