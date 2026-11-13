# scphasedcn

Single-cell copy-number calling for **SPLONGGET**: single-cell long-read
DNA sequencing in which every read carries a 10x Multiome cell barcode
and is sequenced on Oxford Nanopore. scphasedcn combines read depth with
read-phased alleles (each molecule covering a phased heterozygous SNP
votes for a haplotype), per cell and region:

- **normal tissue**: rare, cell-private gains, losses and copy-neutral
  LOH, with thresholds calibrated on null cells (1 false call per 1,000
  null cell-regions), per cell and in cell pools;
- **tumours**: clone discovery, allele-specific clone copy-number
  profiles, likelihood assignment of every cell, per-cell inherited /
  confirmed calls, and input for MEDICC2 trees.

The repository holds two packages, both named scphasedcn, that are used
together:

| package | language | does |
|----|----|----|
| `scphasedcn` (Python, `python/`) | Python (pysam) | BAM -\> extraction directory: per-cell 100 kb bin counts, a fragment table with haplotypes, cell QC |
| `scphasedcn` (R, root) | R | extraction directory -\> calls, clones, tables |

## Install

Container (Docker / Apptainer / Singularity; R, Python, both packages
and htslib):

``` bash
apptainer pull docker://ghcr.io/ljwharbers/scphasedcn:0.1.1
apptainer exec scphasedcn_0.1.1.sif scphasedcn --help
```

Conda (from a checkout):

``` bash
git clone https://github.com/ljwharbers/scphasedcn && cd scphasedcn
conda env create -f environment.yml && conda activate scphasedcn      # or micromamba
pip install ./python && R CMD INSTALL .
bash tools/smoke_test.sh                                            # end-to-end check on simulated data
```

To reproduce results exactly, use the container of the release: per-cell
calls are identical across installations, but the chromatin-kNN pooled
calls depend slightly on the BLAS and `irlba` versions.

From source into existing installations: `pip install ./python` (Python
\>= 3.10) and, in R (\>= 4.2),
`install.packages("pak"); pak::pak("ljwharbers/scphasedcn")`
(Bioconductor dependencies: rhdf5, Rsamtools). A `tabix` on the PATH
speeds up reading the fragment table.

## Quick start

``` bash
# 1. phased heterozygous SNPs (Clair3 + longphase on the pooled BAM; see vignette("scphasedcn-inputs"))
scphasedcn select-hets --vcf sample.longphase.vcf.gz --loci 1000G_snvs_hg38.vcf.gz --out het_phased.vcf.gz
#    normal tissue (recommended): join the blocks into chromosome-wide haplotypes with a population phase of the sample
#    (e.g. Beagle 5 + 1000 Genomes panel, one VCF per chromosome), then extract with het_hybrid.vcf.gz
scphasedcn hybrid-phase --blocks het_phased.vcf.gz --population 'popphased_{chrom}.vcf.gz' --out het_hybrid.vcf.gz

# 2. extraction (all chromosomes in parallel, then merge; tumour mode: --snps het_phased.vcf.gz)
scphasedcn run --bam sample.bam --barcodes barcodes.tsv.gz --fasta GRCh38.fa \
    --snps het_hybrid.vcf.gz --phased --out results/sample/extract --jobs 8

# 3. calling
CALL=$(Rscript -e 'cat(system.file("scripts", "call_scn.R", package = "scphasedcn"))')   # in the container: call_scn.R
Rscript $CALL --dir results/sample/extract --out results/sample/caller --sex female --mode normal --null-cells select
```

or in R:

``` r

library(scphasedcn)
res = call_scn("results/sample/extract", mode = "tumour", sex = "male")
write_scphasedcn(res, "results/sample/caller")
```

Documentation: <https://ljwharbers.github.io/scphasedcn/>, or in R
[`vignette("scphasedcn")`](https://ljwharbers.github.io/scphasedcn/articles/scphasedcn.md)
(overview and a run on simulated data),
[`vignette("scphasedcn-inputs")`](https://ljwharbers.github.io/scphasedcn/articles/scphasedcn-inputs.md)
(input requirements, phasing recipe, extraction, cell metadata),
[`vignette("scphasedcn-outputs")`](https://ljwharbers.github.io/scphasedcn/articles/scphasedcn-outputs.md)
(result tables).

Scope: GRCh38 / hg38 with `chr` names; Linux and macOS.

## Citation

Please cite the SPLONGGET paper (reference to be added on publication)
and the package: `citation("scphasedcn")`.

## License

MIT (`LICENSE.md`). The bundled ENCODE blacklist (Amemiya et al. 2019)
and centromere coordinates (Battenberg) are credited in
`python/src/scphasedcn/resources/README.md`.
