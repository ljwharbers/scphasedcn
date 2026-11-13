# Inputs: from a SPLONGGET BAM to the extraction directory

## What you need

| input | required | notes |
|----|----|----|
| DNA BAM | yes | coordinate sorted and indexed, GRCh38 / hg38 with `chr` names, the cell barcode in the `CB` tag; duplicates flagged or removed |
| cell whitelist | yes | one barcode per line (a trailing `-1` is stripped) |
| reference FASTA | yes | the FASTA the reads were aligned to (GC content per bin), with its `.fai` |
| phased heterozygous SNPs | for the allele channel; required in tumour mode | VCF with `GT` phased (`0|1`) and a `PS` phase-set tag, e.g. Clair3 + longphase, selected with `scphasedcn select-hets` |
| sex | yes | `female` or `male` (chrX / chrY copy number) |
| cell metadata | no | cell types (normal tissue), cells known to be copy-number neutral, groups; see below |
| known events | no | regions reported as their own calls (e.g. a known trisomy; `name=chr:start-end`, 0-based half-open, an open end = the chromosome end); normal mode calls them per cell, in pools only with `pool_events` (off by default, experimental) |

The genome is hg38 only: the bundled chromosome sizes, ENCODE blacklist
v2 and centromere coordinates are hg38, and so are the chrX centromere
and acrocentric arm defaults of the caller.

Reads: the extractor keeps primary, non-duplicate records with MAPQ \>=
30 whose `CB` is in the whitelist; each record is one fragment.
SPLONGGET read names of the form `<CB>_#<uuid>_<sign><piece>of<npieces>`
let it recognise pieces of one raw read; other read names are accepted
(each record is then its own read).

Whitelists in RNA-barcode space (10x Multiome GEX barcodes) are
translated to the ATAC-space barcodes of the BAM with `--barcode-map` (a
TSV with columns `atac_barcode` and `rna_barcode`).

## 1. Phasing

The allele channel needs heterozygous SNPs phased within the sample.
Call and phase them on the pooled reads of all cells (the “pseudobulk”
BAM); one recipe that works for ONT R10.4.1 data:

``` bash
# 1. a clean BAM: primary, mapped, not duplicate, MAPQ >= 30
samtools view -@ 16 -b -F 0xd04 -q 30 -o sample.clean.bam sample.bam $(echo chr{1..22} chrX)
samtools index -@ 16 sample.clean.bam

# 2. Clair3 (docker://hkubal/clair3; pick the model of your basecaller, e.g. r1041_e82_400bps_sup_v520)
run_clair3.sh --bam_fn=sample.clean.bam --ref_fn=GRCh38.fa --threads=32 --platform=ont \
    --model_path=/models/r1041_e82_400bps_sup_v520 --ctg_name=$(echo chr{1..22} chrX | tr ' ' ',') \
    --sample_name=sample --output=clair3/

# 3. longphase on the PASS heterozygous biallelic SNVs
bcftools view -f PASS -v snps -m2 -M2 -g het -Oz -o pass_het.vcf.gz clair3/merge_output.vcf.gz
tabix -p vcf pass_het.vcf.gz
longphase phase -s pass_het.vcf.gz -b sample.clean.bam -r GRCh38.fa -o sample.longphase --ont -t 32
bgzip sample.longphase.vcf && tabix -p vcf sample.longphase.vcf.gz

# 4. select the sites the extractor uses
scphasedcn select-hets --vcf sample.longphase.vcf.gz --out het_phased.vcf.gz \
    --loci 1000G_snvs_hg38.vcf.gz --stats het_phased.select_stats.tsv
```

`select-hets` keeps PASS biallelic heterozygous SNVs with GQ \>= 10 and
DP between 20 and the 99th percentile of the candidates (collapsed
repeats have excess depth). There is deliberately no allele-frequency
window: it would remove sites in regions of LOH or amplification in
tumours. `--loci` (recommended) intersects the sites with a population
SNV panel, e.g. the 1000 Genomes sites on GRCh38 (position, REF and ALT
must match), which removes most systematic long-read genotyping errors;
chrX is kept without the intersect by default (`--no-loci-contigs`).

longphase phase blocks are not oriented relative to each other. The
caller handles that: normal-tissue mode orients each block of a cell on
the pseudobulk of the other cells (leave-one-out), tumour mode per
clone. Population (reference-panel) phasing that joins the blocks into
chromosome-wide haplotypes is not needed. It is used when present (one
phase set per chromosome), but statistical phasing flips a fraction of
the block junctions, which tumour mode then has to undo.

## 2. Extraction (`scphasedcn`)

``` bash
scphasedcn run --bam sample.bam --barcodes barcodes.tsv.gz --fasta GRCh38.fa \
    --snps het_phased.vcf.gz --phased --out results/sample/extract --jobs 8
```

`run` extracts every chromosome in parallel (`--jobs`; a deep chromosome
takes a few GB of memory) and merges. On a cluster, run
`scphasedcn extract --chrom <chrom>` as one job per chromosome and
`scphasedcn merge` after them, with the same options. A sample yaml
keeps the paths in one place (`--config sample.yaml --sample <key>`,
also read by the caller’s command line):

``` r
cat(readLines(system.file("extdata", "sample_template.yaml", package = "scphasedcn")), sep = "\n")
# Sample yaml shared by the extractor (`scphasedcn ... --config <yaml> --sample <key>`) and the caller
# (`call_scn.R --config <yaml> --sample <key>`). Copy it, fill in your paths, delete what you do not use.
# Paths are absolute. The keys shown are the documented ones; unknown keys are ignored.
reference:
  fasta: /path/to/GRCh38/genome.fa                 # required: hg38 with "chr" names, indexed (.fai); gives the GC per bin
  # centromeres: /path/to/centromeres.tsv          # optional: default = the hg38 table shipped with both packages
  # blacklist_v2: /path/to/hg38-blacklist.v2.bed   # optional: default = the ENCODE blacklist v2 shipped with the Python extractor
  # trf: /path/to/simple_repeats.trf.bed           # optional: tandem repeats (fraction per bin; omitted = no TRF mask)
outputs:
  root: /path/to/results                           # per-sample outputs: <root>/<sample>/extract, <root>/<sample>/caller/<tag>
samples:
  my_sample:
    type: normal                                   # normal | tumour (or tumor; the caller's default --mode)
    sex: female                                    # female | male
    bam: /path/to/my_sample.dna.bam                # SPLONGGET DNA BAM (coordinate sorted, indexed, CB tag per read)
    barcodes: /path/to/cell_barcodes.tsv.gz        # cell whitelist, one barcode per line (a trailing -1 is stripped)
    # barcode_map: /path/to/barcode_mapping.tsv    # only if the whitelist is in RNA-barcode space: columns
    #                                              # atac_barcode, rna_barcode (translated to the BAM's ATAC space)
    phased_vcf: /path/to/het_phased.vcf.gz         # phased heterozygous SNPs (GT with |, PS tag), e.g. from
    #                                              # Clair3 + longphase + `scphasedcn select-hets`
    # chroms: [chr1, chr2, ..., chr22, chrX]       # default: chr1-22, chrX, chrY, chrM (`run`: those in the BAM)
    # peaks: /path/to/peaks.bed                    # optional: ATAC peaks (FRiP in the QC table only)
    events:                                        # optional: known events, reported as their own regions
      # chr5_gain: "chr5"
      # 20q_gain: "chr20:28000000-"
```

Reference files: the FASTA is required. The blacklist, centromeres and
chromosome sizes default to the bundled hg38 tables; a tandem-repeat BED
(`--trf`) is optional (it adds the repeat fraction of each bin, used by
the bin filter).

## 3. The extraction directory

`call_scn(dir)` reads three files of the directory:

- `bins_100kb.h5ad` (AnnData): `X` = clean fragments per cell (rows) x
  100 kb bin (columns). `var` per bin: `chrom`, `start`, `end` (0-based,
  half-open), `gc`, `frac_blacklist`, `frac_trf`, `is_centromere`,
  `peak_bp`. `obs` per cell: the QC counts below.
- `fragments.tsv.gz` + `.tbi` (bgzipped, tabix-indexed, one header line
  starting with `#`): one row per fragment with columns
  `chrom, start, end, cell, strand, mapq, flags, len, piece, npieces, hap, nsnp, ps, read_id`.
  `hap` is the fragment’s haplotype from the phased SNPs it covers (0
  none, 1 haplotype A, 2 haplotype B, 3 ambiguous) and `ps` the phase
  set of its first phased SNP (-1 none). The caller collapses a cell’s
  fragments closer than `mol_gap` (2 kb) into one molecule with one
  vote.
- `cells_qc.tsv`: per cell `n_pass` (passing records), `n_clean`
  (outside the blacklist), `frac_mapq60`, `frac_chrM`, … The caller’s
  cell QC
  ([`cell_qc()`](https://ljwharbers.github.io/scphasedcn/reference/cell_qc.md))
  keeps cells with at least 5,000 clean fragments, at least half of them
  at MAPQ 60 and at most 20 % on chrM.

Also written: `alleles_dna.h5ad` (allele counts per cell x SNP),
`summary.json`, and `partials/` (the per-chromosome pieces).

The package’s simulator writes the same layout:

``` r

library(scphasedcn)
sim = simulate_scn(n_cells = 40, seed = 1)
d = file.path(tempdir(), "sim_extract")
invisible(write_sim_dataset(sim, d))
b = read_bins_h5ad(file.path(d, "bins_100kb.h5ad"))
dim(b$counts)
#> [1] 2100   40
head(b$bins)
#>                   bin  chrom start   end    gc frac_blacklist    frac_trf
#>                <char> <char> <num> <num> <num>          <num>       <num>
#> 1:      chr1:0-100000   chr1 0e+00 1e+05  0.41              0 0.022501392
#> 2: chr1:100000-200000   chr1 1e+05 2e+05  0.41              0 0.013028650
#> 3: chr1:200000-300000   chr1 2e+05 3e+05  0.41              0 0.035913346
#> 4: chr1:300000-400000   chr1 3e+05 4e+05  0.41              0 0.020442569
#> 5: chr1:400000-500000   chr1 4e+05 5e+05  0.41              0 0.008772530
#> 6: chr1:500000-600000   chr1 5e+05 6e+05  0.41              0 0.007845243
#>    is_centromere    arm n_peaks peak_bp bin_id
#>           <lgcl> <char>   <num>   <num>  <int>
#> 1:         FALSE      p       2    1000      1
#> 2:         FALSE      p       5    2500      2
#> 3:         FALSE      p       8    4000      3
#> 4:         FALSE      p       2    1000      4
#> 5:         FALSE      p       2    1000      5
#> 6:         FALSE      p       2    1000      6
head(read_fragments(file.path(d, "fragments.tsv.gz")))
#>     chrom start   end      cell strand  mapq flags   len piece npieces    hap
#>    <char> <num> <num>    <char> <char> <int> <int> <int> <int>   <int> <char>
#> 1:   chr1  2492  2716 cell_0012      -    60     2   224     1       1      U
#> 2:   chr1  4419  4866 cell_0003      -    60     2   447     1       1      B
#> 3:   chr1  4461  4937 cell_0004      +    60     2   476     1       1      B
#> 4:   chr1  4473  4895 cell_0037      -    60     2   422     1       1      B
#> 5:   chr1  4517  5391 cell_0007      -    60     2   874     1       1      B
#> 6:   chr1  4558  5720 cell_0034      -    60     2  1162     1       1      A
#>     nsnp    ps read_id
#>    <int> <int>   <int>
#> 1:     0    -1   41642
#> 2:     1  1005    6369
#> 3:     1  1005   12005
#> 4:     1  1005  140313
#> 5:     1  1005   24041
#> 6:     1  1005  131030
```

## 4. Cell metadata

Everything here is optional.
[`make_cell_meta()`](https://ljwharbers.github.io/scphasedcn/reference/make_cell_meta.md)
builds the table (`--cell-meta` of the command line) from the
extractor’s `cells_qc.tsv` and your labels:

- `cell_type`: cell types from an RNA (or ATAC) annotation of the same
  barcodes. Normal-tissue mode builds the depth reference of a cell from
  the other cells of its type, which removes cell-type chromatin effects
  on read depth; without cell types, all cells form one reference.
- `flat`: cells known to be copy-number neutral (e.g. from a matched
  assay). Normal-tissue mode uses them as its null cells with
  `null_cells = "features"`; without them use `null_cells = "select"`
  (the caller selects null cells from the data). Tumour mode fits its
  input depth and allelic null on them; its normal reference pool is the
  normal clone the caller finds and cleans itself, not these cells.
- `group`: tumour mode, groups whose pseudobulk orients the phase blocks
  (`groups = "group"`), e.g. spatial regions or RNA clusters.

``` r

qc = data.table::fread(file.path(d, "cells_qc.tsv"))
types = data.table::data.table(cell = qc$cell, cell_type = sim$scn$cells$cell_type)
m = make_cell_meta(file.path(d, "cells_qc.tsv"), cell_types = types)
head(m)
#>         cell qc_pass cell_type   flat  group
#>       <char>  <lgcl>    <char> <lgcl> <char>
#> 1: cell_0001   FALSE    tumour     NA   <NA>
#> 2: cell_0002   FALSE    normal     NA   <NA>
#> 3: cell_0003   FALSE    normal     NA   <NA>
#> 4: cell_0004    TRUE    tumour     NA   <NA>
#> 5: cell_0005   FALSE    normal     NA   <NA>
#> 6: cell_0006   FALSE    normal     NA   <NA>
data.table::fwrite(m, file.path(tempdir(), "cell_meta.tsv"), sep = "\t")
```

Multiplets (barcodes holding more than one nucleus) inflate calls; pass
a per-cell multiplet probability as `multiplets` (a table `cell, prob`)
to leave them out of the null and the clone discovery.
