# scphasedcn: single-cell copy-number calling from long-read single-cell DNA data

The caller
[`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
runs on an extraction directory written by the companion Python
extractor (`scphasedcn run`): per-cell region features
([`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)),
then a normal-tissue mode
([`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md))
or a tumour mode
([`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md):
clone discovery, clone states, cell re-assignment);
[`write_scphasedcn()`](https://ljwharbers.github.io/scphasedcn/reference/write_scphasedcn.md)
writes the result tables. Supporting functions: readers for the
extraction outputs
([`read_bins_h5ad()`](https://ljwharbers.github.io/scphasedcn/reference/read_bins_h5ad.md),
[`read_fragments()`](https://ljwharbers.github.io/scphasedcn/reference/read_fragments.md)),
bins and arms
([`annotate_arms()`](https://ljwharbers.github.io/scphasedcn/reference/annotate_arms.md),
[`bin_filter()`](https://ljwharbers.github.io/scphasedcn/reference/bin_filter.md)),
the depth references behind
[`normalise_depth()`](https://ljwharbers.github.io/scphasedcn/reference/normalise_depth.md),
beta-binomial allelic primitives, the MEDICC2 export
([`medicc2_pieces()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_pieces.md),
[`medicc2_table()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_table.md)),
plot helpers and a generative simulator
([`simulate_scn()`](https://ljwharbers.github.io/scphasedcn/reference/simulate_scn.md),
[`write_sim_dataset()`](https://ljwharbers.github.io/scphasedcn/reference/write_sim_dataset.md))
used by the examples and tests. Start with
[`vignette("scphasedcn")`](https://ljwharbers.github.io/scphasedcn/articles/scphasedcn.md).

## Details

Conventions: matrices are bins x cells (`dgCMatrix` for counts, dense
for ratios); coordinates are 0-based half-open as in BED / the Python
extraction layer; `=` assignment and data.table throughout.

## See also

Useful links:

- <https://github.com/ljwharbers/scphasedcn>

- <https://ljwharbers.github.io/scphasedcn/>

- Report bugs at <https://github.com/ljwharbers/scphasedcn/issues>

## Author

**Maintainer**: Luuk Harbers <luuk.harbers@vib.be>
