# Plot the tracks of one cell

Three stacked tracks on a concatenated genome axis:

- depth log2 ratio (grey: 1 Mb windows; black: `agg` windows) with the
  expected log2 of the fitted segments;

- oriented B-allele fraction (`agg` windows, point size by the number of
  molecule votes) with the expected major / minor fractions (none on
  haploid chromosomes);

- a ribbon of the fitted nA:nB states (labelled when longer than 20 Mb),
  with segments confirmed in the cell outlined.

## Usage

``` r
plot_cell_tracks(
  tracks,
  cell,
  ylim_log2 = c(-2, 2),
  min_hap = 3,
  raster = TRUE,
  base_size = 11
)
```

## Arguments

- tracks:

  [`cell_tracks()`](https://ljwharbers.github.io/scphasedcn/reference/cell_tracks.md)
  output

- cell:

  cell barcode

- ylim_log2:

  limits of the log2 track (points outside are squished onto the border)

- min_hap:

  minimal number of molecule votes of a BAF point

- raster:

  rasterise the 1 Mb points (ggrastr)

- base_size:

  base font size

## Value

a patchwork object

## Details

Segments are coloured by state label; inherited segments are dashed,
confirmed and clone (neutral) segments solid. Needs ggplot2 and
patchwork (ggrastr, when installed, rasterises the 1 Mb points; when it
cannot rasterise here, e.g. without a working Cairo, the points are
drawn unrasterised, with one message per session).
