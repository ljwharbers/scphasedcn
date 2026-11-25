# Write the tracks of many cells to a multi-page PDF

`per_page` cells per page
([`plot_cell_tracks()`](https://ljwharbers.github.io/scphasedcn/reference/plot_cell_tracks.md)),
with `grDevices::cairo_pdf(onefile = TRUE)`
([`grDevices::pdf()`](https://rdrr.io/r/grDevices/pdf.html) when R has
no cairo). The cells are checked before the file is opened; only `file`
(and its directory) is written, and only the device opened here is
closed.

## Usage

``` r
plot_cells_pdf(
  tracks,
  file,
  cells = tracks$meta$cell,
  per_page = 4,
  width = 14,
  height = 4 * per_page,
  ...
)
```

## Arguments

- tracks:

  [`cell_tracks()`](https://ljwharbers.github.io/scphasedcn/reference/cell_tracks.md)
  output

- file:

  PDF path

- cells:

  cells (default: all cells of the tracks)

- per_page:

  cells per page

- width, height:

  page size (inches)

- ...:

  passed to
  [`plot_cell_tracks()`](https://ljwharbers.github.io/scphasedcn/reference/plot_cell_tracks.md)

## Value

`file`, invisibly
