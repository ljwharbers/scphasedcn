# The regions table of the caller

Arms
([`arm_table()`](https://ljwharbers.github.io/scphasedcn/reference/arm_table.md)),
whole chromosomes (first to last bin), tiles of `params$tile_sizes` (one
per
[`window_ids()`](https://ljwharbers.github.io/scphasedcn/reference/window_ids.md)
window, named by
[`tile_name()`](https://ljwharbers.github.io/scphasedcn/reference/tile_name.md),
`end = min(start + size, chromosome end)`) and named events, on the
chromosomes that have usable bins (chrY only for males). Regions without
a usable bin are dropped. Columns follow
[`vignette("scphasedcn-outputs", package = "scphasedcn")`](https://ljwharbers.github.io/scphasedcn/articles/scphasedcn-outputs.md):
`n_bins` / `n_use` count all / usable bins under the containment rule;
tiles are `eligible` when `n_bins >= 0.9 * size / bin width` and
`n_use >= 0.7 * n_bins` (the eligibility rule of the spike-in
calibration), other regions when they have a usable bin; `null_region` =
arm, chrom and tile regions on autosomes not in `null_exclude`
(chromosome or region names); `allelic_ok` is FALSE on male chrX / chrY;
`sd_class` by
[`sd_class_of()`](https://ljwharbers.github.io/scphasedcn/reference/sd_class_of.md).

## Usage

``` r
region_set(
  bins,
  centromeres,
  classes = c("arm", "chrom", "5Mb", "10Mb", "20Mb", "50Mb"),
  events = NULL,
  sex = "female",
  null_exclude = character(),
  use = NULL,
  params = caller_defaults()
)
```

## Arguments

- bins:

  bins table (chrom, start, end)

- centromeres:

  Battenberg centromere table or path
  ([`arm_table()`](https://ljwharbers.github.io/scphasedcn/reference/arm_table.md))

- classes:

  region classes to build (`arm`, `chrom`, tile labels); `event` is
  added when `events` is given

- events:

  NULL, a data.table `region, chrom, start, end` (0-based half-open; NA
  end = chromosome end),
  [`cli_events()`](https://ljwharbers.github.io/scphasedcn/reference/cli_events.md)
  output or a named vector of region strings (`"chr20:28000000-"`,
  `"chr7"`;
  [`parse_region()`](https://ljwharbers.github.io/scphasedcn/reference/parse_region.md)).
  An event without a name is named after its region string (a bare
  chromosome `"chr7"`: `"event|chr7"`). An event that yields no region
  (a chromosome without usable bins, or no whole usable bin inside) is
  dropped with a warning.

- sex:

  "female" or "male" (case-insensitive, or "F" / "M"; any other value
  stops)

- null_exclude:

  chromosomes (or region names) whose regions never contribute null rows

- use:

  logical bin filter (NULL = all bins)

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

## Value

regions data.table (passes
[`check_regions()`](https://ljwharbers.github.io/scphasedcn/reference/check_tables.md))
