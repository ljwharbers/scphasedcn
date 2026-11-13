# Caller parameters from the caller config

[`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
updated from a YAML file (default `inst/config/caller.yaml` of the
package: one key per default, named vectors as maps) and then from
`overrides`. Values are coerced to the type of the default (numeric,
logical, character; names kept for named vectors). Empty strings and
nulls keep the default (e.g. `tabix: ""` = the `tabix` on the PATH); an
empty list (`[]`) is an empty vector (e.g. `allelic_exclude: []`).
Unknown keys in the file are ignored with a warning; unknown keys in
`overrides` are added as-is.

## Usage

``` r
caller_params(file = NULL, overrides = list())
```

## Arguments

- file:

  YAML file (NULL = the package's `caller.yaml`; a missing file gives
  the defaults)

- overrides:

  named list of values that replace the file's / default values

## Value

named list (the keys of
[`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md),
plus `tumour` / `normal` when the file has those blocks)

## Details

Two optional blocks hold the arguments of the two callers: `tumour:`
(tumour-mode options, any
[`clone_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/clone_defaults.md)
key) and `normal:` (arguments of
[`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
other than the inputs
[`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
sets, e.g. `n_draws`, `pool_loo`). When a block is present the result
carries it as `$tumour` / `$normal`, filled up to the full defaults
([`clone_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/clone_defaults.md),
the defaults of
[`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md));
[`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
passes it to the caller of its mode (its `caller_args` take precedence)
and records the values used in the result's `params`, so the
`params.yaml` of
[`write_scphasedcn()`](https://ljwharbers.github.io/scphasedcn/reference/write_scphasedcn.md)
reads back with this function.
