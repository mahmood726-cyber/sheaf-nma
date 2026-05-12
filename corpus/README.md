# SheafNMA real-NMA corpus

Source: R `netmeta` package built-in datasets, exported via
`corpus/export_netmeta.R`. Each CSV uses the sheafnma standard contrast
schema: `study, treat1, treat2, effect, se`.

See `data/_manifest.json` for per-dataset row counts and the netmeta
package version used at export time.

## Reproducing

```sh
Rscript corpus/export_netmeta.R
```

This wipes / rewrites `data/*.csv` and `data/_manifest.json`.

## Excluded datasets

Datasets the script skips:
- Format unrecognised (neither contrast-level `TE/seTE` nor arm-based
  `treatment/event/n` or `treatment/mean/sd/n`).
- Fewer than 3 usable contrast rows after dropping rows with non-finite
  `TE`/`seTE` or zero `se`.

Inspect the script stdout when re-running to see exact reasons.
