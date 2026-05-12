# SheafNMA

Cellular-sheaf inconsistency localization for network meta-analysis.

**Status:** v0.2.0 (shipped 2026-05-12). Real-NMA-corpus validation + simulation power study + methods paper draft.

## What it does

SheafNMA treats treatment contrasts as graph sections, builds a precision-weighted coboundary operator `F` and sheaf Laplacian `L = FᵀF`, then reports the **Global Inconsistency Index** (`GII = Σ r²`) alongside **per-edge residuals** `r_e` that localize inconsistency to specific treatment comparisons. One algebraic object yields both the global statistic and per-edge diagnostics — complementing node-splitting under a unified framework.

## v0.2 headline

> Across 13 published Cochrane / textbook NMAs (R `netmeta` built-ins spanning diabetes, smoking cessation, thrombolytics, acupuncture, antimanic agents), the global design-by-treatment χ² failed to reject in all 13. **SheafNMA flagged ≥1 localized inconsistency edge in 5 of those 13 (38.5%)** at the Youden-frozen threshold of 2.6912 — networks that look globally consistent but harbor at least one localized inconsistency the global test misses.

Simulation rigor (1000 reps × 36 cells):

- Type-I error at `τ_inc = 0`: **0.000** across all cells
- Edge-level sensitivity at `τ_inc = 0.5`: mean **0.981**, specificity mean **0.945**

## Repository layout

```
sheafnma/                   Python package
├── core.py                 coboundary, Laplacian, GII, per-edge residuals (pure NumPy)
├── io.py                   standard contrast-CSV loader (study, treat1, treat2, effect, se)
├── simulate.py             xoshiro128** RNG + planted-inconsistency network generator
├── comparators.py          DBT-χ² (Higgins 2012) + Bucher (1997) loop closure
└── r_bridge.py             Rscript subprocess wrapper for netmeta::netsplit + dataset loader

corpus/                     real-NMA corpus
├── export_netmeta.R        one-shot exporter (handles 7 netmeta dataset formats)
├── README.md
└── data/                   13 CSVs + _manifest.json

analysis/                   orchestrators (results CSVs gitignored)
├── run_real_corpus.py      apply sheaf + 3 comparators to each NMA
├── power_simulation.py     36-cell × 1000-rep grid with held-out Youden threshold
└── generate_figures_v2.py  paper figures

paper/
├── sheafnma_v2.md          ~3500-word RSM-format manuscript (1190 words drafted; [TODO] prose blocks)
├── methods_note.md         156-word / 7-sentence E156 body for #153
└── figures/                fig_power_grid.png, fig_corpus_agreement.png, fig_gii_forest.png

tests/                      33 tests, all green on master
docs/superpowers/
├── specs/2026-05-12-sheafnma-v0.2-real-corpus-validation-design.md
└── plans/2026-05-12-sheafnma-v0.2-implementation.md

generate_figures.py         Phase-1 figure script (now a thin wrapper around sheafnma.*)
sheaf_nma.html              Phase-1 interactive dashboard (unchanged in v0.2)
```

## Install + run

```bash
git clone https://github.com/mahmood726-cyber/sheaf-nma.git
cd sheaf-nma
python -m pip install -e .[test]
pytest tests/ -q                          # 33 passing

# Real-corpus run (requires R + netmeta)
python -m analysis.run_real_corpus
python -m analysis.power_simulation --reps 1000
python -m analysis.generate_figures_v2
```

R dependency: `netmeta` (tested against 3.2.0) at `C:\Program Files\R\R-4.5.2\bin\Rscript.exe` (or `Rscript` on PATH). Tests in `tests/test_r_bridge.py` auto-skip when R is unavailable.

## Methodology

| Statistic | What it tells you | Where in code |
|---|---|---|
| `gii(network)` | Global Inconsistency Index — sum of squared precision-weighted residuals | `sheafnma.core.gii` |
| `edge_residuals(network)` | Per-edge residuals `r_e = d̃_e − (F x*)_e` | `sheafnma.core.edge_residuals` |
| `design_by_treatment_chi2` | Higgins 2012 RSM §3 global Wald χ² | `sheafnma.comparators.design_by_treatment_chi2` |
| `bucher_loop_closure` | Bucher 1997 per-loop direct-minus-indirect | `sheafnma.comparators.bucher_loop_closure` |
| `run_netsplit` | Dias 2010 node-splitting via `netmeta::netsplit` | `sheafnma.r_bridge.run_netsplit` |

## Versions

- **v0.1.0** (2026-04-07): cellular-sheaf math + single-file HTML dashboard + E156 #153 micro-paper draft on a synthetic 12-contrast network.
- **v0.2.0** (2026-05-12): refactored into `sheafnma/` package, validated on 13 real netmeta NMAs, 36-cell × 1000-rep simulation, RSM-format methods paper draft. **E156 #153 CURRENT BODY replaced** with the v0.2 netmeta-corpus framing.

## Limitations

- Single-dimensional stalks only — no multi-arm covariance correction (deferred to v0.3).
- Threshold-selection convention (Youden's J on held-out cell) may bias coverage on networks unlike the training cell.
- Gaussian residual assumption.

## License

MIT. See `LICENSE`.

## Citation (preferred form, pending Crossref DOI)

> Ahmad, M. (2026). *Cellular-sheaf inconsistency localization in network meta-analysis: validation on 13 published NMAs* (v0.2.0). github.com/mahmood726-cyber/sheaf-nma · tag `v0.2.0`.
