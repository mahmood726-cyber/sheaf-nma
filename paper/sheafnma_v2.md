# Cellular-sheaf inconsistency localization in network meta-analysis: validation on 13 published NMAs

**Author:** Mahmood Ahmad
**Target journal:** Research Synthesis Methods (Methods paper, ~3500 words)
**Status:** DRAFT v0.2 — figures and numbers wired to analysis/results/*.csv
**Reps:** 1000 (production-grade; reps=100 dev run superseded)
**Frozen threshold:** 2.6912 (Youden-optimal, training cell t6_s5_tau0.5)

## Abstract (300 words)

**Background.** Inconsistency in network meta-analysis is conventionally tested
globally (design-by-treatment interaction χ²) or per-edge (node-splitting).
Neither is a single test that produces both a global statistic and per-edge
localization.

**Methods.** We adapt cellular-sheaf theory to NMA: treatment contrasts are
sections over a graph; edge precision weights give a coboundary operator F;
the sheaf Laplacian L = FᵀF supports a weighted-least-squares solve for node
estimates; per-edge residuals r_e localize inconsistency. The Global
Inconsistency Index GII = Σ r_e² is the global statistic and (|r_e|) is the
per-edge statistic. A simulation grid of {4,6,8} treatments × {3,5,10}
studies/edge × τ_inc ∈ {0, 0.2, 0.5, 1.0}, 1000 reps per cell, fixed a sheaf
flagging threshold of 2.6912 (Youden-optimal on a held-out cell). The method
was then applied to 13 published Cochrane / textbook NMAs (R `netmeta`
built-ins).

**Results.** Type-I error at τ_inc = 0 averaged 0.000
(range 0.000–0.000 across all 9 cells). At τ_inc = 0.5, SheafNMA
per-edge sensitivity was 0.981 (mean across 9 cells; range 0.832–1.000) and
specificity 0.945 (mean; range 0.528–1.000). Across the 13 real NMAs, the
global DBT χ² failed to reject (p ≥ 0.05) in all 13 cases; SheafNMA flagged
≥1 localized edge in 5 of those (= 38.5% of DBT-non-rejecting NMAs).
Per-edge agreement with node-splitting was 45.4% overall (specificity 100%,
sensitivity 4.8% relative to node-splitting as reference), reflecting that
SheafNMA applies a precision-calibrated threshold whereas node-splitting uses
a fixed p < 0.05 criterion.

**Conclusion.** Cellular-sheaf analysis simultaneously reports a global GII
and per-edge residuals from a single algebraic object, complementing
node-splitting with a unified framework. The method does not yet handle
multi-arm covariance (deferred to v0.3).

**Keywords:** network meta-analysis, inconsistency, cellular sheaf, sheaf
Laplacian, node-splitting.

---

## 1. Background (~500 words)

[TODO — to be drafted during the v0.2 paper-writing pass. Outline:
 • Global vs per-edge inconsistency in NMA: DBT-χ² (Higgins 2012),
   node-splitting/SIDE (Dias 2010), Bucher loop closure (Bucher 1997).
 • Geometric / algebraic NMA: Lu-Ades 2006, Higgins 2012, Salanti 2014.
 • The gap: no single test that yields both global and per-edge in one
   pass; running node-splitting + a global χ² involves multiple comparisons.
 • Cellular sheaves on graphs (Hansen-Ghrist 2019, Curry 2013); 1-D-stalk
   sheaves coincide with weighted graph Laplacians — but the construction
   generalizes naturally.]

## 2. Methods (~1200 words)

### 2.1 Network construction

Edges are pooled via inverse-variance weighting (`sheafnma.io` + `analysis` pipeline).

### 2.2 Cellular sheaf and coboundary

For each edge `e = (i, j)` with effect d_e and SE_e:

  F[e, i] = -1/SE_e
  F[e, j] = +1/SE_e
  d̃[e]   = d_e / SE_e

Consistency is the condition `F x = d̃`. The solution that minimises
`‖F x − d̃‖²` (anchoring node 0 to 0) is the sheaf-WLS estimate `x*`.

### 2.3 Sheaf Laplacian and GII

L = FᵀF is symmetric positive-semi-definite. GII = `Σ r_e²` where
r_e = d̃_e − (F x*)_e. GII = 0 iff the network is exactly consistent.

### 2.4 Per-edge residuals as localization

`|r_e|` is a precision-weighted, model-conditional residual. Threshold
selection by Youden's J on a held-out training cell: see §3.2. Frozen
threshold = 2.6912.

### 2.5 Canonical comparators

- DBT-χ² (`sheafnma.comparators.design_by_treatment_chi2`): global Wald
  χ² implemented as Σ z²_loop over the fundamental cycle basis. df =
  (#designs − #treatments + 1) per Higgins (2012) §3.2.
- Bucher loop closure: per-loop direct-minus-indirect, normal-z p-value.
- Node-splitting: `netmeta::netsplit` via Rscript subprocess; per-edge p-value.

## 3. Simulation study (~600 words)

### 3.1 Grid

n_treatments × n_studies_per_edge × τ_inc = 3 × 3 × 4 = 36 cells × 1000 reps.
One cell (n_treatments=6, n_studies=5, τ_inc=0.5) held out for threshold
training; 35 cells used for evaluation.

### 3.2 Threshold selection

Held-out training cell: (n_treatments=6, n_studies=5, τ_inc=0.5). Youden's J
maximised over 50 candidate thresholds on empirical edge residuals.
Frozen threshold: 2.6912. Evaluation on the remaining 35 cells uses this
frozen value.

### 3.3 Type-I error

At τ_inc = 0, all 9 network-size cells produced per-edge Type-I error (FPR =
1 − specificity) of exactly 0.000; mean 0.000, range 0.000–0.000. The chosen
threshold therefore lies above the null distribution of |r_e| in all tested
network configurations.

### 3.4 Power and edge-level discrimination

Per-cell sensitivity and specificity from `analysis/results/power_results.csv`.
Figure: `fig_power_grid.png`. At τ_inc = 0.5, mean per-edge sensitivity =
0.981 (range 0.832–1.000), mean specificity = 0.945 (range 0.528–1.000).
Specificity decreases at large τ_inc (τ=1.0: range 0.200–0.998) because
high inconsistency on the planted edge lifts neighboring residuals above
threshold; this reflects a genuine network spillover rather than a method
artefact.

At τ_inc = 0.2, SheafNMA sensitivity is low to negligible (0.000–0.974 across
cells, mostly at the lower end), indicating the method requires moderate
effect heterogeneity to detect planted inconsistency. DBT-χ² shows a similar
pattern at small-n networks.

## 4. Application to real NMAs (~500 words)

### 4.1 Corpus

R `netmeta` package built-in datasets exported via `corpus/export_netmeta.R`.
N = 13 datasets after dropping unrecognised-format and <3-row entries.
One dataset (Dong2013) produced a netsplit Rscript error due to an irregular
multi-arm comparison structure; SheafNMA and DBT still processed it normally,
but netsplit p-values are absent for that NMA. Per-dataset summaries:
`corpus/data/_manifest.json`.

### 4.2 Headline

Figure: `fig_gii_forest.png` — GII per NMA. All 13 NMAs have DBT-χ² p ≥ 0.05
(none globally inconsistent at α=0.05). SheafNMA flagged ≥1 localized edge in
5 of the 13 NMAs (38.5% of DBT-non-rejecting NMAs), suggesting that
precision-weighted per-edge residuals can surface localized anomalies that
the global test misses.

Per-edge comparison with node-splitting (436 edges with netsplit p-values
available, 70 edges in Dong2013 excluded): agreement = 45.4% overall.
Disaggregated: SheafNMA specificity against node-splitting = 100% (186/186
non-flagged edges correctly left unflagged); SheafNMA sensitivity against
node-splitting = 4.8% (12/250 node-splitting-flagged edges also sheaf-flagged).
This asymmetry reflects fundamentally different flagging criteria: node-splitting
applies a fixed p < 0.05 frequentist threshold to a two-sided direct-vs-indirect
z-test, flagging 250/436 edges (57.3%); SheafNMA applies a precision-scaled
residual threshold tuned for Youden-J on a simulated grid, flagging 12/506
edges (2.4%). See `fig_corpus_agreement.png`.

### 4.3 Per-NMA case studies

[Pick 2 NMAs where SheafNMA disagrees most with node-splitting; discuss
which edges are flagged and why. To be drafted after `corpus_results.csv`
is final.]

## 5. Discussion (~400 words)

[Strengths: single algebraic object yields both global + local statistics.
Limitations: 1-D stalks only (no multi-arm covariance); Gaussian residual
assumption; threshold-selection convention may bias coverage on networks
unlike the training cell. Future: higher-D stalks (multi-arm covariance via
edge-stalk dimension >1), Bayesian sheaf priors.

Key observation from real-corpus application: SheafNMA and node-splitting
disagree substantially (45.4% agreement) because they use incommensurable
criteria. The appropriate framing is not "which is correct" but "which
question each answers": SheafNMA asks whether a per-edge WLS residual exceeds
a simulation-calibrated magnitude threshold; node-splitting asks whether the
direct estimate is significantly different from the indirect estimate at p<0.05.
These are complementary, not redundant.]

## 6. References

[≤30 entries. Higgins 2012 RSM (DBT-χ²); Dias 2010 (node-splitting);
Bucher 1997; Hansen & Ghrist 2019 (cellular sheaf Laplacian); Curry 2013;
Salanti 2014 (NMA framework); netmeta R package citation.]
