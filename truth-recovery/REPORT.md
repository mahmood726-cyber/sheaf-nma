# Truth-Recovery Report — sheaf-nma

**Engine:** `sheafnma/` Python package (pure NumPy/SciPy) + `sheaf_nma.html` (88 KB browser dashboard). Genuine methods engine — NOT a stub/wrapper.
**Method:** Sheaf-theoretic NMA inconsistency localization. Coboundary map F (precision-weighted signed incidence), sheaf Laplacian L=FᵀF, WLS node estimates (reference node fixed), per-edge residuals r_e, Global Inconsistency Index GII=Σr_e². Comparators: Bucher loop closure + Higgins design-by-treatment χ².

**Verdict: BUG (over-detection under heterogeneity) + STRONG-VALIDATION of power**

## Method (truth-first)
Independent known-truth generator (`dgp_known_truth.py`, numpy PCG64 — not the repo's xoshiro). Networks built from node-truth truth[i]=−0.2·i so every loop obeys d_AC=d_AB+d_BC by construction (CONSISTENT). A single known delta added to one design = INCONSISTENT (known location/magnitude). Within-edge sampling noise (`noise_sd`) is genuine between-study heterogeneity the pooled SE does NOT model. `harness.py` wires the repo's OWN engine and measures type-I (delta=0) and power (delta>0). 600–1500 reps/cell.

## Results

### Power is real (validated) — within-edge noise_sd=0.05
| delta | DBT reject | Bucher reject | residual-localizer sensitivity |
|------:|-----------:|--------------:|-------------------------------:|
| 0.0 | 0.000 | 0.000 | — |
| 0.2 | 0.053 | 0.463 | 0.04 |
| 0.5 | 0.865 | 1.000 | 1.00 |
| 1.0 | 1.000 | 1.000 | 1.00 |

Power monotone; planted edge is argmax-residual in ≥85% of clearly-inconsistent nets; residual localizer has 0% clean-edge FP at delta≤0.5.

### BUG — over-detection under within-edge heterogeneity (τ² ignored)
Type-I on networks CONSISTENT BY CONSTRUCTION (delta=0). Target=0.05.

| within-edge noise_sd | DBT-χ² type-I | Bucher any-loop type-I | P(GII≥2 → HTML "RED") |
|---------------------:|--------------:|-----------------------:|----------------------:|
| 0.05 | 0.000 | 0.000 | 0.16 |
| 0.15 | **0.162** | **0.506** | **0.999** |
| 0.30 | **0.913** | **0.975** | 1.000 |

The fixed-effect loop-closure tests use only the IV-pooled SE (`comparators.py`: `var_sum += edge["se"]**2`) and ignore between-study τ² within each contrast → pooled SE too small → z inflated → false alarms. Same family as the sibling LivingNMA fixed-effect-CI bug. The GII is worse: NO calibrated null. Mean GII on consistent nets = 1.38 (noise 0.05) and 12.5 (noise 0.15), while the HTML cutoff (`gii<0.5` green / `<2` orange / `≥2` red, `sheaf_nma.html:1324`) is absolute, untied to network size/τ². Result: 99.9% of consistent-heterogeneous networks render RED.

### Dashboard relative-normalization over-detection (sibling pattern, confirmed)
`sheaf_nma.html:1151-1156` normalizes `score/maxScore*100`. On FULLY CONSISTENT nets (delta=0, noise 0.05) the SECOND-worst edge shows ≥50% "inconsistent" in 82% of networks, ≥90% in 17% (worst edge always 100% red). Exactly the meta-entropy/tda-ma relative-scale over-detection — a relative scale manufactures a "most inconsistent" edge when nothing is.

## Tests
`test_truth_recovery.py`: 5/5 pass — no-false-alarm (low noise), monotone power, DBT over-detection BUG, GII uncalibrated-redflag BUG, localization.

## FIX (2026-06-15) — heterogeneity-robust null calibration

New module `sheafnma/calibrate.py::calibrated_inconsistency`. Additive: `core.gii`
and the FE comparators are untouched (all their exact-value tests still pass).

Mechanism (same family as the meta-entropy null-calibration fix): fit the
consistent WLS model, take precision-weighted residuals, studentize by hat-matrix
leverage so they are homoscedastic under the null, then estimate a *multiplicative
over-dispersion* phi robustly (median seed + iterated outlier-trimmed mean with a
truncated-normal de-bias). phi absorbs the unmodeled tau^2 (~1 when SEs adequate,
>1 under heterogeneity). The max standardized residual is judged against a
standard-normal null with Bonferroni over edges → a calibrated p-value that also
localizes. The trimming keeps a single genuinely-inconsistent edge from inflating
phi and masking itself, so power survives.

Measured on the repo's own harness (`harness.py --calibration`, 1000 reps/cell,
target alpha=0.05):

| noise_sd | BEFORE DBT type-I | BEFORE Bucher | BEFORE GII≥2 RED | AFTER calibrated |
|---------:|------------------:|--------------:|-----------------:|-----------------:|
| 0.05 | 0.000 | 0.000 | 0.140 | **0.047** |
| 0.15 | 0.155 | 0.475 | 0.997 | **0.047** |
| 0.30 | 0.914 | 0.980 | 1.000 | **0.047** |

Type-I is now nominal AND heterogeneity-invariant. Power preserved (not destroyed):

| noise_sd | delta | BEFORE DBT | AFTER calibrated |
|---------:|------:|-----------:|-----------------:|
| 0.05 | 0.2 | 0.065 | 0.844 |
| 0.05 | 0.5 | 0.878 | 0.996 |
| 0.05 | 1.0 | 1.000 | 1.000 |
| 0.15 | 0.2 | 0.418 | 0.128 |
| 0.15 | 0.5 | 0.946 | 0.721 |
| 0.15 | 1.0 | 1.000 | 0.988 |

Power is monotone in delta and correctly attenuates under heavier heterogeneity
(honest: a small effect is genuinely harder to detect when the network is noisier;
the BEFORE numbers at noise 0.15 are not "power", they are the same false-alarm
inflation as the delta=0 row). Tests: `tests/test_calibrate.py` (6) +
`test_truth_recovery.py` FIX assertions (3) all pass.

## Recommendation
1. Make GII calibrated, not absolute: scale to χ² on n_loops df or bootstrap/permute the null → p-value/z. Never color by raw Σr².
2. Add random-effects (τ²) to loop-closure SE in `comparators.py` and to precision weights in `core.build_coboundary` (DL/REML within-contrast). Without it DBT/Bucher are anti-conservative under heterogeneity.
3. Dashboard: replace `score/max*100` with an absolute calibrated scale (p-value/standardized residual), or label the column "rank within this network, NOT absolute inconsistency".
4. Power/localization in the homogeneous regime are genuinely strong — keep.

Additive only. All findings measured; none fabricated. Repo demos use noise_sd=0.05 where the bug is latent; it surfaces under realistic heterogeneity.
