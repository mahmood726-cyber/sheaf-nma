# SheafNMA: Sheaf-Theoretic NMA Consistency Analysis — Design Spec

**Date:** 2026-03-26
**Author:** Mahmood Ahmad
**Target:** Single-file HTML app, ~2.5K lines, `C:\Models\SheafNMA\sheaf_nma.html`
**Paper target:** JRSS Series A or Biometrics
**Phase:** 1 (Lean MVP — sheaf Laplacian + inconsistency heatmap)

---

## 1. Purpose

SheafNMA is the world's first browser-based tool for localizing inconsistency in network meta-analysis using cellular sheaf theory. Traditional methods (node-splitting, net heat plot) detect *whether* inconsistency exists but provide limited insight into *where* it concentrates. The sheaf Laplacian's eigendecomposition pinpoints which treatment comparisons drive inconsistency and decomposes it into independent modes.

**Problem it solves:** After running an NMA, the reviewer asks "is there inconsistency, and if so, which comparisons are driving it?" Current tools give a global test (Q_inconsistency) or edge-by-edge node-splitting (multiple tests, no decomposition). Sheaf theory provides a principled mathematical framework that answers both questions simultaneously.

---

## 2. Architecture

Single-file HTML app, all JS + CSS embedded, no external dependencies. Sheaf-only diagnostic tool — users bring contrast-level NMA data, we localize inconsistency.

### Pipeline

```
Contrast data (treat1, treat2, effect, SE, study_id)
  → Build treatment network graph
  → Pool multi-study edges (inverse-variance)
  → Construct cellular sheaf (stalks + restriction maps)
  → Coboundary map F (precision-weighted)
  → Sheaf Laplacian L = F^T F
  → Solve WLS: x* = argmin ||Fx - d||²
  → Edge residuals: r_e = d_e/SE_e - (x*_i - x*_j)/SE_e
  → Edge inconsistency scores: |r_e|²
  → Eigendecomposition of L
  → Render: network heatmap + spectral view
```

### 4 Tabs

1. **Data** — demo datasets, CSV import, network preview
2. **Network** — interactive graph with edges colored by inconsistency (hero)
3. **Spectrum** — eigenvalue scree plot + eigenvector loadings (Advanced)
4. **Export** — R code, CSV, PNG, JSON

---

## 3. Sheaf Engine

### 3.1 Network Construction

- Parse contrast data: nodes = unique treatments, edges = unique (treat1, treat2) pairs
- Orient edges consistently: alphabetical order (treat1 < treat2)
- Multi-study edges: pool via inverse-variance weighting
  - d_pooled = sum(w_i * d_i) / sum(w_i), where w_i = 1/SE_i²
  - SE_pooled = 1/sqrt(sum(w_i))
- Store: node list, edge list with pooled effect and SE, per-edge study count

### 3.2 Cellular Sheaf

- **Node stalks:** R¹ per node (the "true" treatment effect relative to reference)
- **Edge stalks:** R¹ per edge (the direct comparison estimate)
- **Restriction maps:** For edge e = (i,j): restriction from node i is +1, from node j is -1
- **Sheaf condition:** d_e should equal effect_i - effect_j if the network is consistent

### 3.3 Coboundary Map

- Matrix F with dimensions (n_edges × n_nodes)
- For edge e connecting nodes i and j: F[e, i] = +1/SE_e, F[e, j] = -1/SE_e
- Precision-weighted so edges with smaller SE have more influence
- Observed data vector: d = [d_e / SE_e] for each edge

### 3.4 Sheaf Laplacian

- L = F^T × F (n_nodes × n_nodes, symmetric positive semi-definite)
- This is the precision-weighted graph Laplacian derived from sheaf theory
- For NMA with 1-dimensional stalks, it coincides with the standard weighted Laplacian — but the sheaf derivation generalizes to higher-dimensional stalks (Phase 2)

### 3.5 Inconsistency Computation

- **Network estimates:** x* = (F^T F)^{-1} F^T d = L^{-1} F^T d (pseudoinverse for rank-deficient L)
  - Use: fix reference treatment (node 0) at x*_0 = 0, solve reduced (n-1) × (n-1) system
- **Edge residuals:** r_e = d_e/SE_e - F[e,:] × x*
- **Edge inconsistency score:** s_e = r_e² (squared residual)
- **Normalize:** s_e / max(s_e) × 100 for 0-100% color scale
- **Global inconsistency index:** GII = sum(non-zero eigenvalues of L) / (n_nodes - 1)

### 3.6 Eigendecomposition

- QR algorithm with Householder reflections for symmetric matrices
- For n ≤ 20 treatments (typical NMA), computation is instant
- λ₁ = 0 always (connected graph); number of zero eigenvalues = connected components
- Non-zero eigenvalues = independent inconsistency modes
- Eigenvector v_k of eigenvalue λ_k: large |v_k[i]| means node i participates in the k-th inconsistency mode

---

## 4. Demo Datasets

### 4.1 Smoking Cessation (24 trials, 4 treatments)

Treatments: A = No contact, B = Self-help, C = Individual counselling, D = Group counselling.

Contrast-level data from Hasselblad 1998 / Lu-Ades 2004 (log OR scale). 24 studies across 6 possible edges (not all observed directly).

Expected: B-vs-D edge flagged as highest inconsistency (known from published node-splitting analyses).

### 4.2 Simulated with Planted Inconsistency (5 treatments, 12 contrasts)

True effects: A=0, B=-0.3, C=-0.5, D=-0.8, E=-1.0. Generate 10 consistent contrasts (effect = true_i - true_j + N(0, 0.1²)), SE ~ U(0.1, 0.3). Plant 2 inconsistent contrasts on A-D edge: shift by +0.5. Ground truth known.

Expected: sheaf localizes inconsistency to A-D edge. Largest eigenvector loads on nodes A and D.

Validation: planted edge ranked #1 by inconsistency score.

---

## 5. UI Design

### Tab 1 — Data

- Dropdown: "Smoking Cessation", "Simulated (Planted)", "Custom"
- Custom: textarea for CSV paste + file upload
- Expected columns: study, treat1, treat2, effect, se
- Network preview: simple node-edge diagram with study counts per edge
- "Analyze Consistency" primary button

### Tab 2 — Network (hero)

- Force-directed graph on Canvas (Fruchterman-Reingold)
- Nodes: circles labeled with treatment names, size ∝ degree
- Edges: colored green (consistent) → yellow → red (inconsistent), thickness ∝ precision
- Hover edge: tooltip with direct effect [CI], inconsistency score, study count
- Hover node: tooltip with estimated effect [CI], degree
- Color legend: gradient bar from green (0%) to red (100%)
- Below graph: ranked edge table sorted by inconsistency score (most inconsistent first)

### Tab 3 — Spectrum (Advanced)

- Left: eigenvalue scree plot (bar chart). Zero eigenvalues in grey, non-zero in blue. Clickable.
- Right: eigenvector loadings for selected eigenvalue. Bar chart of |v_i| per node.
- GII displayed prominently: "Global Inconsistency Index = X.XX"
- Annotation: "K non-zero eigenvalues = K independent inconsistency modes"

### Tab 4 — Export

- R code: igraph + manual Laplacian + netmeta comparison
- CSV: edge inconsistency scores, node estimates
- PNG: network heatmap, scree plot (canvas-to-blob, revoke after use)
- JSON: full results

### Standard Features

- Dark mode (CSS custom properties)
- Responsive layout
- Keyboard accessible (ARIA tablist/tab/tabpanel)
- escapeHtml() on all user input
- Seeded PRNG (xoshiro128**)

---

## 6. Testing Strategy

### Unit Tests (in-browser, 8 categories)

1. Network construction: correct node/edge counts, multi-study pooling
2. Coboundary map: correct dimensions, correct signs (+1/-1)
3. Laplacian: symmetric, positive semi-definite (all eigenvalues ≥ 0)
4. Eigenvalues: at least one zero (connected graph), non-negative
5. Eigenvectors: orthogonal (v_i · v_j ≈ 0 for i ≠ j)
6. Consistent network: all residuals ≈ 0, all eigenvalues ≈ 0 (except structural)
7. Simulated: planted inconsistent edge ranked #1
8. Smoking cessation: produces valid output (no NaN, no crashes)

### R Validation

- Compare edge inconsistency scores against netmeta node-splitting results
- Compare GII against Q_inconsistency from netmeta
- Compare network estimates against netmeta pooled effects

---

## 7. Phase 2 Roadmap (not in MVP)

- Node-splitting comparison (side-by-side with sheaf)
- Permutation-based p-values per edge
- Heat diffusion animation (inconsistency propagation visualization)
- Higher-dimensional stalks (multi-outcome NMA)
- Sheaf cohomology (higher-order obstructions to consistency)
- Integration with Component NMA tool
- WebR validation tier

---

## 8. Success Criteria

1. Simulated dataset: planted inconsistent edge ranked #1 by sheaf score
2. Smoking cessation: valid output, no crashes, inconsistency scores plausible
3. Eigendecomposition correct: eigenvalues non-negative, eigenvectors orthogonal
4. All unit tests pass
5. Network graph renders interactively with <200ms on 4-treatment network
6. Single HTML file, no external dependencies, works offline
7. GII correlates with Q_inconsistency for smoking cessation
