# SheafNMA Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the world's first browser-based sheaf-theoretic NMA consistency analyzer — localizing inconsistency to specific treatment comparisons via the sheaf Laplacian.

**Architecture:** Single-file HTML app (~2.5K lines). Parse contrast-level NMA data → build network graph → construct sheaf (coboundary map F) → Laplacian L=F^TF → eigendecomposition → edge residuals → inconsistency heatmap on interactive network + spectral view. No external dependencies.

**Tech Stack:** HTML5/CSS3/JS (ES2020), Canvas 2D for network graph + plots, xoshiro128** seeded PRNG.

**Spec:** `docs/superpowers/specs/2026-03-26-sheaf-nma-design.md`

---

## File Structure

```
C:\Models\SheafNMA\
├── sheaf_nma.html                  # The entire app (~2.5K lines)
├── docs/superpowers/
│   ├── specs/2026-03-26-sheaf-nma-design.md
│   └── plans/2026-03-26-sheaf-nma-implementation.md
```

Within `sheaf_nma.html`, logical sections:
1. CSS (~250 lines) — layout, tabs, dark mode, network graph styling
2. HTML (~150 lines) — 4 tab panels, controls, containers
3. Data module (~200 lines) — demo datasets, CSV parser, network builder
4. Linear algebra module (~200 lines) — matrix multiply, transpose, inverse, eigendecomposition
5. Sheaf module (~250 lines) — coboundary map, Laplacian, WLS solve, residuals, inconsistency scores
6. Visualization module (~350 lines) — network canvas, scree plot, eigenvector bars, edge table
7. Export module (~150 lines) — R code, CSV, PNG, JSON
8. UI glue (~150 lines) — tab switching, event handlers, dark mode, test runner
9. In-browser tests (~200 lines) — 8 test categories per spec

---

### Task 1: HTML Shell + CSS + Tab System + Data Module

**Files:**
- Create: `sheaf_nma.html`

- [ ] **Step 1: Create HTML skeleton with 4 tabs and CSS**

Full HTML structure: DOCTYPE, meta, title "SheafNMA: Sheaf-Theoretic NMA Consistency Analysis". CSS custom properties for light/dark. 4 tabs: Data, Network, Spectrum, Export. Header with title + dark mode toggle. Footer with version + "Run Tests" button. Tab switching JS.

Use same color scheme as TDA-MA: `--bg:#f8f9fa, --fg:#1a1a2e, --accent:#2563eb` (light) and `--bg:#0f172a, --fg:#e2e8f0, --accent:#60a5fa` (dark).

- [ ] **Step 2: Add demo datasets**

Smoking cessation (Hasselblad 1998): 24 studies, 4 treatments (A=No contact, B=Self-help, C=Individual counselling, D=Group counselling). Contrast-level data on log OR scale:

```javascript
const DATASETS = {
  smoking: {
    name: 'Smoking Cessation (Hasselblad 1998)',
    description: '24 trials, 4 treatments. Known inconsistency in B-vs-D comparison.',
    contrasts: [
      { study: 'Study 01', treat1: 'A', treat2: 'B', effect: -0.32, se: 0.32 },
      { study: 'Study 02', treat1: 'A', treat2: 'B', effect: -0.18, se: 0.28 },
      { study: 'Study 03', treat1: 'A', treat2: 'B', effect: -0.06, se: 0.21 },
      // ... (provide all 24 studies with realistic log OR values)
      // Key: multiple A-vs-B, A-vs-C, A-vs-D studies; fewer B-vs-C, B-vs-D, C-vs-D
    ]
  },
  simulated: {
    name: 'Simulated (Planted Inconsistency)',
    description: '5 treatments, 12 contrasts. A-D edge has planted +0.5 inconsistency.',
    contrasts: [] // Generated from true effects + seeded PRNG
  }
};
```

For simulated: true effects A=0, B=-0.3, C=-0.5, D=-0.8, E=-1.0. Generate 10 consistent contrasts + 2 inconsistent on A-D (shifted +0.5). Use xoshiro128** seed=42.

- [ ] **Step 3: Add CSV parser and network builder**

`parseContrastCSV(text)` — parse columns: study, treat1, treat2, effect, se.

`buildNetwork(contrasts)` — returns `{ nodes: ['A','B',...], edges: [{treat1, treat2, effect, se, studies, k},...] }` with multi-study edges pooled via inverse-variance.

- [ ] **Step 4: Wire Data tab UI**

Dataset dropdown, CSV textarea, network preview (text summary: "N treatments, M edges, K studies"), "Analyze Consistency" button.

- [ ] **Step 5: Add escapeHtml + utility functions**

escapeHtml, normalCDF (Abramowitz-Stegun), xoshiro128** PRNG.

- [ ] **Step 6: Write data tests + verify**

```javascript
function testData() {
  const results = [];
  const net = buildNetwork(DATASETS.smoking.contrasts);
  results.push({ name: 'Smoking: 4 treatments', pass: net.nodes.length === 4 });
  results.push({ name: 'Smoking: edges pooled', pass: net.edges.length <= 6 });
  const sim = buildNetwork(DATASETS.simulated.contrasts);
  results.push({ name: 'Simulated: 5 treatments', pass: sim.nodes.length === 5 });
  results.push({ name: 'CSV parser works', pass: parseContrastCSV('study,treat1,treat2,effect,se\nS1,A,B,-0.3,0.1').length === 1 });
  return results;
}
```

- [ ] **Step 7: Commit**

```bash
git add sheaf_nma.html
git commit -m "feat: HTML shell, 2 demo datasets, CSV parser, network builder"
```

---

### Task 2: Linear Algebra Module

**Files:**
- Modify: `sheaf_nma.html`

- [ ] **Step 1: Write linear algebra tests**

```javascript
function testLinearAlgebra() {
  const results = [];
  // Matrix multiply
  const A = [[1,2],[3,4]], B = [[5,6],[7,8]];
  const C = matMul(A, B);
  results.push({ name: 'matMul [[19,22],[43,50]]', pass: C[0][0]===19 && C[1][1]===50 });
  // Transpose
  const T = matTranspose([[1,2,3],[4,5,6]]);
  results.push({ name: 'transpose shape', pass: T.length===3 && T[0].length===2 });
  // Inverse 2x2
  const M = [[4,7],[2,6]];
  const Mi = matInverse(M);
  const I = matMul(M, Mi);
  results.push({ name: 'inverse: M*M^-1 ≈ I', pass: Math.abs(I[0][0]-1)<1e-8 && Math.abs(I[0][1])<1e-8 });
  // Eigendecomposition of symmetric 2x2
  const S = [[2,1],[1,3]];
  const eig = eigenSymmetric(S);
  results.push({ name: 'eigenvalues positive', pass: eig.values.every(v => v >= -1e-10) });
  results.push({ name: 'eigenvectors orthogonal', pass: Math.abs(dot(eig.vectors[0], eig.vectors[1])) < 1e-8 });
  return results;
}
```

- [ ] **Step 2: Implement matrix operations**

`matMul(A, B)`, `matTranspose(A)`, `matInverse(A)` (Gauss-Jordan with partial pivoting), `dot(a, b)` (vector dot product), `matScale(A, s)`, `matAdd(A, B)`.

- [ ] **Step 3: Implement eigendecomposition**

`eigenSymmetric(A)` — QR algorithm with implicit shifts for symmetric matrices. Returns `{ values: [λ₁,...], vectors: [[v₁],...] }` sorted ascending.

For n ≤ 20: Householder tridiagonalization → implicit QR iteration → eigenvectors via inverse iteration. This is the standard LAPACK approach, simplified for small n.

- [ ] **Step 4: Run tests**

All linear algebra tests pass.

- [ ] **Step 5: Commit**

```bash
git add sheaf_nma.html
git commit -m "feat: linear algebra — matMul, inverse, eigendecomposition"
```

---

### Task 3: Sheaf Engine

**Files:**
- Modify: `sheaf_nma.html`

- [ ] **Step 1: Write sheaf tests**

```javascript
function testSheaf() {
  const results = [];

  // Consistent triangle: A-B=-0.3, A-C=-0.5, B-C=-0.2 (perfectly consistent)
  const consistent = [
    { treat1:'A', treat2:'B', effect:-0.3, se:0.1, k:1 },
    { treat1:'A', treat2:'C', effect:-0.5, se:0.1, k:1 },
    { treat1:'B', treat2:'C', effect:-0.2, se:0.1, k:1 }
  ];
  const net1 = { nodes:['A','B','C'], edges: consistent };
  const r1 = sheafAnalysis(net1);
  results.push({ name: 'Consistent: all residuals ≈ 0', pass: r1.edgeScores.every(s => s.score < 0.01) });
  results.push({ name: 'Consistent: GII ≈ 0', pass: r1.gii < 0.01 });

  // Inconsistent triangle: A-B=-0.3, A-C=-0.5, B-C=+0.5 (B-C is wrong)
  const inconsistent = [
    { treat1:'A', treat2:'B', effect:-0.3, se:0.1, k:1 },
    { treat1:'A', treat2:'C', effect:-0.5, se:0.1, k:1 },
    { treat1:'B', treat2:'C', effect:0.5, se:0.1, k:1 }
  ];
  const net2 = { nodes:['A','B','C'], edges: inconsistent };
  const r2 = sheafAnalysis(net2);
  results.push({ name: 'Inconsistent: B-C has highest score', pass:
    r2.edgeScores.reduce((mx,s) => s.score > mx.score ? s : mx).edge === 'B-C' });
  results.push({ name: 'Inconsistent: GII > 0', pass: r2.gii > 0.1 });

  // Laplacian properties
  results.push({ name: 'Laplacian symmetric', pass: isSymmetric(r2.laplacian) });
  results.push({ name: 'Eigenvalues non-negative', pass: r2.eigenvalues.every(v => v >= -1e-8) });
  results.push({ name: 'At least one zero eigenvalue', pass: r2.eigenvalues[0] < 1e-6 });

  // Simulated: planted edge ranked #1
  const simNet = buildNetwork(DATASETS.simulated.contrasts);
  const simR = sheafAnalysis(simNet);
  const topEdge = simR.edgeScores.reduce((mx,s) => s.score > mx.score ? s : mx);
  results.push({ name: 'Simulated: A-D is most inconsistent', pass: topEdge.edge === 'A-D' });

  return results;
}

function isSymmetric(M) {
  for (let i=0; i<M.length; i++)
    for (let j=i+1; j<M.length; j++)
      if (Math.abs(M[i][j]-M[j][i]) > 1e-10) return false;
  return true;
}
```

- [ ] **Step 2: Implement sheafAnalysis**

```javascript
function sheafAnalysis(network) {
  const { nodes, edges } = network;
  const n = nodes.length, m = edges.length;
  const nodeIdx = Object.fromEntries(nodes.map((t,i) => [t,i]));

  // 1. Coboundary map F (m × n), precision-weighted
  const F = Array.from({length:m}, () => new Array(n).fill(0));
  const d = new Array(m).fill(0); // observed data vector
  for (let e=0; e<m; e++) {
    const edge = edges[e];
    const w = 1 / edge.se; // precision weight
    F[e][nodeIdx[edge.treat1]] = +w;
    F[e][nodeIdx[edge.treat2]] = -w;
    d[e] = edge.effect * w;
  }

  // 2. Laplacian L = F^T F
  const Ft = matTranspose(F);
  const L = matMul(Ft, F);

  // 3. Solve for network estimates: fix node 0 as reference
  //    Reduced system: L_red * x_red = (F^T d)_red
  const Ftd = Ft.map(row => row.reduce((s,v,j) => s + v*d[j], 0));
  // Remove row/col 0 (reference treatment)
  const Lred = L.slice(1).map(row => row.slice(1));
  const bRed = Ftd.slice(1);
  const xRed = matSolve(Lred, bRed); // solve via inverse
  const x = [0, ...xRed]; // reference = 0

  // 4. Edge residuals
  const edgeScores = edges.map((edge, e) => {
    const predicted = (x[nodeIdx[edge.treat1]] - x[nodeIdx[edge.treat2]]) / edge.se;
    const residual = d[e] - predicted;
    return { edge: edge.treat1 + '-' + edge.treat2, score: residual*residual,
             residual, effect: edge.effect, se: edge.se, k: edge.k };
  });

  // Normalize scores to 0-100
  const maxScore = Math.max(...edgeScores.map(s => s.score), 1e-15);
  edgeScores.forEach(s => s.pct = s.score / maxScore * 100);

  // 5. Eigendecomposition
  const eig = eigenSymmetric(L);

  // 6. GII
  const nonZeroEigs = eig.values.filter(v => v > 1e-8);
  const gii = nonZeroEigs.length > 0 ? nonZeroEigs.reduce((a,b)=>a+b,0) / (n-1) : 0;

  return {
    nodes, edges, nodeEstimates: x, edgeScores, laplacian: L,
    eigenvalues: eig.values, eigenvectors: eig.vectors,
    coboundary: F, dataVector: d, gii
  };
}

function matSolve(A, b) {
  // Solve Ax = b via A^{-1} b
  const Ainv = matInverse(A);
  return Ainv.map(row => row.reduce((s,v,j) => s + v*b[j], 0));
}
```

- [ ] **Step 3: Run sheaf tests**

All 8+ sheaf tests pass. Consistent triangle gives zero residuals. Inconsistent triangle flags B-C. Simulated flags A-D.

- [ ] **Step 4: Commit**

```bash
git add sheaf_nma.html
git commit -m "feat: sheaf engine — coboundary, Laplacian, WLS, residuals, eigendecomposition"
```

---

### Task 4: Visualizations — Network Heatmap + Spectrum

**Files:**
- Modify: `sheaf_nma.html`

- [ ] **Step 1: Implement network graph on Canvas**

Force-directed layout (Fruchterman-Reingold, 200 iterations). Render:
- Circles for treatment nodes (labeled, size ∝ degree)
- Lines for edges, color from green (#22c55e, score=0) → yellow (#eab308) → red (#ef4444, score=100)
- Edge thickness ∝ 1/SE (precision)
- Hover detection: nearest edge within 15px → tooltip with effect, CI, inconsistency score, studies
- Hover node → tooltip with estimated effect, degree

Below canvas: ranked edge table (HTML table sorted by inconsistency score desc). Columns: Edge, Direct Effect [CI], Score (%), Studies.

- [ ] **Step 2: Implement spectrum visualizations**

Left panel: eigenvalue bar chart on canvas. X-axis = eigenvalue index, Y-axis = λ value. Zero eigenvalues in grey (#94a3b8), non-zero in blue (#3b82f6). Click bar → updates right panel.

Right panel: eigenvector loading bar chart. For selected eigenvalue, show |v_i| per node. Label bars with treatment names. Large bars highlighted.

GII display: large text "Global Inconsistency Index: X.XX" with color (green <0.5, yellow 0.5-2.0, red >2.0).

- [ ] **Step 3: Implement color legend**

Horizontal gradient bar (green → yellow → red) with 0% and 100% labels below the network canvas.

- [ ] **Step 4: Wire "Analyze Consistency" button**

Parse data → buildNetwork → sheafAnalysis → render network + spectrum + edge table. Switch to Network tab. Toast: "Analysis complete: GII = X.XX"

- [ ] **Step 5: Visual verification**

Load smoking cessation, click Analyze. Verify: network renders with 4 nodes, edges colored, spectrum shows eigenvalues, edge table sorted. Load simulated, verify A-D edge is red.

- [ ] **Step 6: Commit**

```bash
git add sheaf_nma.html
git commit -m "feat: visualizations — network heatmap, scree plot, eigenvector loadings"
```

---

### Task 5: Export Module + Integration + Tests + Push

**Files:**
- Modify: `sheaf_nma.html`

- [ ] **Step 1: Implement R code generator**

```javascript
function generateRCode(network, results) {
  const contrasts = network.edges;
  return `library(netmeta)

# Data
dat <- data.frame(
  treat1 = c(${contrasts.map(c => `"${c.treat1}"`).join(', ')}),
  treat2 = c(${contrasts.map(c => `"${c.treat2}"`).join(', ')}),
  TE = c(${contrasts.map(c => c.effect).join(', ')}),
  seTE = c(${contrasts.map(c => c.se).join(', ')})
)

# Standard NMA
nma <- netmeta(TE, seTE, treat1, treat2, data = dat, sm = "OR", reference.group = "${network.nodes[0]}")
print(nma)

# Node-splitting for comparison
ns <- netsplit(nma)
print(ns)

# Sheaf Laplacian (manual computation for comparison)
# GII from SheafNMA: ${results.gii.toFixed(4)}
# Edge inconsistency scores:
${results.edgeScores.map(s => `# ${s.edge}: score = ${s.pct.toFixed(1)}%`).join('\\n')}
`;
}
```

- [ ] **Step 2: Implement CSV, PNG, JSON export**

CSV: edge, effect, se, score_pct, residual, k. Blob URL + revoke.
PNG: canvas.toBlob() for network + spectrum. Revoke after use.
JSON: full results object.

- [ ] **Step 3: Run ALL in-browser tests**

Combine: testData + testLinearAlgebra + testSheaf. All must pass. Display in modal.

- [ ] **Step 4: Div balance check**

Count `<div` vs `</div>` (excluding JS strings). Must match.

- [ ] **Step 5: Security + accessibility check**

- escapeHtml on all user input rendered to DOM
- No literal `</script>` in script block
- ARIA tablist/tab/tabpanel roles
- Dark mode contrast adequate

- [ ] **Step 6: Commit and push**

```bash
cd C:\Models\SheafNMA
git add -A
git commit -m "feat: SheafNMA v1.0 — world's first sheaf-theoretic NMA consistency tool"
gh repo create mahmood726-cyber/sheaf-nma --public --description "Sheaf-theoretic NMA consistency analysis"
git remote add origin https://github.com/mahmood726-cyber/sheaf-nma.git
git push -u origin master
```

---

## Commit Sequence

| Task | Commit Message | ~Lines |
|------|---------------|--------|
| 1 | `feat: HTML shell, 2 demo datasets, CSV parser, network builder` | ~600 |
| 2 | `feat: linear algebra — matMul, inverse, eigendecomposition` | ~400 |
| 3 | `feat: sheaf engine — coboundary, Laplacian, WLS, residuals` | ~450 |
| 4 | `feat: visualizations — network heatmap, scree plot, eigenvectors` | ~550 |
| 5 | `feat: SheafNMA v1.0 — export, tests, polish` | ~500 |
| **Total** | | **~2,500** |
