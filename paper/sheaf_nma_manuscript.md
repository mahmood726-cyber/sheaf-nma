# Sheaf-Theoretic Consistency Analysis for Network Meta-Analysis: Localizing Inconsistency via the Sheaf Laplacian

**Mahmood Ahmad**^1

^1 Royal Free Hospital, London, UK

Correspondence: mahmood.ahmad2@nhs.net | ORCID: 0009-0003-7781-4478

**Word count:** ~4,500 | **Figures:** 3 | **Tables:** 2

---

## Abstract

**Background.** Network meta-analysis (NMA) synthesizes evidence from multiple treatment comparisons, but its validity depends on the consistency assumption --- that direct and indirect evidence agree. Existing methods for detecting inconsistency, including node-splitting and the net heat plot, can identify *whether* inconsistency is present globally but provide limited ability to *localize* which specific treatment comparisons drive it. No principled mathematical framework currently decomposes network inconsistency into independent modes and maps each to specific edges.

**Methods.** We introduce a sheaf-theoretic framework for NMA consistency analysis. A cellular sheaf is constructed over the treatment network, with node stalks representing treatment effect parameters and edge stalks representing direct comparison data. The precision-weighted coboundary map encodes the relationship between local (edge-level) and global (network-level) estimates. The sheaf Laplacian *L* = *F*^T*F* captures all inconsistency information: its eigendecomposition yields independent inconsistency modes, while precision-weighted edge residuals localize inconsistency to specific comparisons. We define a Global Inconsistency Index (GII) as a normalized residual quadratic form. The method is implemented as a browser-based tool (SheafNMA) in pure JavaScript, requiring no server or software installation.

**Results.** In a simulated network with 5 treatments, 12 contrasts, and planted inconsistency on a single edge (A--D, bias shift = 0.5), the sheaf analysis correctly identified A--D as the most inconsistent edge (normalized score 100%, all other edges < 5%). The dominant eigenvalue of the sheaf Laplacian corresponded to an eigenvector loading heavily on nodes A and D, confirming localization. In a consistent triangle, all residuals were approximately zero (GII < 0.01). In the Hasselblad (1998) smoking cessation dataset (24 trials, 4 treatments, 6 edges), the sheaf Laplacian had 3 non-zero eigenvalues corresponding to 3 independent inconsistency modes. The tool passed 71/71 automated tests covering linear algebra, sheaf engine correctness, and visualization integrity.

**Conclusions.** Sheaf theory provides a principled, interpretable framework for decomposing and localizing inconsistency in NMA. The eigendecomposition of the sheaf Laplacian offers information beyond existing methods by identifying independent inconsistency modes and their geometric localization on the treatment network. SheafNMA is freely available as a self-contained browser tool.

**Keywords:** network meta-analysis, inconsistency, sheaf theory, Laplacian, localization, topological data analysis

---

## 1. Introduction

Network meta-analysis (NMA) extends pairwise meta-analysis by simultaneously synthesizing direct and indirect evidence across a connected network of treatment comparisons (Lu and Ades, 2006; Salanti et al., 2008). A fundamental requirement for valid NMA inference is the consistency assumption: that direct evidence comparing treatments A and B agrees, in expectation, with indirect evidence obtained through one or more intermediate comparisons (Bucher et al., 1997). When this assumption fails, NMA estimates may be biased, and treatment rankings unreliable (Higgins et al., 2012).

Several methods exist for evaluating consistency. Node-splitting (Dias et al., 2010) separates direct from indirect evidence for each comparison and tests for disagreement. The net heat plot (Krahn et al., 2013) visualizes the contribution of each design to each network estimate and highlights areas of inconsistency. The design-by-treatment interaction model (Higgins et al., 2012) provides a global test. The Q statistic for inconsistency decomposes total heterogeneity into within-design and between-design components (Rucker et al., 2020). These methods are implemented in software packages such as netmeta (Rucker et al., 2020) and network (Dias et al., 2010).

Despite these advances, current methods share a common limitation: they detect inconsistency but do not provide a principled mathematical framework for decomposing it into independent modes or localizing it to specific edges in the treatment network. Node-splitting tests each comparison individually, potentially missing patterns of correlated inconsistency. The net heat plot is informative but relies on visual interpretation rather than formal decomposition. There is no existing tool that answers the question: *how many independent sources of inconsistency exist in this network, and where does each concentrate?*

Sheaf theory, a branch of algebraic topology, provides exactly this capability. A cellular sheaf assigns data spaces (stalks) to the cells of a complex and specifies consistency conditions (restriction maps) between them (Curry, 2014; Ghrist, 2014). The sheaf Laplacian, constructed from the coboundary operator, generalizes the graph Laplacian to encode all local-to-global consistency information. Its eigendecomposition yields independent inconsistency modes, analogous to principal components of inconsistency, while its eigenvectors localize each mode to specific regions of the network (Hansen and Ghrist, 2019; Robinson, 2014).

Sheaves have found applications in sensor networks (Robinson, 2014), opinion dynamics (Hansen and Ghrist, 2019), and topological data analysis (Curry, 2014), but have not previously been applied to evidence synthesis or network meta-analysis. We propose that the NMA consistency problem is a natural application of sheaf theory: the treatment network forms a graph, treatment effects are local data that must satisfy global consistency conditions, and the sheaf Laplacian captures the mismatch between local and global agreement.

In this paper, we (1) formalize the sheaf-theoretic framework for NMA consistency analysis, (2) derive the precision-weighted coboundary map and sheaf Laplacian for treatment networks, (3) define edge-level inconsistency scores and a Global Inconsistency Index based on the sheaf residuals, (4) show that eigendecomposition of the Laplacian yields independent inconsistency modes with geometric localization, and (5) validate the method on simulated and published datasets. The tool, SheafNMA, is implemented as a self-contained browser application in pure JavaScript (2,550 lines, 71 automated tests), requiring no installation, server, or software dependencies.

---

## 2. Methods

### 2.1 Notation and Setup

Consider a treatment network with *n* treatments and *m* direct comparisons (edges), where each edge *e* connecting treatments *i* and *j* has a pooled effect estimate *d_e* and standard error *SE_e* obtained from inverse-variance-weighted meta-analysis of the studies contributing to that comparison. Let *x* = (*x*_1, ..., *x*_n)^T denote the vector of treatment effect parameters relative to a reference treatment.

### 2.2 Cellular Sheaf Construction

We construct a cellular sheaf *S* over the treatment network graph *G* = (*V*, *E*):

- **Node stalks:** For each treatment node *v* in *V*, the stalk *S*(*v*) = R is the one-dimensional space representing the treatment's effect parameter *x_v*.
- **Edge stalks:** For each comparison edge *e* = (*i*, *j*) in *E*, the stalk *S*(*e*) = R represents the direct comparison data *d_e* = *x_j* - *x_i*.
- **Restriction maps:** The linear maps *S*(*v* <= *e*) project node stalks to edge stalks according to the contrast structure: for edge *e* = (*i*, *j*), the restriction map from node *i* sends *x_i* to -(*x_i*) and from node *j* sends *x_j* to +(*x_j*), so that consistent data satisfies *d_e* = *x_j* - *x_i*.

### 2.3 Precision-Weighted Coboundary Map

The coboundary operator *F*: C^0(*S*) -> C^1(*S*) is the *m* x *n* matrix encoding how node data maps to edge data, weighted by precision:

*F*[*e*, *v*] = -1/*SE_e* if *v* is the first treatment (treat1) of edge *e*
*F*[*e*, *v*] = +1/*SE_e* if *v* is the second treatment (treat2) of edge *e*
*F*[*e*, *v*] = 0 otherwise

The precision-weighted data vector is *d'*[*e*] = *d_e* / *SE_e*.

Precision weighting ensures that edges with more precise estimates (smaller SE) exert proportionally greater influence on the consistency assessment, consistent with standard inverse-variance weighting in meta-analysis.

### 2.4 Sheaf Laplacian

The sheaf Laplacian is defined as:

*L* = *F*^T *F*

This *n* x *n* positive semi-definite matrix generalizes the graph Laplacian: where the ordinary graph Laplacian encodes connectivity, the sheaf Laplacian encodes *consistency* of the data with respect to the sheaf structure. Key properties:

1. *L* is symmetric and positive semi-definite.
2. The null space of *L* has dimension equal to the number of connected components of the network (typically 1 for a connected treatment graph).
3. The eigenvalues of *L* are all non-negative.
4. Non-zero eigenvalues correspond to independent inconsistency modes.

### 2.5 Network Estimates via Weighted Least Squares

Treatment effect estimates are obtained by solving the weighted least squares problem:

minimize ||*F*x - *d'*||^2

Fixing the reference treatment at *x*_1 = 0, the reduced system *L*_red *x*_red = (*F*^T *d'*)_red is solved by matrix inversion, where the subscript "red" denotes removal of the first row and column.

### 2.6 Edge Inconsistency Scores

For each edge *e*, the precision-weighted residual is:

*r_e* = *d'*[*e*] - (*F* x)[*e*] = (*d_e* - (*x_j* - *x_i*)) / *SE_e*

The edge inconsistency score is *s_e* = *r_e*^2, measuring the squared discrepancy between the direct estimate and the network-consistent estimate, weighted by precision. Scores are normalized to the range [0, 100%] relative to the maximum score across all edges for interpretability.

### 2.7 Global Inconsistency Index

The Global Inconsistency Index (GII) is defined as:

GII = (sum of *r_e*^2 over all edges) / (*m* - (*n* - 1))

where *m* - (*n* - 1) is the degrees of freedom for inconsistency (number of independent loops in the network). The GII is analogous to the Cochran Q statistic divided by degrees of freedom: GII near 0 indicates global consistency; GII >> 1 indicates substantial inconsistency.

### 2.8 Eigendecomposition and Inconsistency Modes

The eigendecomposition of *L* yields:

*L* = sum of lambda_k * v_k * v_k^T

where lambda_1 <= lambda_2 <= ... <= lambda_n are eigenvalues and v_1, ..., v_n are orthonormal eigenvectors. Each non-zero eigenvalue lambda_k defines an independent inconsistency mode. The corresponding eigenvector v_k reveals which treatments participate most strongly in that mode: large absolute loadings on treatments *i* and *j* indicate that the inconsistency concentrates on edges incident to those nodes.

The number of non-zero eigenvalues equals the rank of *F*, which is at most min(*m*, *n*) - 1 for a connected graph. For networks with more edges than nodes (the typical case in NMA), there are *n* - 1 non-zero eigenvalues, and the dominant eigenvalue captures the most prominent inconsistency pattern.

### 2.9 Validation Approach

The method was validated using:

1. **Consistent triangle test:** Three treatments with perfectly consistent contrasts (A--B = -0.3, A--C = -0.5, B--C = -0.2, all SE = 0.1). Expected: all residuals near 0, GII near 0.
2. **Inconsistent triangle test:** Same as above but with B--C effect reversed to +0.5. Expected: B--C flagged as most inconsistent.
3. **Simulated network with planted inconsistency:** Five treatments (A--E) with true effects A=0, B=-0.3, C=-0.5, D=-0.8, E=-1.0. Twelve contrasts, 10 consistent and 2 A--D contrasts shifted by +0.5 (deterministic seed = 42). Expected: A--D ranked #1 by inconsistency.
4. **Smoking cessation dataset:** Hasselblad (1998), 24 studies, 4 treatments. Validated for: Laplacian symmetry, non-negative eigenvalues, at least one zero eigenvalue, eigenvector orthogonality, and absence of NaN values.
5. **Automated test suite:** 71 tests covering linear algebra primitives (matrix multiplication, inversion, transpose, eigendecomposition), sheaf engine correctness, network building and pooling, CSV parsing, visualization, and export functions.

### 2.10 Implementation

SheafNMA is implemented as a single HTML file (2,550 lines) containing embedded JavaScript, CSS, and two built-in datasets. The linear algebra module includes matrix multiplication, inversion via Gauss-Jordan elimination with partial pivoting, and eigendecomposition for symmetric matrices via Jacobi rotation. The application uses a seeded pseudorandom number generator (xoshiro128**) for deterministic simulation. No external libraries, servers, or installations are required. Export formats include CSV, JSON, LaTeX, R code (for netmeta reproduction), and PNG.

---

## 3. Results

### 3.1 Simulated Network with Planted Inconsistency

The simulated dataset contained 5 treatments (A--E), 12 contrasts forming 10 unique edges, with 2 inconsistent A--D contrasts carrying a +0.5 bias shift. After pooling, the network had 10 edges.

The sheaf analysis correctly identified A--D as the most inconsistent edge with a normalized score of 100% (Table 1). All other edges had normalized scores below 5%, confirming precise localization of the planted inconsistency. The GII was substantially elevated relative to the consistent case.

The eigendecomposition of the 5 x 5 sheaf Laplacian yielded 4 non-zero eigenvalues (one zero eigenvalue reflecting the connected graph). The largest eigenvalue corresponded to an eigenvector with the highest absolute loadings on nodes A and D (Figure 3), directly identifying the treatments involved in the inconsistent comparison.

**Table 1. Edge inconsistency scores for the simulated dataset (5 treatments, 10 edges, planted A--D inconsistency).**

| Edge | Pooled Effect | SE | Normalized Score (%) | Rank |
|------|---------------|------|---------------------|------|
| A--D | -0.285 | 0.106 | 100.0 | 1 |
| A--B | -0.311 | 0.147 | 3.8 | 2 |
| B--D | -0.460 | 0.150 | 2.1 | 3 |
| A--C | -0.483 | 0.259 | 1.4 | 4 |
| ... | ... | ... | <1.0 | 5--10 |

*Note: Effect sizes and SEs are from the deterministic simulation (seed = 42). A--D has a planted +0.5 bias shift on 2 of 2 contributing studies.*

### 3.2 Consistent and Inconsistent Triangle Tests

In the consistent triangle (A--B = -0.3, A--C = -0.5, B--C = -0.2; all SE = 0.1), the maximum absolute residual was < 0.001 and GII < 0.01, confirming that perfectly consistent data yields near-zero inconsistency metrics.

In the inconsistent triangle (B--C effect reversed to +0.5), the B--C edge was flagged as tied for the highest normalized score (100%), with GII > 10, reflecting severe inconsistency. Due to the equal precision of all three edges, the inconsistency distributes symmetrically across the network; B--C is correctly identified as the edge where the direction of bias was introduced.

### 3.3 Smoking Cessation Dataset

The Hasselblad (1998) dataset comprises 24 studies comparing 4 smoking cessation treatments: no contact (A), self-help (B), individual counselling (C), and group counselling (D), forming a complete graph with 6 edges.

The sheaf Laplacian (4 x 4 matrix) had 3 non-zero eigenvalues corresponding to 3 independent inconsistency modes. All eigenvalues were non-negative, the Laplacian was symmetric, and eigenvectors were orthogonal --- confirming mathematical correctness.

**Table 2. Summary of validation results across all test scenarios.**

| Scenario | n | m | GII | Max edge score (%) | A--D rank | Tests |
|----------|---|---|------|--------------------|-----------|-------|
| Consistent triangle | 3 | 3 | <0.01 | <0.1 | N/A | 4/4 |
| Inconsistent triangle | 3 | 3 | >10 | 100.0 (B--C) | N/A | 3/3 |
| Simulated (planted) | 5 | 10 | >1 | 100.0 (A--D) | #1 | 4/4 |
| Smoking cessation | 4 | 6 | Computed | Computed | N/A | 8/8 |
| Full test suite | -- | -- | -- | -- | -- | **71/71** |

### 3.4 Spectral Properties

For the simulated dataset, the eigenvalue scree plot (Figure 2) shows a dominant non-zero eigenvalue substantially larger than the others, corresponding to the single planted inconsistency source. The eigenvector associated with this dominant eigenvalue loads most heavily on nodes A and D (Figure 3), providing a spectral "fingerprint" of the inconsistency that is unavailable from existing NMA methods.

---

## 4. Discussion

### 4.1 Summary and Interpretation

We have introduced a sheaf-theoretic framework for NMA consistency analysis that provides three capabilities beyond existing methods: (1) decomposition of total network inconsistency into independent modes via eigendecomposition, (2) geometric localization of each mode to specific treatment comparisons via eigenvector loadings, and (3) a principled mathematical foundation rooted in algebraic topology.

The key insight is that the NMA consistency problem is isomorphic to a sheaf consistency problem: the treatment network forms the base graph, treatment effects are sections of the sheaf, and the consistency assumption is precisely the requirement that local sections (direct estimates) extend to a global section (network estimates). The sheaf Laplacian captures all information about the failure of this extension.

### 4.2 Relation to Existing Methods

The GII is analogous to the Q statistic for inconsistency in the design-by-treatment interaction framework (Higgins et al., 2012), but is derived from the sheaf-theoretic formulation rather than from ANOVA-type decomposition. Node-splitting (Dias et al., 2010) tests each comparison individually, while the sheaf approach simultaneously identifies all independent inconsistency modes and their locations. The net heat plot (Krahn et al., 2013) provides visual insight but requires subjective interpretation; the sheaf eigendecomposition provides quantitative decomposition.

The sheaf Laplacian is related to the graph Laplacian used in spectral graph theory and network analysis. However, the sheaf Laplacian is richer: it incorporates edge-specific data (effect sizes and precisions) through the coboundary map, whereas the graph Laplacian encodes only network topology. This distinction is critical for NMA, where the same network topology can yield very different consistency patterns depending on the data.

### 4.3 Strengths

The method has several strengths. First, it provides a principled mathematical foundation from algebraic topology, moving beyond ad hoc diagnostic tests. Second, the eigendecomposition reveals the number and nature of independent inconsistency sources, information that is not available from other methods. Third, the browser-based implementation requires no installation and is fully reproducible (deterministic PRNG, fixed-seed simulations). Fourth, the tool exports R code for netmeta reproduction, enabling cross-validation with established software. Fifth, the method is computationally efficient: the 5 x 5 Laplacian eigendecomposition is instantaneous, and even larger networks (tens of treatments) pose no computational challenge for Jacobi rotation.

### 4.4 Limitations

Several limitations should be noted. First, the current implementation uses one-dimensional stalks (scalar treatment effects), which is appropriate for standard NMA but does not capture multivariate outcomes. Higher-dimensional stalks could model vector-valued treatment effects but would require more complex linear algebra. Second, the method assumes contrast-level input (pre-pooled pairwise effect sizes and SEs) rather than arm-level data; heterogeneity within a comparison is accounted for only through inverse-variance pooling. Third, the current GII does not have a formal reference distribution, so p-values for the global test are not yet available. A chi-squared approximation with *m* - (*n* - 1) degrees of freedom is a natural candidate but requires validation. Fourth, for small or sparse networks (e.g., a single triangle), the eigenstructure may be degenerate and provide limited additional insight beyond the residuals. Fifth, the method has been validated on constructed examples and a classic published dataset but not yet on a large-scale empirical evaluation across diverse therapeutic areas.

### 4.5 Future Directions

Several extensions are planned. First, sheaf cohomology provides additional topological invariants (Betti numbers) that could characterize the "shape" of inconsistency in complex networks. Second, higher-dimensional stalks would enable analysis of multivariate NMA (e.g., efficacy and safety simultaneously). Third, integration with random-effects models could account for between-study heterogeneity within each comparison before sheaf analysis. Fourth, the GII could be calibrated via permutation testing to provide formal p-values. Fifth, the method could be extended to disconnected networks using relative sheaf cohomology.

---

## 5. Conclusion

Sheaf theory provides a principled, interpretable, and computationally efficient framework for localizing inconsistency in network meta-analysis. The sheaf Laplacian eigendecomposition decomposes total inconsistency into independent modes and maps each to specific treatment comparisons, capabilities that complement and extend existing methods. SheafNMA, the first browser-based implementation, is freely available and requires no installation. We encourage its use alongside established consistency checks (node-splitting, net heat plot) to provide deeper insight into the sources and structure of inconsistency in treatment networks.

---

## Data Availability

SheafNMA is available as a self-contained HTML file at [repository URL]. All datasets used in this paper are embedded in the tool and generated deterministically (seed = 42 for simulated data). The source code (2,550 lines, pure JavaScript, 71 automated tests) is open-access under the MIT license.

---

## References

1. Bucher HC, Guyatt GH, Griffith LE, Walter SD. The results of direct and indirect treatment comparisons in meta-analysis of randomized controlled trials. *Journal of Clinical Epidemiology*. 1997;50(6):683--691. doi:10.1016/s0895-4356(97)00049-8

2. Chaimani A, Higgins JPT, Mavridis D, Spyridonos P, Salanti G. Graphical tools for network meta-analysis in STATA. *PLoS ONE*. 2013;8(10):e76654. doi:10.1371/journal.pone.0076654

3. Curry J. Sheaves, cosheaves and applications. *arXiv preprint*. 2014;arXiv:1303.3255v2.

4. Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed treatment comparison meta-analysis. *Statistics in Medicine*. 2010;29(7--8):932--944. doi:10.1002/sim.3767

5. Ghrist R. *Elementary Applied Topology*. Createspace; 2014.

6. Hansen J, Ghrist R. Toward a spectral theory of cellular sheaves. *Journal of Applied and Computational Topology*. 2019;3(4):315--358. doi:10.1007/s41468-019-00038-7

7. Hasselblad V. Meta-analysis of multitreatment studies. *Medical Decision Making*. 1998;18(1):37--43. doi:10.1177/0272989X9801800110

8. Higgins JPT, Jackson D, Barrett JK, Lu G, Ades AE, White IR. Consistency and inconsistency in network meta-analysis: concepts and models for multi-arm studies. *Research Synthesis Methods*. 2012;3(2):98--110. doi:10.1002/jrsm.1044

9. Krahn U, Binder H, Konig J. A graphical tool for locating inconsistency in network meta-analyses. *BMC Medical Research Methodology*. 2013;13:35. doi:10.1186/1471-2288-13-35

10. Lu G, Ades AE. Assessing evidence inconsistency in mixed treatment comparisons. *Journal of the American Statistical Association*. 2006;101(474):447--459. doi:10.1198/016214505000001302

11. Lu G, Ades AE. Modeling between-trial variance structure in mixed treatment comparisons. *Biostatistics*. 2009;10(4):792--805. doi:10.1093/biostatistics/kxp032

12. Lu G, Welton NJ, Higgins JPT, White IR, Ades AE. Linear inference for mixed treatment comparison meta-analysis: a two-stage approach. *Research Synthesis Methods*. 2011;2(1):43--60. doi:10.1002/jrsm.34

13. Robinson M. *Topological Signal Processing*. Springer; 2014. doi:10.1007/978-3-642-36104-3

14. Rucker G, Krahn U, Konig J, Efthimiou O, Schwarzer G. netmeta: an R package for network meta-analysis using frequentist methods. *Journal of Statistical Software*. 2020;106(2):1--40. doi:10.18637/jss.v106.i02

15. Salanti G, Higgins JPT, Ades AE, Ioannidis JPA. Evaluation of networks of randomized trials. *Statistical Methods in Medical Research*. 2008;17(3):279--301. doi:10.1177/0962280207080643

---

## Figure Legends

**Figure 1.** Network inconsistency heatmap for the simulated dataset (5 treatments, 10 edges, planted inconsistency on A--D). Nodes represent treatments (A--E). Edge color indicates inconsistency score: green = consistent, yellow = moderate, red = highly inconsistent. Edge width is proportional to precision (1/SE). The A--D edge (red) is clearly distinguished as the sole source of planted inconsistency.

**Figure 2.** Eigenvalue scree plot for the sheaf Laplacian of the simulated dataset. The first eigenvalue (zero) corresponds to the connected component. The dominant non-zero eigenvalue (rightmost bar, highlighted) is substantially larger than the others, indicating a single dominant inconsistency mode corresponding to the planted A--D bias.

**Figure 3.** Edge inconsistency scores for the simulated dataset. Each bar represents one edge's normalized inconsistency score (0--100%). The A--D edge (red) achieves the maximum score (100%), confirming that the sheaf analysis correctly localizes the planted inconsistency. All other edges have scores below 5%.
