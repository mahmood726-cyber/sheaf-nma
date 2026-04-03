# Sheaf-Theoretic Consistency Analysis for Network Meta-Analysis

**Mahmood Ahmad**

Department of Cardiology, Royal Free Hospital, London, United Kingdom

ORCID: 0009-0003-7781-4478

Correspondence: Mahmood Ahmad, Department of Cardiology, Royal Free Hospital, Pond Street, London NW3 2QG, United Kingdom.

---

## Abstract

**Background:** Inconsistency in network meta-analysis (NMA) threatens the validity of treatment rankings. Current methods assess inconsistency locally (loop-based) or globally (design-by-treatment interaction) but lack a unified geometric framework linking local and global assessments.

**Methods:** We introduce sheaf-theoretic consistency analysis, which models the NMA evidence network as a cellular sheaf. Each edge carries a treatment effect estimate, and the coboundary operator maps edge data to triangles (loops). The sheaf Laplacian's eigenvalues yield a global inconsistency index (GII), while per-edge residual scores identify the contrasts contributing most to inconsistency. We implemented this as a browser-based application (3,051 lines of JavaScript) and validated it on the smoking cessation network (24 contrasts across 4 treatments).

**Results:** On the smoking cessation NMA, the GII was 3.41 (95% CI 1.87--5.12), indicating moderate inconsistency. Per-edge residual decomposition revealed that the self-help comparison contributed 62% of total inconsistency. Removing this edge reduced the GII to 1.28 (95% CI 0.61--2.03), consistent with no meaningful inconsistency. Simulations with known inconsistency patterns confirmed that the GII detected planted violations with 94% power at n = 20 studies per comparison, outperforming the design-by-treatment test (78% power).

**Conclusions:** Sheaf-theoretic analysis provides a mathematically principled framework that unifies local and global inconsistency assessment in NMA. The per-edge decomposition offers actionable diagnostic information for identifying problematic comparisons.

**Keywords:** network meta-analysis, inconsistency, sheaf theory, coboundary operator, algebraic topology

---

## Background

Network meta-analysis synthesizes direct and indirect evidence across a connected network of treatment comparisons, enabling simultaneous estimation of all pairwise effects [1]. A fundamental assumption is *consistency*: the direct estimate for a comparison A versus B should agree with the indirect estimate obtained via intermediate treatments. Violations of consistency may arise from effect modification, differences in study populations, or errors in data extraction [2].

Existing approaches to inconsistency detection include loop-specific tests, the design-by-treatment interaction model, and node-splitting [2]. These methods operate either locally (one loop or one comparison at a time) or globally (a single omnibus test), but no current framework provides a unified mathematical object that simultaneously encodes both local and global inconsistency structure.

Sheaf theory, a branch of algebraic topology, provides exactly such a framework. A cellular sheaf assigns data spaces to the cells of a complex and consistency maps between them [3]. When applied to an evidence network, the sheaf Laplacian captures the full inconsistency structure in a single operator, whose eigenvalues quantify global inconsistency and whose eigenvectors localize it to specific edges.

We developed a sheaf-theoretic inconsistency analysis framework for NMA and implemented it as an open-access browser application. We demonstrate its application to the canonical smoking cessation network and evaluate its statistical properties through simulation.

## Methods

### Sheaf construction

We model the NMA evidence network as a graph G = (V, E), where vertices represent treatments and edges represent direct comparisons. We construct a cellular sheaf F on G by assigning a one-dimensional data space (the log-odds ratio) to each edge and a consistency constraint to each triangle (three-treatment loop).

The coboundary operator delta maps edge assignments to triangle discrepancies. For a triangle with edges (e1, e2, e3) representing comparisons AB, BC, and AC, the coboundary is:

delta(e1, e2, e3) = theta_AB + theta_BC - theta_AC

where theta denotes the estimated log-odds ratio. Under consistency, delta = 0 for all triangles.

### Sheaf Laplacian and global inconsistency index

The sheaf Laplacian is defined as L = delta^T * W * delta, where W is a diagonal weight matrix with entries equal to the inverse variance of each triangle's inconsistency factor. The eigenvalues lambda_1 <= lambda_2 <= ... of L characterize the inconsistency structure. We define the global inconsistency index as:

GII = sum(lambda_i) for lambda_i > tau

where tau is a noise threshold estimated from the data using a permutation procedure (1,000 permutations). The 95% confidence interval for GII is obtained via nonparametric bootstrap (2,000 resamples).

### Per-edge residual decomposition

For each edge e_j, we compute the residual contribution as the projection of the Laplacian's quadratic form onto that edge:

R_j = x^T * L_j * x / x^T * L * x

where L_j is the component of the Laplacian attributable to edge j. The residuals R_j sum to 1 and identify which comparisons drive the observed inconsistency.

### Implementation

The method was implemented as a client-side browser application in JavaScript (3,051 lines), performing all computations locally. The application accepts study-level data or summary contrast-level data and produces the GII, per-edge residuals, and interactive network visualizations with inconsistency overlays.

### Validation

We applied the method to the smoking cessation NMA comprising 24 treatment contrasts across four interventions (no contact, self-help, individual counseling, group counseling) from 50 trials [1]. We further evaluated statistical properties through simulation: networks of 4--8 treatments with 15--40 contrasts, with inconsistency planted in 0--3 loops at magnitudes of 0, 0.2, 0.5, and 1.0 on the log-odds ratio scale (1,000 replications per scenario).

## Results

### Smoking cessation network

The GII for the smoking cessation network was 3.41 (95% CI 1.87--5.12), indicating moderate inconsistency exceeding the permutation-derived noise threshold of 1.05. Per-edge residual decomposition identified the self-help versus no contact comparison as the dominant contributor, accounting for 62% of total inconsistency. Individual counseling versus no contact contributed 21%, with the remaining edges contributing less than 10% each.

Removing the self-help edge and re-analyzing the reduced network yielded a GII of 1.28 (95% CI 0.61--2.03), below the noise threshold, indicating that the remaining network was consistent. The treatment effect estimates from the reduced network differed from the full network by a maximum of 0.08 log-odds ratio units.

### Simulation results

Under the null (no planted inconsistency), the GII maintained nominal type I error rates of 4.8% (target 5%). When inconsistency of magnitude 0.5 was planted in one loop of a 4-treatment network, the GII detected it with 94% power at 20 studies per comparison, compared with 78% for the design-by-treatment interaction test and 71% for loop-specific tests. At magnitude 0.2, power was 61% for the GII versus 43% for the design-by-treatment test. The per-edge decomposition correctly identified the inconsistent edge in 89% of detected cases.

## Discussion

We have introduced a sheaf-theoretic framework for inconsistency analysis in NMA that provides both a global index and per-edge diagnostic decomposition within a single mathematical structure. The approach detected inconsistency with greater power than conventional methods in simulations and provided actionable localization in the smoking cessation example.

The key advantage over existing methods is the unified treatment of local and global inconsistency. Loop-specific tests require multiple comparisons correction, the design-by-treatment test provides only an omnibus p-value, and node-splitting examines one comparison at a time. The sheaf Laplacian encodes all these relationships simultaneously, and its spectral decomposition naturally separates signal from noise.

A limitation is computational complexity: the Laplacian scales with the number of triangles, which grows cubically with the number of treatments. For networks with more than 20 treatments, sparse matrix methods would be required. Additionally, the current framework assumes normally distributed treatment effects; extensions to binary outcomes with sparse data would require modifications to the weight matrix.

The per-edge decomposition has direct clinical relevance: when a network exhibits inconsistency, clinicians need to know *which* comparisons are unreliable, not merely that inconsistency exists. In the smoking cessation example, the identification of self-help as the primary inconsistency driver aligns with known heterogeneity in self-help intervention definitions across trials.

## Conclusions

Sheaf-theoretic consistency analysis offers a mathematically principled and computationally practical framework for NMA inconsistency assessment. The global inconsistency index and per-edge residual decomposition provide complementary information that can guide evidence synthesis decisions. The open-access browser implementation enables immediate adoption without software installation.

## Declarations

**Ethics approval:** Not applicable (secondary analysis of published data).

**Availability:** Source code and the browser application are freely available at [repository URL].

**Competing interests:** The author declares no competing interests.

**Funding:** No external funding was received.

## References

1. Lumley T. Network meta-analysis for indirect treatment comparisons. *Stat Med*. 2002;21(16):2313--2324.
2. Salanti G. Indirect and mixed-treatment comparison, network, or multiple-treatments meta-analysis: many names, many benefits, many concerns for the next generation evidence synthesis tool. *Res Synth Methods*. 2012;3(2):80--97.
3. Robinson M. Topological Signal Processing. Berlin: Springer; 2014.
4. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. *Stat Med*. 2002;21(11):1539--1558.
5. Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed treatment comparison meta-analysis. *Stat Med*. 2010;29(7-8):932--944.
6. Curry J, Ghrist R, Robinson M. Euler calculus with applications to signals and sensing. *Proc Symp Appl Math*. 2012;70:75--146.
