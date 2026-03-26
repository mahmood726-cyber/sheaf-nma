# Cover Letter

**To:** The Editor, *Journal of the Royal Statistical Society: Series A (Statistics in Society)*

**Date:** 26 March 2026

**Re:** Submission of manuscript "Sheaf-Theoretic Consistency Analysis for Network Meta-Analysis: Localizing Inconsistency via the Sheaf Laplacian"

---

Dear Editor,

I am pleased to submit the enclosed manuscript for consideration for publication in JRSS Series A.

**The problem.** Network meta-analysis (NMA) is widely used to compare multiple treatments simultaneously, but its validity depends on the consistency assumption. Current methods for detecting inconsistency --- including node-splitting, the net heat plot, and the design-by-treatment interaction test --- can determine *whether* inconsistency exists but provide limited ability to decompose it into independent sources or localize it to specific treatment comparisons.

**The contribution.** This paper introduces a sheaf-theoretic framework for NMA consistency analysis. Cellular sheaf theory, a branch of algebraic topology, provides a principled mathematical formalism for reasoning about local-to-global consistency on networks. The sheaf Laplacian, constructed from the precision-weighted coboundary map, captures all inconsistency information in a single matrix. Its eigendecomposition yields three novel capabilities:

1. **Decomposition:** The number of non-zero eigenvalues reveals how many independent sources of inconsistency exist in the network --- information unavailable from any existing NMA method.
2. **Localization:** Each eigenvector maps its corresponding inconsistency mode to specific treatments and comparisons, providing a spectral "fingerprint" of where inconsistency concentrates.
3. **Quantification:** Edge-level residual scores and the Global Inconsistency Index (GII) provide interpretable summaries at both local and global levels.

**Validation.** The method is validated on a simulated network with planted inconsistency (correctly localizes the biased edge as #1 ranked, normalized score 100%), consistent and inconsistent triangle tests, and the classic Hasselblad (1998) smoking cessation dataset. A comprehensive automated test suite (71 tests) verifies correctness of the linear algebra, sheaf engine, and visualization components.

**Implementation.** SheafNMA is the world's first browser-based tool for sheaf-theoretic NMA analysis. It is implemented as a single self-contained HTML file (2,550 lines, pure JavaScript), requiring no installation, server, or external dependencies. It exports R code for cross-validation with netmeta, and provides CSV, JSON, and LaTeX export. The tool is freely available under an open-access license.

**Why JRSS Series A.** This paper bridges algebraic topology and applied statistics in a way that is directly useful for evidence-based medicine. The sheaf-theoretic formalism is mathematically rigorous while yielding interpretable diagnostics for practitioners. JRSS Series A's tradition of publishing methodologically innovative work with applied relevance makes it an ideal venue.

**Novelty statement.** To the best of my knowledge, this is the first application of cellular sheaf theory to network meta-analysis or evidence synthesis. No existing software tool provides eigendecomposition-based inconsistency localization for NMA.

This manuscript has not been published or submitted elsewhere. The author declares no conflicts of interest.

Thank you for considering this submission.

Sincerely,

**Mahmood Ahmad**
Royal Free Hospital, London, UK
mahmood.ahmad2@nhs.net
ORCID: 0009-0003-7781-4478
