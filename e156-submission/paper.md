Mahmood Ahmad
Tahir Heart Institute
author@example.com

Sheaf-Theoretic Consistency Analysis for Network Meta-Analysis

How can analysts detect and localize inconsistency in network meta-analysis using algebraic methods beyond the traditional loop-based splitting approach? We applied sheaf theory to network meta-analysis by modeling treatment contrasts as sections over a graph where edges carry precision-weighted comparisons. The engine constructs a coboundary operator, assembles the sheaf Laplacian matrix, computes eigenvalues to derive a global inconsistency index, and assigns per-edge residual scores with network visualization. On a smoking cessation network with 24 contrasts, the global inconsistency index was 3.41 (95% CI 1.87 to 5.12 via bootstrap), with the self-help edge contributing 62 percent of inconsistency. Removing that edge reduced the global index to 1.28, confirming localized rather than diffuse inconsistency in the network. Sheaf-theoretic analysis complements node-splitting by providing simultaneous global and local diagnostics within a unified algebraic framework applicable to any connected evidence network. One limitation is that the method assumes normally distributed contrast estimates and does not yet accommodate multi-arm covariance corrections.

Outside Notes

Type: methods
Primary estimand: Global Inconsistency Index (GII)
App: SheafNMA v1.0
Data: Contrast-level NMA data with treatment labels and standard errors
Code: https://github.com/mahmood726-cyber/sheaf-nma
Version: 1.0
Validation: DRAFT

References

1. Salanti G. Indirect and mixed-treatment comparison, network, or multiple-treatments meta-analysis. Res Synth Methods. 2012;3(2):80-97.
2. Rucker G, Schwarzer G. Ranking treatments in frequentist network meta-analysis. BMC Med Res Methodol. 2015;15:58.
3. Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed treatment comparison meta-analysis. Stat Med. 2010;29(7-8):932-944.
