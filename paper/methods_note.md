# E156 Methods Note — SheafNMA v0.2

**E156 Entry:** #153 (update — supersedes v0.1 synthetic-network framing)
**Type:** methods
**Primary Estimand:** Global Inconsistency Index (GII)

## 156-word body

Can treating NMA contrasts as algebraic graph sections yield a global inconsistency statistic and per-edge localization from one object? We applied cellular-sheaf theory to thirteen published Cochrane and textbook NMAs from the R netmeta package, spanning diabetes, thrombolytics, acupuncture, and antimanic agents. The engine builds a precision-weighted coboundary operator and sheaf Laplacian, solves a weighted least-squares system, and reports the Global Inconsistency Index and per-edge residual scores. Sheaf-residual flagging identified localized inconsistency in five of thirteen NMAs where the design-by-treatment chi-squared failed to reject at 0.05, indicating 38.5 percent of apparently-consistent networks harbored at least one inconsistent edge. Across a 36-cell simulation grid with one thousand replications, Type-I error averaged zero and sensitivity reached 0.981 at planted-inconsistency magnitude 0.5. Per-edge agreement with node-splitting was 45.4 percent, reflecting a deliberately conservative Youden-optimal threshold of 2.69. Sheaf analysis emits global and local diagnostics from one algebraic object, complementing node-splitting, though single-dimensional stalks do not yet handle multi-arm covariance.
