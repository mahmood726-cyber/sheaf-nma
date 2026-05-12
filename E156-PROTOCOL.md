# E156 Protocol — `sheaf-nma`

This repository is the source code and dashboard backing an E156 micro-paper on the [E156 Student Board](https://mahmood726-cyber.github.io/e156/students.html).

---

## `[153]` Sheaf-Theoretic Consistency Analysis for Network Meta-Analysis

**Type:** methods  |  ESTIMAND: Global Inconsistency Index (GII)  
**Data:** Contrast-level NMA data with treatment labels and standard errors

### 156-word body

Can treating NMA contrasts as algebraic graph sections yield a global inconsistency statistic and per-edge localization from one object? We applied cellular-sheaf theory to thirteen published Cochrane and textbook NMAs from the R netmeta package, spanning diabetes, thrombolytics, acupuncture, and antimanic agents. The engine builds a precision-weighted coboundary operator and sheaf Laplacian, solves a weighted least-squares system, and reports the Global Inconsistency Index and per-edge residual scores. Sheaf-residual flagging identified localized inconsistency in five of thirteen NMAs where the design-by-treatment chi-squared failed to reject at 0.05, indicating 38.5 percent of apparently-consistent networks harbored at least one inconsistent edge. Across a 36-cell simulation grid with one thousand replications, Type-I error averaged zero and sensitivity reached 0.981 at planted-inconsistency magnitude 0.5. Per-edge agreement with node-splitting was 45.4 percent, reflecting a deliberately conservative Youden-optimal threshold of 2.69. Sheaf analysis emits global and local diagnostics from one algebraic object, complementing node-splitting, though single-dimensional stalks do not yet handle multi-arm covariance.

### Submission metadata

```
Corresponding author: Mahmood Ahmad <mahmood.ahmad2@nhs.net>
ORCID: 0000-0001-9107-3704
Affiliation: Tahir Heart Institute, Rabwah, Pakistan

Links:
  Code:      https://github.com/mahmood726-cyber/sheaf-nma
  Protocol:  https://github.com/mahmood726-cyber/sheaf-nma/blob/main/E156-PROTOCOL.md
  Dashboard: https://mahmood726-cyber.github.io/sheaf-nma/

References (topic pack: network meta-analysis):
  1. Rücker G. 2012. Network meta-analysis, electrical networks and graph theory. Res Synth Methods. 3(4):312-324. doi:10.1002/jrsm.1058
  2. Lu G, Ades AE. 2006. Assessing evidence inconsistency in mixed treatment comparisons. J Am Stat Assoc. 101(474):447-459. doi:10.1198/016214505000001302

Data availability: No patient-level data used. Analysis derived exclusively
  from publicly available aggregate records. All source identifiers are in
  the protocol document linked above.

Ethics: Not required. Study uses only publicly available aggregate data; no
  human participants; no patient-identifiable information; no individual-
  participant data. No institutional review board approval sought or required
  under standard research-ethics guidelines for secondary methodological
  research on published literature.

Funding: None.

Competing interests: MA serves on the editorial board of Synthēsis (the
  target journal); MA had no role in editorial decisions on this
  manuscript, which was handled by an independent editor of the journal.

Author contributions (CRediT):
  [STUDENT REWRITER, first author] — Writing – original draft, Writing –
    review & editing, Validation.
  [SUPERVISING FACULTY, last/senior author] — Supervision, Validation,
    Writing – review & editing.
  Mahmood Ahmad (middle author, NOT first or last) — Conceptualization,
    Methodology, Software, Data curation, Formal analysis, Resources.

AI disclosure: Computational tooling (including AI-assisted coding via
  Claude Code [Anthropic]) was used to develop analysis scripts and assist
  with data extraction. The final manuscript was human-written, reviewed,
  and approved by the author; the submitted text is not AI-generated. All
  quantitative claims were verified against source data; cross-validation
  was performed where applicable. The author retains full responsibility for
  the final content.

Preprint: Not preprinted.

Reporting checklist: PRISMA 2020 (methods-paper variant — reports on review corpus).

Target journal: ◆ Synthēsis (https://www.synthesis-medicine.org/index.php/journal)
  Section: Methods Note — submit the 156-word E156 body verbatim as the main text.
  The journal caps main text at ≤400 words; E156's 156-word, 7-sentence
  contract sits well inside that ceiling. Do NOT pad to 400 — the
  micro-paper length is the point of the format.

Manuscript license: CC-BY-4.0.
Code license: MIT.

SUBMITTED: [ ]
```


---

_Auto-generated from the workbook by `C:/E156/scripts/create_missing_protocols.py`. If something is wrong, edit `rewrite-workbook.txt` and re-run the script — it will overwrite this file via the GitHub API._
