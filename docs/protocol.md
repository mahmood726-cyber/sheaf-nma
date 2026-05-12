        # Sheaf Nma: Protocol Registration

        **Author:** Mahmood Ahmad, Royal Free Hospital, London, UK
        **ORCID:** 0009-0003-7781-4478
        **Registration Date:** 2026-03-28 (v0.1) · 2026-05-12 (v0.2 amendment)
        **Current version:** v0.2.0
        **Repository:** https://github.com/mahmood726-cyber/sheaf-nma.git

        ## Objective

        SheafNMA applies sheaf-theoretic topology to network meta-analysis, detecting and localising inconsistency in treatment networks using higher-order algebraic structures. v0.2 validates the method on real published NMAs against canonical inconsistency tests.

        ## Methods

        This project employs a deterministic, TruthCert-certified pipeline for evidence synthesis. All analytical choices are pre-specified in this protocol document prior to data analysis. The implementation uses versioned code with fixed random seeds where applicable, structured input schemas, and automated validation against reference outputs. Provenance is recorded via hash-linked evidence locators in accordance with the TruthCert framework. Statistical methods follow established meta-analytic guidelines (Cochrane Handbook, PRISMA 2020) and are validated against reference implementations in R (metafor, meta, netmeta) or Python with tolerance ≤ 1×10⁻⁶.

        ## v0.2 amendment

        - Real-NMA-corpus validation against 13 `netmeta` built-in datasets (exported via `corpus/export_netmeta.R`).
        - 36-cell × 1000-rep simulation grid with held-out Youden threshold (`τ_inc = 0.5` training cell at `n_treatments=6, n_studies=5`; frozen threshold 2.6912).
        - Comparators: Higgins (2012) DBT-χ², Bucher (1997) loop closure, Dias (2010) node-splitting via `netmeta::netsplit`.
        - Multi-arm covariance correction and higher-dimensional stalks deferred to v0.3.

        ## Availability

        - Code: https://github.com/mahmood726-cyber/sheaf-nma.git (tag `v0.2.0`)
        - Dashboard: https://mahmood726-cyber.github.io/sheaf-nma/
        - Manuscript draft: `paper/sheafnma_v2.md` (Research Synthesis Methods format)
        - E156 #153 methods note: `paper/methods_note.md` (156 words, 7 sentences)
        - Implementation spec / plan: `docs/superpowers/specs/2026-05-12-*.md`, `docs/superpowers/plans/2026-05-12-*.md`

        ---

        **AI Disclosure Statement**

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI is used as a constrained synthesis engine operating on structured inputs and predefined rules, rather than as an autonomous author. Deterministic components of the pipeline, together with versioned, reproducible evidence capsules (TruthCert), are designed to support transparent and auditable outputs. All results and text were reviewed and verified by the author, who takes full responsibility for the content. The workflow operationalises key transparency and reporting principles consistent with CONSORT-AI/SPIRIT-AI, including explicit input specification, predefined schemas, logged human-AI interaction, and reproducible outputs.
