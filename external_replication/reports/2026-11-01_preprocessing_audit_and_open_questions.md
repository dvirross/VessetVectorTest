# Utah external replication: independent preprocessing audit and open questions

**Date:** 2026-11-01  
**Status:** Preliminary audit; methodological and inferential verification incomplete. **Do not interpret this note as validation of the external replication.**  
**Scope:** Preserve the existing manuscript, analytical data, simulations, and published results unchanged. This is a review record, not a revised analysis.

## Research concern

The Utah replication was originally implemented by automated assistants. The author raised the possibility that its strikingly stable results reflect a problematic preprocessing choice or overly mechanical transfer of decisions from the Marquette analysis. The crucial distinction is between (1) deterministic code reproducibility, (2) scientifically defensible inclusion and segmentation rules, and (3) valid statistical inference. Passing (1) does not establish (2) or (3).

## Findings independently checked against repository data

- Utah source `external_replication/data/raw/CrMcyclelength_share.csv`: **3,324 records, 581 women**. Reimplemented the specified gap-based run selection independently and compared every resulting row with `external_replication/data/sequences.csv`: **2,139/2,139 rows identical**. All recorded lengths matched the inclusive start–end date span. This checks the *implemented* preprocessing rule, not whether it is a good rule.
- Independently recounted pattern starts on the processed sequences using arithmetic definitions, obtaining **Haflaga 41; Dilug 76; Week 35; Week-Dilug 42; Dilug-in-Dilug 41**, the current counts. This does not validate the scientific definitions or the Monte Carlo results.
- The existing Utah pipeline keeps only the **single longest consecutive run with at least five cycles per woman**, breaking ties in favor of the earliest run. **Seven women have more than one eligible run**; consequently **37 cycles in otherwise eligible secondary runs are discarded**. Retaining these would add 37 analyzed cycles without adding women, if all other eligibility rules remain unchanged.
- Examples of qualifying additional runs (original cycle numbers): woman **377: 5–12 (8 cycles), then 14–19 (6)**, with #13 absent from the provided file; woman 89: 2–6 and 8–12 (both 5); 193: 2–6 and 8–12 (both 5); 286: 1–5 and 7–11 (both 5); 381: 3–10 (8) and 12–16 (5); 487: 5–10 (6) and 14–18 (5); 573: 4–10 (7) and 12–17 (6).
- Across the raw deposit, **968 usable recorded cycles belong to runs shorter than five cycles**, under the same date-gap segmentation. These are *not* 968 confirmed eligible observations; they highlight how consequential the eligibility rule is.
- The preprocessing implementation logged **170 date gaps among usable cycles**, including **169 with jumps in cycle numbers and one with no cycle-number jump**. The earlier broader descriptive report lists **186 gaps between consecutive raw records**, which uses a different denominator. Do not conflate the two.
- An illustrative date gap: woman's #377 cycle 12 ends **2011-04-29**, cycle 14 begins **2011-05-26**. Both missing cycle number and non-adjacent dates suggest an omitted record. The source does **not** establish the clinical reason or that the biological/halachic pattern ceased during that interval.
- The Utah report records two historical changes to counting definitions: Week anchor corrected on 2026-09-17; Week-Dilug anchors generalized, with a rerun on 2026-09-22. The original `reports/utah_replication_report.md` is marked *superseded* for affected numbers; rely on the current `results/summary.json` or `comparison.md` rather than obsolete figures.

## Source-data and interpretation concerns

The Utah deposit combines three CrM cohorts (CMFS, TTP, CEIBA). The bundled `external_replication/data/raw/Najmabadi_README20230824.txt` describes the CSV as a **derived dataset**, exported from Access, cleaned in SAS, and released for a different analysis. Its seven variables contain neither a reason for each missing cycle nor complete original diaries. A gap in the *released* sequence cannot automatically be interpreted as a physiological interruption or a follow-up dropout.

Marquette `scripts/filter_raw.py` removed duplicate subject records and subjects with <5 recorded cycles (plus handled a special duplicate), yielding 118 women / 1,554 cycles. Utah's extra *longest-run-only* rule is not equivalent. The original Marquette notebook `scripts/Filtering_original.ipynb` should also be reviewed for authorial judgement not fully captured by a terse code summary.

## Methodological options to test, not silently implement

1. **All qualifying contiguous runs:** retain *all* runs of length >=5, each as a separate countable interval; retain the original woman ID as the cluster key. Never count events across gaps. Compare observed counts, per-woman clustering, and null results with current longest-run-only baseline.
2. **Short-run sensitivity:** quantify support and pattern opportunities by length (minimum 2, 3, or 4 depending on pattern) rather than assume all runs shorter than five are useless. Declare a separate sensitivity estimand; avoid changing the primary inclusion rule after seeing results without labeling it post hoc.
3. **Null model definition:** permuting *within each observed run* conditions on run-specific multisets and is different from permuting *across all observed runs of the same woman* while preserving run lengths/boundaries. Choose and justify the estimand, and perform both as labeled robustness checks if appropriate. Treat multiple runs from one woman as dependent for per-woman analyses, leave-one-woman-out, and resampling.
4. **Missing cycles:** do not impute cycles into the primary exact-pattern count. If investigated, use transparent assumption-driven sensitivity or partial-identification bounds; actual missing-cycle lengths and continuation of a fixed pattern are unknown.
5. **Selection bias:** characterize how gap/short-run exclusion correlates with cohort, follow-up time, cycle lengths, and observed regularity, insofar as available. The observed file cannot independently identify the reason for each omitted cycle.

## Public supplementary data lead — not yet reconciled

A separately released CEIBA dataset is reported at https://hive.utah.edu/records/kfcm2-ynw54 (2025; reported 296 couples, 2,894 cycles, cycle-level detail and codebook). **The prior session could inspect descriptive documentation but could not retrieve/compare its actual Excel data. No IDs were matched, no missing cycles recovered, and no exclusion reasons confirmed.** Retrieve files if accessible and compare *without assuming IDs align*; date/length sequences may be needed for probabilistic linkage, with uncertainty clearly recorded. Check related published supplements and deposited source materials before contacting investigators.

## Outstanding validation

- Recover upstream CEIBA files and other cohort releases; determine whether any of the missing records are observable there and document provenance.
- Revisit the clinical and halachic semantics of counting patterns across gaps and what constitutes “fixed pattern persists” versus an observable event; seek author decision.
- Implement and document a **non-destructive** alternative-preprocessing sensitivity analysis on a separate branch, with unit tests for boundary handling, duplicate runs, woman IDs, and missing cycle numbers.
- Independently reproduce the permutation/null simulations, Mahalanobis test, Monte Carlo p-values, and comparison between the current and all-runs sensitivity datasets. **Not done in this review.**
- Do not update the manuscript's inferential claims until the audit is completed and alternative specifications are reviewed.

## Relevant files

- `scripts/filter_raw.py`, `scripts/Filtering_original.ipynb` (Marquette)
- `external_replication/scripts/prepare_utah.py`, `external_replication/data/raw/CrMcyclelength_share.csv`, `external_replication/data/sequences.csv`
- `external_replication/data/raw/Najmabadi_README20230824.txt`
- `scripts/patterns.py`, `external_replication/scripts/analysis_engine.py`, `external_replication/scripts/run_external.py`
- `external_replication/reports/utah_preprocessing_and_deviations.md`, `external_replication/results/provenance.json`

**Research question to carry forward:** What is the scientifically defensible observation unit and missing-data rule for detecting exact menstrual interval patterns when a woman has multiple disjoint observed runs? This must be resolved before asserting that the Utah replication confirms the paper.
