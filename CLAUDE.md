# CLAUDE.md — VessetVectorTest

Statistical analysis of halachic vesset patterns in menstrual cycle data.
Paper targeting: **Biometrical Journal (Wiley, Q1) — Case Study article type** (decided 2026-09-15 after BJ pre-submission triage; see `revision_reports/bj_round1_triage.md`).

> **Before every session:** ask the user what they want to do. Do not start a series of tasks without explicit instruction. Confirm intent before any chain of edits.

---

## 1. Plan before touching anything

Before writing code or editing LaTeX:
- State what you intend to change and why
- If it affects simulation results, flag it — the rerun is now fast (~1 min for the three nulls; ~20 min for the validation study on 4 cores), but every number in the paper derives from it
- If it affects the paper narrative (primary vs secondary result, Case Study framing), confirm with user first

**Pre-computed results live in `results/`:**
- `results/null_replicates.npz` — H_G, H_iid, H_W replicate count matrices (seed 17, B = 50,000, **1,554 cycles**, 118 women) + per-woman H_W counts (`within_by_woman`) for leave-one-out
- `results/sim_validation.npz` — size calibration (Part A), power vs AR(1)/persistence (Part B), truncation to 12 cycles (Part C)
- `results/summary.json` — every reported statistic, produced by `scripts/postprocess.py`

**Legacy files (do not use for the paper):** `sim_results.npz` and `strat_full_118.pkl` were computed with a boundary bug that dropped the last woman's final cycle (1,553 cycles). They are kept only for provenance. `boot_arr` inside `sim_results.npz` is an old within-woman bootstrap, not H_W.

**Code layout:**
```
scripts/filter_raw.py      archived source file -> data/FilteredData.csv (no 18–54 filter: archive already in range;
                           removes 3 duplicated-record women, women with <5 cycles, 8 duplicate rows of nfp8107)
scripts/patterns.py        counting rules (vectorised + reference loop, equivalence test), null generators,
                           two-sided Monte Carlo p, Holm, Mahalanobis test, load_data (asserts CycleNumber contiguity)
scripts/hcc_rule.py        Haflaga Chozer Chalila (block of >=3 values repeated 3x; window >=9): loop + vectorised
                           counter, equivalence test, 3 nulls B=50,000 seed 17 -> results/hcc_null.npz (SI S3)
scripts/hcc_joint.py       regenerates the 3 null streams (same seed/call order as rerun_nulls.py, verified identical)
                           + HCC as 6th component -> results/hcc_joint.npz; 6-pattern D_M unchanged, p: <.001/<.001/.150
scripts/per_woman_hw.py    H_W within regularity strata (SD terciles/halves) + per-woman concentration checks, from the
                           cached per-woman replicates; run(DATA, NPZ) is called by both notebooks; -> results/per_woman_hw*.json (SI S6)
scripts/rerun_nulls.py     -> results/null_replicates.npz
scripts/sim_validation.py  -> results/sim_validation.npz
scripts/postprocess.py     -> results/summary.json (prints all paper numbers)
scripts/build_notebook.py  -> notebooks/analysis.ipynb (then execute with nbconvert to regenerate figures)
```
Any change to counting logic goes in `scripts/patterns.py` only; the notebook imports it.
External replication on the Utah Creighton Model cohort (2026-09-17) is self-contained in `external_replication/` (see its README).

---

## 2. Git workflow

- Work on the branch designated for the session (currently `claude/gallant-pascal-5zy75j`); earlier sessions used `main`
- Commit and push after every few changes — token limits can truncate sessions
- Commit message format: brief imperative summary + the attribution footer given in the session's system reminder

---

## 3. Biometrical Journal requirements (target journal)

**Submission portal**: https://authors.wiley.com/journal/BIMJ
**Article type**: **Case Study** (4500 words recommended). Do not re-frame as a methodological Research Article: three independent pre-submission reviews agreed the Mahalanobis-from-permutation statistic is an adaptation of known quadratic-form tests, not a new method.
**Peer review**: Single-blind. Submission does **not** need to be anonymous; the `vesset_stat_submission.tex` two-file workflow is retired.

| Requirement | Status |
|---|---|
| No footnotes | ✅ |
| ≤5 keywords, alphabetical | ✅ exchangeability, Halacha, menstrual cycle, permutation test, within-subject null |
| Unstructured abstract | ✅ |
| Data Availability / CoI / Ethics statements | ✅ |
| Author–year (APA-like) references | ✅ `natbib[authoryear,round,sort]` + `plainnat` (`apalike.bst` is not natbib-compatible — do not switch to it) |
| Acknowledgments (AI disclosure + funding) on the **title page**; Patient Consent a separate `\section*` after Ethics | ✅ end matter: Ethics, Patient Consent, Conflict of Interest, Data Availability, References — no Funding section, no end-matter Acknowledgments |
| PDF metadata (title/author) | ✅ `pdftitle`/`pdfauthor` in main and SI |
| ORCID on title page | ✅ |
| Full postal addresses (both affiliations) + corresponding address (Shenkar) | ✅ title page and cover letter |
| Reproducible Research ZIP | prepare at revision: `scripts/`, `results/`, `data/`, notebook |
| Figures as separate files | ✅ PDF + PNG in `paper/figures/` |
| Complete list of figure legends in text | ⬜ add at revision if requested |

Cover letter (`paper/cover_letter.tex`) is written for the Case Study framing and includes all BJ-required declarations.

---

## 4. Notebook: keep in sync with the paper

`notebooks/analysis.ipynb` is generated by `scripts/build_notebook.py` — **edit the builder, not the notebook**, then (this regenerates every figure except Fig 2, which is produced separately by `paper/figures/gen_consort_flow.py`):
```bash
python scripts/build_notebook.py
(cd notebooks && jupyter nbconvert --to notebook --execute --inplace analysis.ipynb)
```

| Paper figure | Generator |
|---|---|
| fig1_cycle_distribution | notebook |
| fig2_consort_flow | `paper/figures/gen_consort_flow.py` (standalone, run from repo root; writes next to itself via `Path(__file__).parent`; source box cites Fehring 2012 dataset + Fehring et al. 2013 trial) |
| fig3_null_distributions | notebook — PMF bars, **two-sided** p labels |
| fig4_zscores | notebook |
| fig5_cdfs | notebook — now **Figure S4 in the SI** (Section S5); no longer in the main text |
| fig6_stratified_comparison | notebook — H_W two-sided p per panel computed from data |
| fig7_dm_nulldist | notebook |
| fig8_chisq_heatmap | notebook — now **Figure S1 in the SI** (values computed, not hardcoded) |
| fig9_joint_3panels | notebook (exploratory 3-pattern subspace) |
| fig10_power | notebook (from `results/sim_validation.npz`) |
| figA1_heaping_diagnostic | notebook — frequency of H values |
| figE2_null_distributions, figE5_comparison_z | `external_replication/scripts/external_figures.py`; copied to `paper/figures/` for SI Section S4 (Figures S2, S3) |

---

## 5. LaTeX: always compile and check

```bash
cd paper
pdflatex -interaction=nonstopmode vesset_stat.tex
bibtex vesset_stat
pdflatex -interaction=nonstopmode vesset_stat.tex
pdflatex -interaction=nonstopmode vesset_stat.tex
grep "^!" vesset_stat.log        # must be empty
```
Non-fatal warnings to ignore: `titlesec` "entered in horizontal mode", `rerunfilecheck`.
Use `\citep{}` / `\citet{}` (author–year), never bare `\cite{}` for parenthetical citations.
Compile the SI separately: `pdflatex vesset_stat_SI.tex` ×2 (no BibTeX needed).

**Bibliography conventions (final pass 2026-09-16):**
- `fehring2012data` = the archived dataset (Fehring, 2012, Marquette e-Publications) — cite for *the data*.
- `fehring2013` = the randomised trial article (Fehring, Schneider, Raviele, Rodriguez & Pruszynski, *Contraception* 88(1):24–30) — cite for *how the data were collected*. Never attribute the 18–54-day filter to Creinin et al.
- `ross2026statistical` (and `anon2026statistical`) is the B.D.D. article, published April 2026 in the triple issue *B.D.D* 38–40, pp. 59–74 (end page inferred from the catalogue contents: the next article starts at p. 75), `note = {In Hebrew}`; formerly the `@unpublished` key `ross2021statistical` (accepted 3 September 2021). No DOI is known; do not invent one.
- `ross2022phd`/`anon2022phd` exact title: *Probabilistic and Statistical Analysis of the Menstrual Cycle with regard to Jewish Religious Laws* (not "in a Halachic Context").
- `king2020` (J. Econometrics 219(2):425–455) replaced the old `king2018` working paper; currently uncited. Bibliography DOI pass done 2026-09-20 (13 DOIs added; Münster co-author **Lone** Schmidt; Harlow first name **Siobán**).
- `ross_planned` uses `year = {in preparation}` (renders "Ross, in preparation").
- Ecochard et al. 2024 co-author is **Marie** Schneider; Münster 1992 first author is **Kirstine** Münster.
- Prose uses `$p$-value` (math-mode p), British spelling throughout.

---

## 6. Terminology — get these right

| Wrong | Correct |
|---|---|
| autocorrelation is "preserved" under H_W | H_W preserves each woman's multiset and n_i; it **destroys** ordering and hence serial dependence |
| H_W shows no temporal structure / no serial dependence | for these exact count statistics, the observed ordering is not unusual relative to random reorderings conditional on each woman's multiset |
| rules "function primarily as detectors of individual regularity" | results are **consistent with** between-woman heterogeneity; this explanation is supported, not uniquely identified |
| Haflaga women's lower SD "validates" the rule | it is expected by construction (Haflaga is defined by equal intervals); descriptive only |
| stratified permutation / bootstrap (for H_W) | within-woman permutation null (without replacement; conditional Monte Carlo test) |
| "pre-specified confirmatory" | "primary" (no dated analysis plan exists; 5-pattern family was fixed after two rules had zero counts) |
| methodological contribution / new test | adaptation of a standard quadratic-form permutation statistic (PERMANOVA / energy-distance lineage) |
| MDC = μ_W + 2.802σ_W, "would require a larger dataset" | **removed**; sensitivity is quantified by the simulation power study |
| Dilug "significant" | Dilug exploratory (H_G p = .085, Holm .34) |
| DiD deficit "detects genuine second-order regularity" | direction noted; mechanism undetermined; Holm p = .19 |
| Bonferroni assumes independence | Holm valid under any dependence |
| establish / prove / demonstrate | suggest / consistent with / provide no evidence that |
| one-sided p in any table/figure | all marginal p are two-sided: min{1, 2·min((b⁺+1)/(B+1), (b⁻+1)/(B+1))} |
| Week-Dilug "= 30" | Week-Dilug anchors are **all** H ≡ 2 (mod 7) = {23, 30, 37, …} (one-weekday progression; `WEEK_DILUG_ANCHORS` in `patterns.py`). H = 30 is the dominant anchor in practice (all 18 Marquette events; 40 of 42 Utah events) — say so, but never define the rule by 30 alone |
| "haflaga value(s)" for the quantity H_k | **interval(s)** H_k = L_k + 1 (defined once in §2.1 as the halachic *haflaga*; "haflaga convention" for the +1 shift is fine). Capitalised **Haflaga** = the pattern (three equal intervals) only |

---

## 7. Current numbers (seed 17, B = 50,000, 1,554 cycles; **Week-Dilug anchors {7n+2}**, 2026-09-22) — from `results/summary.json`

**H_W (within-woman) — primary:**

| Pattern | Obs | μ_W | σ_W | z | p (two-sided) | Holm |
|---|---|---|---|---|---|---|
| Haflaga | 26 | 31.84 | 4.68 | −1.25 | .254 | 1 |
| Dilug | 61 | 67.04 | 7.38 | −0.82 | .460 | 1 |
| Week | 32 | 33.79 | 3.88 | −0.46 | .749 | 1 |
| Week-Dilug | 18 | 21.46 | 3.30 | −1.05 | .370 | 1 |
| Dilug-in-Dilug | 27 | 39.12 | 5.84 | −2.07 | .038 | .192 |
| Joint 5-pattern | D_M = 2.97, p = .115 (split-batch 2.99, .107); LOO range .051–.207 |

**H_G (global permutation):** Haflaga z = 4.50, two-sided p = .0002 (4 of 50,000 ≥ 26), Holm .001; Dilug z = 1.84, p = .085 (Holm .34); Week p = .265; Week-Dilug z = 0.42, p = .774; DiD p = 1.0. Joint D_M = 5.07, p < .001 (7 exceed). H_iid: joint D_M = 4.87 (16 exceed); Week p = .373; Week-Dilug p = .817.
**3-pattern exploratory:** D_M = 4.50 / 4.39 / 2.72 (p = .059 under H_W).
**Pairwise (rate Spearman):** Haflaga×DiD r = .27, p_perm = .004, Holm .038; Haflaga×Week-Dilug r = .23, p = .011, Holm .102; binary check χ² = 9.99 (Holm .042), 7.26 (Holm .118).
**SD comparison:** Haflaga women n = 20, 1.54 (0.59) vs no-pattern n = 35, 3.17 (1.75); U = 133, Welch t = −5.05 df 45.5, CI [−2.29, −0.98], d = −1.13.
**Influence:** 11 women with ≥ 20 cycles hold 20% of cycles and 35% of Haflaga events.
**Validation study:** see `results/sim_validation.npz` and Section 4.6 of the paper.

---

## 8. Self-citations — keep them low

- `ross2022phd` ≤ 5 occurrences; `ross2026statistical` ≤ 2. Describe methods inline rather than citing the thesis.

---

## 9. Key file paths

```
data/RawData.csv               author's copy of the Marquette archive file (1,665 cycles, 159 women, 80 columns)
data/FilteredData.csv          1,554 cycles, 118 women; CycleNumber contiguous within every woman
results/                       cached Monte Carlo output (see §1)
scripts/                       all analysis code (see §1)
paper/vesset_stat.tex          manuscript (Case Study framing, author–year refs)
paper/vesset_stat_SI.tex       Supporting Information (exploratory pairwise analysis; compile separately)
paper/references.bib           bibliography
paper/cover_letter.tex         BJ cover letter + declarations
paper/supplementary_for_review/  Ross_2022_PhD_dissertation_Hebrew.pdf, Ross_2026_BDD_38-40_article_Hebrew.pdf (author's copy of the article published in B.D.D 38–40, 2026, pp. 59–74; review only)
paper/figures/                 figures as PDF + PNG
notebooks/analysis.ipynb       generated notebook; produces fig1, fig3–fig10 and figA1 (Fig 2 comes from paper/figures/gen_consort_flow.py)
revision_reports/              SMMR rounds 1–3; bj_round1_*.md; bj_round2_triage.md; bj_round2_changelog.md
external_replication/          Utah Creighton Model replication (self-contained: data, scripts, results, reports); do not edit from the paper side
```

Raw Marquette archive (not reachable from the remote sandbox): https://epublications.marquette.edu/data_nfp/7/ — the author's copy is `data/RawData.csv` (also in https://github.com/dvirross/PhD with the original `Filtering.ipynb`); `scripts/filter_raw.py` reproduces the filtered file from it. **Provenance wording (2026-09-16):** never say an 18–54-day inclusion criterion was applied by us — the archive already lies in that range; say the script exists.

---

## 10. Status (2026-09-16, BJ pre-submission rounds 1–2 + final editorial pass)

Responded to three AI pre-submission reviews (statistical, journal-fit, literature). Key changes: Case Study reframing and new title; H_W/autocorrelation misstatement corrected; MDC section removed and replaced by simulation validation (size, power, truncation, LOO); Algorithm 1 and Fig 3 made two-sided; Table 2 one-sided values corrected; boundary bug fixed and all nulls rerun on 1,554 cycles; literature on permutation for clustered data and menstrual variability added; claims softened. Full list in `revision_reports/bj_round1_changelog.md`. **Round 2 (targeted):** fourfold→2.3-fold; age statement corrected (trial eligibility 18–42; analysed sample 99/118 with age, 21–43, M 31.8 SD 5.4); 7→5 wording made truthful (no degeneracy claim); size wording 'at/below nominal'; power statements DGP-specific; pairwise analysis moved to SI with follow-up-stratified robustness; abstract opens with the statistical problem; Edgington/Good cited; reproducibility boundary + SHA-256 in DAS. Verdict: READY TO SUBMIT (`bj_round2_changelog.md`). **Final editorial pass (2026-09-16):** bibliography corrected (Fehring 2012 dataset vs Fehring et al. 2013 trial, 5 authors + correct title; Marie Schneider; Kirstine Münster; Ross in preparation); dataset provenance made consistent in Author Note, §2.3, §3.1, Fig 2 (script + caption), Ethics, DAS, SI, cover letter, README; Creinin attribution of the filter removed; Funding and Patient Consent statements added; exact title with "Jewish" everywhere incl. PDF metadata; `natbib sort`; `$p$-value` typography; cover letter's duplicate declarations page removed; `gen_consort_flow.py` now writes `fig2_consort_flow.*`. All three PDFs compile with 0 errors/0 overfull/0 undefined. **Final compliance passes (2026-09-16):** AI disclosure finalised with model names and moved, together with the funding statement, to the title-page Acknowledgments; AI-assisted-code sentence lives in §3.2 (Pattern Detection); `fehring2012data` note reduced to `Dataset (.sav/.csv)`; `gen_consort_flow.py` made portable; DAS says the notebook and the Fig 2 script together produce all figures. Literature-search caveat in §1 replaced by a reference to a documented search (author's dissertation-era search + 12 logged web queries, 2026-09-16, none found): log in `revision_reports/literature_search_log.md`, summary in SI Section S2. Keep "we are not aware of"; do not claim absence. **External replication incorporated (2026-09-20):** Utah Creighton Model cohort (Stanford & Najmabadi 2023, doi:10.7278/S50d-4gxs-s4hj; Najmabadi et al. 2020) — 270 women, 2,139 cycles; counts (41, 76, 35, 42, 41); D_M 10.04 / 9.56 / 1.22 (p .916) [values after the 2026-09-22 Week-Dilug redefinition; at incorporation they were (…40…), 9.99 / 9.53 / 1.20]; Haflaga z +8.33 (H_G), +1.13 (H_W); decomposition gap 6.28 vs 6.76, residual −1.79 vs +1.57. Added: one abstract sentence, §3.6 (Methods), §4.6 + Table (Results), a decomposition paragraph in §5.1, shared limitations in Scope, Ethics + DAS sentences, bib entries `stanford2023data`/`najmabadi2020`, cover-letter sentence; SI Section S4 (Tables S4–S11, Figures S2–S3) and Section S5 (former Fig 5 = Figure S4). Fig 5 removed from the main text. No Fehring number changed; `scripts/patterns.py` and `external_replication/` untouched. **Week-Dilug redefinition (2026-09-22):** anchors generalised from H = 30 to {7n+2} (one-weekday progression; k-weekday progressions would collapse the rule into a two-cycle Haflaga). Marquette observed counts unchanged (all 18 events at 30); Utah Week-Dilug 40 → 42 (2 events at 37). All nulls, validation, HCC joint, per-woman checks and the Utah pipeline rerun; every dependent number updated (Marquette joint D_M 5.07/4.87/2.97, p_W .115; Utah 10.04/9.56/1.22, p_W .916). `external_replication/reports/*.md` carry a superseded-numbers note. **Per-woman checks (2026-09-22):** H_W holds within regularity terciles in both cohorts (all joint p ≥ .20 except Marquette above-median half p = .047, driven by a Dilug *deficit*, Holm .033, not replicated in Utah); women with most events are those H_W expects (Spearman .61 / .70); concentration statistics at null (top-6 share p .11 / .74; Gini p .55 / .74). §4.2, §4.6, §5.1 sentences; SI Section S6 (Tables S12–S13); notebook cells added to both notebooks.

**Outstanding (author decisions):**
- [x] `fehring2006` (Fehring, Schneider & Raviele; JOGNN 35(3):376–384; DOI 10.1111/j.1552-6909.2006.00051.x; PMID 16700687) and `schmalenberger2021` (10 authors as listed; Psychoneuroendocrinology 123:104895; DOI 10.1016/j.psyneuen.2020.104895; PMID 33113391) verified via web search snippets (Crossref/PubMed themselves are blocked from the sandbox)
- [x] Zero-count rule evaluated (2026-09-17): only **one** further rule is interval-computable — *Haflaga Chozer Chalila* (Dilug Chozer Chalila depends on calendar dates → 17 rules = 11 date-dependent + 6 interval-computable). Observed 0; null mean 0.0003 (H_G), 0.0004 (H_iid), 0.028 (H_W); P(≥1) = 0.03%, 0.04%, 2.8%; two-sided p = 1. 6-component joint test on paired replicates: D_M unchanged (5.07/4.87/2.98), p .0005/.0007/.139 (vs .0002/.0003/.115). Reported in §3.1 and SI Section S3 (Tables S2–S3); conclusions unaffected
- [x] Trial eligibility confirmed (dissertation p. 22 + ClinicalTrials.gov NCT00843336): women 18–42 years, cycle lengths 21–42 days, no hormonal contraception in prior 3 months; NCT number now cited in §3.1 and Fig 2
- [x] Patient Consent wording now cites the Marquette archive record: "anonymized for dissemination", "data reuse was agreed to by subjects in the consent form" (data_nfp/7 "Rights, Permissions, & Privacy" field, confirmed verbatim by the author on 2026-09-16). This is reported in the Patient Consent Statement only — do **not** restore the quotation to the `fehring2012data` bib note, which stays `note = {Dataset (.sav/.csv)}`
- [ ] Word count exceeds BJ's recommended 4,500 for a Case Study (≈6,100 Introduction–Conclusion words by `scratchpad/wc.py`-style count, floats excluded; the earlier ≈4,800 figure used a different, unrecorded method). Fig 5 has been moved to the SI (Figure S4); further cuts are an author decision
- [x] Raw→filtered pipeline recovered: `data/RawData.csv` (author's copy of the archive, 1,665 cycles/159 women) + `scripts/Filtering_original.ipynb` (from https://github.com/dvirross/PhD) + `scripts/filter_raw.py` (re-implementation; verified identical output). Dissertation and the B.D.D. article (accepted version; published 2026 in vol. 38–40) are in `paper/supplementary_for_review/` (Hebrew; upload as "supplementary material for review only")
