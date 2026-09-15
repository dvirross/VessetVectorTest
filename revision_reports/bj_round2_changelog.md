# Biometrical Journal Case Study — targeted pre-submission revision (round 2) — change log

Scope: implement only the items in `bj_round2_triage.md` (Groups 1–2). No study redesign, no new statistical methods, no new simulation families. Manuscript at start: `7192759`.

## A. Files changed
- `paper/vesset_stat.tex` — main manuscript (edits listed in B).
- `paper/vesset_stat_SI.tex` — **new** Supporting Information (exploratory pairwise analysis moved here).
- `paper/vesset_stat.pdf`, `paper/vesset_stat_SI.pdf` — recompiled.
- `paper/cover_letter.tex` / `.pdf` — one sentence noting the SI.
- `scripts/postprocess.py` — added the follow-up-stratified permutation and partial Spearman check for the two structurally expected pairs (written to `results/summary.json`), so the SI robustness numbers are reproducible from released code.
- `results/summary.json` — regenerated (new key `pairwise_followup_robustness`; all other values unchanged).
- `README.md` — SHA-256 of the filtered file; reproducibility-boundary statement.
- `CLAUDE.md` — file paths and status.
- `vesset_overleaf_main.zip` — rebuilt (now includes the SI source).

## B. Exact change log

| # | Issue | Location | Previous | Revised | Reason |
|---|---|---|---|---|---|
| 1 | "fourfold excess" numerically wrong | Conclusion | "driven by a fourfold excess of Haflaga events" | "driven by Haflaga events at 2.3 times their global-null expectation (26 observed versus 11.3 expected; z = 4.50)" | 26/11.31 = 2.30; 4.50 is the z-score. Only occurrence (searched whole text). |
| 2 | Age statement inconsistent and not describing the analysed sample | §3.1 Dataset; §5.7 Scope | §3.1: "Women were aged 21–43 years (M = 30.6, SD = 5.7)"; §5.7: "NFP users aged 18–42" | §3.1: "The source trial enrolled women aged 18–42 years (Fehring et al. 2013); in the analysed sample age is recorded for 99 of the 118 women and ranges from 21 to 43 years (M = 31.8, SD = 5.4)". §5.7: "(trial eligibility 18–42 years; Fehring et al. 2013)" | Verified from `data/FilteredData.csv`: Age present for 99/118 women, min 21, max 43, mean 31.8, SD 5.4. The previous M/SD (30.6/5.7) did not describe the analysed file. 18–42 is the trial's eligibility range as stated in the previous manuscript text citing Fehring et al. (2013). |
| 3 | 7→5 family: unsupported technical claims | §3.1 | "…a zero-count component cannot be included in a Mahalanobis statistic (its null covariance is degenerate at the observed sequence lengths), and adding two marginal tests with observed count 0 would not alter any conclusion below…" | "Of the seven rule types computable from cycle-length data, the two repeating-cycle rules … had no observed occurrences. The five-pattern family analysed here was therefore selected after inspection of the observed data and should be regarded as data-informed rather than pre-specified; the two zero-occurrence rules were not implemented in the analysis code and were not included in the multivariate analysis or in the Holm family." | Degeneracy and "no conclusion changes" were asserted without implemented analysis. History kept truthful (selection after inspection). §5.7 "were dropped" → "were not analysed after zero occurrences were observed". |
| 4 | "nominal size" overstated | Abstract; §4.5 Size; Conclusion | "confirms nominal size for all tests" / "Every procedure holds its nominal level" / "every test used holds its nominal size" | "type-I error at or below the nominal level for the marginal and Holm procedures and approximately nominal for the joint test" (all three places; §4.5 adds "conservative because the counts are discrete") | Table of empirical size shows marginals .029–.048. |
| 5 | Power statements mechanism-free | §4.5 final paragraph; §5.7 Statistical; Conclusion | "Had within-woman ordering produced an excess of roughly five Haflaga events…detected it about half the time; an excess of 14 or more would almost certainly have been detected." | "Under the AR(1) alternative at ρ = 0.25, which added about five Haflaga events on average, the joint test rejected in 51% of datasets; at ρ = 0.5 (about 14 additional events) it rejected in 99.5%. Under the persistence alternative, one exact repetition in ten cycles was detected in 98%…These figures do not transfer to arbitrary forms of temporal dependence." Limitations/Conclusion: "simulated exact-repetition mechanism…simulated whole-day AR(1)…under the simulated alternatives". | Power tied to the two simulated DGPs. |
| 6 | "The reason is…" | §4.2 | "The reason is that women with narrow cycle-length multisets…; conditioning on the multiset absorbs this." | "A principal contributor to this shift is that women with narrow cycle-length multisets…; conditioning on each woman's multiset preserves this feature, whereas the global nulls erase it." | Consistent with Discussion ("not uniquely identified"). |
| 7 | Pairwise analysis in main text | §3.5 (methods, ~600 words), §4.4 (results), Table 4, Fig. 9 | Full methods rationale, results, table and heatmap in main text | Moved to SI §S1 (methods and rationale, results, Table S1, Figure S1, follow-up robustness, interpretation). Main text keeps a one-paragraph methods pointer (§3.4) and a one-paragraph results summary (§4.4) with the headline numbers. | Secondary exploratory material; addresses R-Stat validity concern and R-Fit focus concern together. |
| 8 | Follow-up dependence of pairwise permutation null | SI §S1.2, §S1.4; main §3.4, §5.7 | Acknowledged qualitatively only | SI states the null is marginal, gives Spearman(nᵢ, rate) per pattern (.25 Haflaga … 0.00 Week-Dilug), partial r_s given nᵢ (.254 / .239) and follow-up-quintile-stratified permutation p (.0046 / .0097 vs unrestricted .0036 / .0113); labelled robustness diagnostic, not primary. | Reproduced by `scripts/postprocess.py` (seed 17, B = 50,000). |
| 9 | Abstract opened with the application | Abstract | Two opening sentences on Halacha | Opens with the statistical problem (exchangeability null for rule-defined counts in clustered sequences), then the Case Study application, then the specific rules | Case Study positioning; title unchanged. |
| 10 | Foundational randomisation references not cited | §1 ¶3 | Only clustered/conditional-permutation literature cited | Added "Randomisation tests derive their validity from the exchangeability structure they impose (Edgington & Onghena 2007; Good 2005)" | Entries were already in the .bib; support the sentence they are attached to. |
| 11 | Reproducibility boundary | Data Availability Statement; README | "cached replicates…allow every reported number to be reproduced" | Two paragraphs: what is reproducible (with SHA-256 of the filtered file), what is not (preprocessing script not preserved), the two inclusion criteria, contiguity, and an explicit statement that pre-filtering counts are unavailable | Transparency; no reconstruction attempted. |
| 12 | Power table half-rounding vs code output | Table (power) | +5.0, .19, .67, +24.2, .17 | +5.1, .18, .66, +24.1, .16 | Aligned with `numpy` rounding of the stored results so the table is reproduced exactly by the code. |
| 13 | Algorithm 1 unreferenced | §3.5 Simulation design | — | "(Algorithm 1)" pointer added | Cross-reference completeness. |

## C. Analyses changed
**No primary statistical analysis was changed.** Null replicates, joint tests, marginal tests, Holm adjustments, simulation study, truncation and leave-one-out are untouched. One robustness diagnostic (partial Spearman given nᵢ; permutation stratified by follow-up quintile, seed 17, B = 50,000) was added to `postprocess.py` for the exploratory pairwise analysis in the SI; it does not enter the main text's inference.

## D. Supporting Information
Moved from main: §3.5 pairwise methods/rationale; §4.4 pairwise results; Table 4 (pairwise Spearman); Figure 9 (heatmap). Added: §S1.4 follow-up robustness (values above); statement that the binary check addresses metric choice, not follow-up dependence. Main-text numbering after the move: Figures 1–10 (heatmap removed; 3-panel scatter is Fig. 8, power Fig. 9, heaping Fig. 10), Tables 1–5.

## E. Consistency audit
Checked programmatically against `results/summary.json` and `results/sim_validation.npz`: 118 women / 1,554 cycles (11 and 7 occurrences, consistent); 5–45 cycles; 1,289 cycles after truncation; all H_G/H_iid/H_W means, SDs, z, two-sided p, Holm p; D_M 5.08/4.87/2.93 with exceedances 6/17/6,285; split-batch 2.95/.119; 3-pattern 4.50/4.39/2.68 (p .067); LOO .055–.224; SD comparison; influence 20%/35%; pairwise .27/.004/.038 and .23/.011/.102; binary 9.99/.042, 7.26/.118; size ranges .029–.048, Holm .032–.037, joint .048–.053 (.007–.011); power table cell-by-cell; truncation results. Terminology: `\HG`/`\HI`/`\HW` macros used throughout (23/17/33 occurrences); no "stratified permutation" or "bootstrap" misuse for H_W; "fourfold", "nominal size", "degenerate", "would not alter", "The reason is", 21–43 alone, 30.6/5.7 all absent. All `\ref` targets exist; no dangling references to the moved section. Float order sequential. Main text compiles with 0 errors, 0 overfull boxes (31 pages); SI compiles with 0 errors (5 pages). Body word count ≈ 4,800 (Introduction–Conclusion, excluding floats).

Corrections made as a result of the audit: item 12 (five half-rounded power cells) and item 13.

Unresolved uncertainties: the trial eligibility range 18–42 is taken from the previous manuscript text and its citation of Fehring et al. (2013); the source paper could not be fetched from this environment to re-verify it.

## F. Remaining known limitations
- Preprocessing from the raw Marquette archive is not reproducible (script not preserved; pre-filtering counts unavailable). Disclosed in §3.1 and the Data Availability Statement.
- The two zero-occurrence rules are not implemented; the five-pattern family is data-informed. Disclosed.
- Single non-Jewish NFP cohort, probably overlapping Ecochard et al. (2024). Disclosed.
- Power is characterised only for the two simulated mechanisms. Disclosed.
- Body length ≈ 4,800 words vs the 4,500 recommended for a Case Study.

## G. Submission readiness
**READY TO SUBMIT TO BIOMETRICAL JOURNAL** — with one author glance recommended (not blocking): confirm that "18–42 years" is the enrolment eligibility stated in Fehring et al. (2013).
