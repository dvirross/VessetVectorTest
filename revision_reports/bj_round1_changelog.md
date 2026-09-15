# Biometrical Journal pre-submission revision — change log (Stage 6)

Companion to `bj_round1_triage.md` (Stage 1 matrix). Reviews: R-Stat (ChatGPT, statistical), R-Fit (Gemini, journal fit), R-Lit (Perplexity, literature/novelty). Item numbers refer to the triage matrix.

## A. Issues resolved

| # | Issue | Resolution |
|---|---|---|
| 1 | Gap-bridging after cycle exclusions | Verified: `CycleNumber` is contiguous 1..n_i for all 118 women; stated in §3.1 and asserted in `scripts/patterns.load_data`. Documented that the raw→filtered script is not archived. |
| 2 | "Autocorrelation preserved under H_W" (wrong) | §2.3 rewritten: H_W preserves multiset and length, destroys ordering and hence serial dependence; H_W defined precisely in §3.3 as a conditional Monte Carlo test. |
| 3 | Over-interpretation of non-rejection | Abstract, §4.2, §5.1, §6 rewritten: heterogeneity is "consistent with / supported, not uniquely identified"; explicit list of what H_W does not show. SD comparison relabelled descriptive (definitional). |
| 4, 16 | Methodological novelty overclaimed; title | New title (null-model choice … Case Study); "methodological contribution" removed everywhere; statistic described as a routine adaptation of PERMANOVA/energy-type quadratic forms; article type Case Study stated on title page and in cover letter. |
| 5 | Same replicates for Σ̂ and reference | Split-batch D_M reported on real data (H_W: 2.95/.119 vs 2.93/.126) and size simulation shows plug-in ≡ split-batch (Table 4). |
| 6 | MDC section flawed; Table 5 mislabelled | Equation and table removed; all "would require a larger dataset" inferences removed. |
| 7 | No power/validation simulation | New §3.6 (design) and §4.6 (results): size under all three nulls; power vs AR(1) ρ∈{0,.25,.5,.75} and persistence q∈{0,.1,.2,.3}; new Fig 10. |
| 8 | Algorithm 1 one-sided | Fixed: two-sided formula with truncation at 1; identical to Eq. (1). |
| 9 | Fig 3 stale one-sided p (Dilug .0432) | Regenerated as PMF bars with two-sided p-values from current replicates; caption matches. |
| 10 | 7→5 family data-dependent; "pre-specified confirmatory" | §3.1 states plainly the family was fixed after zero counts were observed and why zero-count components cannot enter D_M; "pre-specified confirmatory" → "primary"; §5.2 notes nothing was pre-registered. |
| 11 | Unequal follow-up | Truncation to 12 cycles (joint p = .29, all z < 0), leave-one-woman-out range (.055–.224), share of events from ≥20-cycle women (35% of Haflaga) — §4.6, Table 3. |
| 12 | DiD deficit over-read | "Mechanism not identified; exploratory; Holm p = .19". |
| 13 | Abstract "correlations confirm association" | Removed; containment described as arithmetic, associations exploratory. |
| 14 | Vague generalisability | §5.6 gives the transferable lesson (null-model contrast) with concrete examples; no new claims. |
| 15 | Literature | Added Friedrich–Brunner–Pauly 2017, Lee & Braun 2012, Berrett et al. 2020, Hemerik & Goeman 2018, Phipson & Smyth 2010, Harlow et al. 2000, Fehring et al. 2006, Bull et al. 2019, Schmalenberger et al. 2021, Anderson 2001 (already in bib, now cited). |
| 20 | **Boundary bug**: stored nulls used 1,553 cycles | Fixed; all three nulls rerun (B = 50,000, seed 17, 1,554 cycles); every table, figure and in-text number updated. |
| 21 | Missing bib entry, stale README/notebook | `fehring2006` added; README rewritten; notebook rebuilt from `scripts/build_notebook.py` (all stale one-sided / "conservative lower bound" text gone; undefined variables gone). |
| 22 | Numeric references | natbib author–year (`plainnat`); all `\cite` → `\citep`/`\citet`. |
| — | **Table 2 mixed conventions** (self-identified) | Week/Week-Dilug/DiD p-values were one-sided (.133/.345/.543); now two-sided (.265/.688/1.00) with Holm column. |

## B. Issues rejected or not acted upon, and why

- **R-Fit: fit a Bayesian hierarchical state-space model instead.** Orthogonal to the question (exact rule-defined events, not prediction). Mentioned as future parametric null.
- **R-Lit citations with wrong attribution** (conditional permutation test attributed to Hemerik & Goeman; diary paper to "Cole"; CLD metric to Schmalenberger; "Drikvandi et al. 2013 Biometrical Journal"). Correct sources cited instead; the Drikvandi item not cited because its bibliographic details could not be reconciled offline.
- **R-Lit: Apple Women's Health Study 2023.** Relevant but details not verifiable offline; not added.
- **R-Stat optional: simulation comparison vs max-T / alternative covariance regularisation.** Qualitative comparison retained; a head-to-head power comparison of joint statistics is a methods-paper task inconsistent with the Case Study framing.
- **R-Stat optional: sensitivity to 18–54 d window.** Impossible with the archived filtered file; stated as a limitation.
- **R-Stat: include the two zero-count rules in the primary family.** Not done: a zero-count component makes Σ̂ singular and the rules need 6–9-cycle multi-value windows; the exclusion is disclosed as data-informed. Left as an author decision (see E).

## C. New analyses performed

1. **Corrected null rerun** — `scripts/rerun_nulls.py`, B = 50,000, seed 17, 1,554 cycles. Changes vs stored values are at the second decimal (e.g. H_W Haflaga mean 31.85→31.84; DiD p .040→.038, Holm .198→.192; H_W joint D_M 2.95→2.93, p .121→.126; H_G joint 5.10→5.08).
2. **Split-batch Mahalanobis** on real data under each null.
3. **Size study** (Part A): 10,000 fresh datasets × 3 nulls; marginal .029–.048, Holm .032–.037, joint .047–.053 (α=.05), .007–.011 (α=.01).
4. **Power study** (Part B): 8 scenarios × 400 datasets × 2,000 inner permutations. Joint power: AR(1) ρ=.25 → .51, ρ=.5 → .995, ρ=.75 → 1; persistence q=.1 → .98, q≥.2 → 1. Size checks at ρ=0/q=0: .06 (SE .012).
5. **Truncation to first 12 cycles** (Part C): joint D_M 2.49, p = .29; all z < 0.
6. **Leave-one-woman-out** (118 analyses): joint p ∈ [.055, .224]; DiD nominal p < .05 in 114/118, Holm-nonsignificant throughout.
7. **Follow-up influence**: ≥20-cycle women (n = 11) hold 20% of cycles, 35% of Haflaga events; Spearman(n_i, Haflaga events) = .28.
8. **Pairwise analyses recomputed** on corrected per-woman counts (Haflaga×DiD r = .27, p = .004, Holm .038; binary χ² 9.99, Holm .042).

## D. Claims weakened or changed

- "Methodological contribution … formal treatment in future work" → "routine adaptation of quadratic-form permutation tests; no novelty claimed".
- "Autocorrelation is preserved under H_W" → "H_W destroys ordering; a non-rejection says nothing about serial dependence in general".
- "Rules function primarily as empirical detectors of individual regularity" → "results are consistent with between-woman heterogeneity; supported, not uniquely identified".
- "Rules correctly identify women whose cycles are intrinsically stable" (SD comparison as validation) → descriptive; expected by construction.
- "All counts below MDC → larger dataset required" → removed; replaced by simulation-based statement of what would and would not have been detected (≈5 extra Haflaga events: 50%; ≥14: near-certain).
- "DiD deficit consistent with detecting genuine second-order regularity" → "direction noted; mechanism undetermined; exploratory".
- "Pre-specified confirmatory" → "primary; family fixed after zero counts observed; nothing pre-registered".
- Table 2 p-values for Week / Week-Dilug / DiD corrected from one-sided to two-sided.
- Haflaga exceedance "2 of 50,000" → "4 of 50,000" (new random stream, corrected data); H_G joint exceedance 0 → 6 of 50,000.
- Generalisability paragraph narrowed to one concrete transferable lesson.
- Cover letter: "three methodological contributions" removed; Case Study requested.

## E. Important unresolved risks

1. **Raw→filtered pipeline not archived.** Contiguity is verified, but a referee can still ask how 159→118 women and 1,649→1,554 cycles were reached. Author should locate the original filtering script (PhD repo?) or re-derive it from the Marquette archive before the Reproducible Research ZIP.
2. **Data-informed family (7→5).** Disclosed, but a strict referee may still object. Implementing the two Chozer Chalila rules under all three nulls (expected count ≈ 0) would close the point; needs the author's precise halachic definitions.
3. **Two bibliographic entries added offline** (`fehring2006` JOGNN 35(3):376–384; `schmalenberger2021` 123:104895) — verify volume/pages before submission. All other new entries are standard and confident.
4. **Single dataset, non-Jewish NFP population, probable overlap with Ecochard et al. cohort** — unchanged, disclosed.
5. **Power against smooth dependence.** ~50% at ρ = .25; weaker dependence is undetectable. Stated in §4.6/§5.7; a referee may still regard the non-rejection as weak evidence — the paper now says so itself.
6. **Length.** Body ≈ 4,900 words vs BJ's recommended 4,500 for a Case Study; Fig 5 and Fig 9 are candidates for a supplement.
7. **Legacy files** `sim_results.npz`/`strat_full_118.pkl` remain in the repo for provenance; a reader could load them by mistake. Documented in CLAUDE.md and README; consider moving to `legacy/`.

## F. Another independent review round?

**Yes — warranted and recommended**, for two reasons: (i) the reframing is substantial (title, positioning, primary/secondary structure, new results section), so a fresh reader is needed to judge whether the Case Study narrative is coherent rather than patched; (ii) the new simulation section has not been reviewed by anyone and its interpretation (what the non-rejection does and does not exclude) is now load-bearing.

## Assessment

**READY FOR NEW BLIND REVIEW.**

Rationale: the three issues that affected the validity of reported numbers (boundary bug, mixed one/two-sided conventions, MDC misinterpretation) are fixed and every number is regenerated from a single scripted pipeline; the two wrong statements (autocorrelation preserved; methodological novelty) are corrected; the interpretive claims are now bounded by a simulation study that quantifies sensitivity; and the article type matches the contribution. The open items in E are disclosure/provenance risks and author decisions, not defects that would make an independent evaluation premature.
