# Biometrical Journal Case Study — pre-submission triage, round 2 (triage only; no edits made)

Manuscript: `paper/vesset_stat.tex` at commit `7192759` (32 pp.).
Reviews: **R-Stat** = ChatGPT independent referee report (Major Revision; "not this exact version"); **R-Fit** = Gemini editorial/referee review (Minor Revision; "submit now"); **R-Lit** = Perplexity literature/novelty audit, two passes (both "No, reframe first"; the second pass is largely a plan and a set of clarifying questions rather than an audit).

Every disputed statement was checked against the current LaTeX source by line number. One numerical check was run on the data solely to judge the severity of R-Stat's pairwise-permutation objection (nothing was added to the manuscript, code or results).

---

## Stage 1–3. Consolidated triage matrix

| ID | Consolidated issue | Raised by | Exact reviewer concern | Location (current ms.) | Category | Validity | Evidence from current manuscript | Severity | Desk review? | Stat. validity? | Interpretation? | Must fix? | Recommended action | Effort | Conf. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 | Pairwise permutation p-values do not condition on follow-up | R-Stat B1 | Both rates depend on nᵢ, so women are not identically distributed; permuting one rate vector is "not a valid randomization test"; the Holm-significant Haflaga×DiD may be anti-conservative; binary check has the same flaw | §3.5 (ll. 453–510), §4.4, Table 4, Fig. 9 | Statistical validity / Exploratory | **PARTLY VALID** | The permutation null is *marginal* independence of the two rates; it does not condition on nᵢ, and §3.5 says so ("does not remove it entirely"). Check on data: Spearman(nᵢ, Haflaga rate) = .25, (nᵢ, DiD rate) = .10; partial r_s(Haf×DiD ∣ nᵢ) = .254 vs raw .269; permutation **stratified by follow-up quintile** p = .0047 vs unrestricted .0036 (Holm ≈ .05). The confound exists in principle and is immaterial in fact. "Must be fixed" is overstated; "make descriptive only" is unnecessary. | MINOR (as a reported result); MODERATE (as an unqualified claim) | No | Marginally | Marginally | **NO** (but see G2) | One sentence stating that the permutation tests marginal independence, that nᵢ correlates .25 with the Haflaga rate, and that a follow-up-stratified permutation gives essentially the same p (.005). Optionally demote the whole section to SI (C2). | Trivial wording (+ one already-computed number) | High |
| C2 | Pairwise section dilutes focus / length | R-Fit #1, E, H; R-Stat G | Move §3.5 and §4.4 to Supporting Information to sharpen the null-model lesson and respect the 4,500-word guide | §3.5, §4.4, Table 4, Fig. 9 | Structure / focus | **OPTIONAL IMPROVEMENT** | Body ≈ 4,900 words; §3.5 alone is ~600. The section is explicitly labelled exploratory and secondary. Moving it costs nothing inferentially. | MINOR | Slightly (length) | No | No | NO | Move §3.5 + §4.4 + Table 4 + Fig. 9 to SI, leaving a 3-sentence pointer in Results. Resolves C1 and C2 together. | Structural revision (low) | Medium |
| C3 | 7→5 family: degeneracy justification not auditable | R-Stat B2 | Code implements only five rules; the claims that zero-count components make Σ̂ degenerate and that seven-rule Holm "would not alter any conclusion" are asserted, not shown | §3.1 ll. 296–304; §5.2; DAS | Reproducibility / Multiple testing | **PARTLY VALID** | l. 301: "(its null covariance is degenerate at the observed sequence lengths)" and l. 302: "adding two marginal tests with observed count 0 would not alter any conclusion" — neither is backed by code or results. The *disclosure* is adequate (Rule G); the *technical justification* is not verified. | MODERATE | Possibly (reproducibility emphasis) | No (5-pattern inference unaffected) | No | **YES — wording**; implementation NO | Reword to what is known: the two rules were not implemented; their observed counts are zero by inspection; they were excluded before the null analyses; consequences for the five-pattern tests: none. Delete the degeneracy claim unless demonstrated. Implementing the rules under the three nulls is an author decision (needs exact halachic definitions). | Trivial wording; implementation = Code/data verification (optional) | High |
| C4 | Raw→filtered pipeline not reproducible | R-Stat B3 | 118 women/1,554 cycles cannot be regenerated from the public archive with released code; contiguity check covers one failure mode only | §3.1 ll. 289–294; Fig. 2; DAS | Reproducibility / Data | **VALID (already disclosed)** | l. 292: "The script that produced the filtered file from the raw archive was not preserved". The gap is real and the manuscript admits it; the archive is not reachable from this environment to reconstruct. | MODERATE | Possibly | No | No | **UNCERTAIN** | Minimum: add SHA-256 of `data/FilteredData.csv`, the deterministic inclusion rules as a numbered list, and step counts (raw women/cycles → filtered) once the author confirms the raw totals. Full reconstruction: author task, Group 3. | Code/data verification (author) | High |
| C5 | Age-range inconsistency | R-Stat B4a | Methods "21–43" vs Limitations "18–42" | l. 280 vs l. 1023 | Reporting / consistency | **VALID** | Confirmed: l. 280 "Women were aged 21--43 years" (computed from the analysed file); l. 1023 "NFP users aged 18--42" (the trial's eligibility range from Fehring et al.). Different quantities, presented as if the same. | MODERATE (credibility) | Yes (easy target) | No | No | **YES** | State both explicitly once: eligibility 18–42 at enrolment (Fehring 2013); observed 21–43 in the analysed sample. Author to confirm the Age variable's timing. | Trivial wording | High |
| C6 | "fourfold excess" is wrong | R-Stat B4b | 26 vs 11.31 is 2.3×, not 4×; 4.50 is the z-score | Conclusion l. 1059 | Reporting / consistency | **VALID** | Confirmed: "driven by a fourfold excess of Haflaga events". 26/11.31 = 2.30. | MODERATE (headline number) | Yes | No | Yes | **YES** | "…driven by Haflaga events at 2.3 times their global-null expectation (26 vs 11.3; z = 4.50)". | Trivial wording | High |
| C7 | "Nominal size for all tests" overstated | R-Stat D1, E | Marginal tests reject 2.9–4.8% at α = .05 — conservative, not nominal | Abstract l. 118; §4.6 l. 809; Conclusion l. 1062 | Reporting / Simulation | **VALID (minor)** | Table 4 itself shows .029–.048 for marginals; text says "confirms nominal size for all tests" / "holds its nominal size". Literal contradiction with own table. | MINOR | Slightly | No | Slightly | **YES** (trivial) | "type-I error controlled at or below nominal for the marginal and Holm procedures and approximately nominal for the joint test" in all three places. | Trivial wording | High |
| C8 | Five-event power statement generalises beyond the simulated DGP | R-Stat D2 | 51% is specific to rounded Gaussian AR(1) at ρ = .25, not to any mechanism adding five Haflaga events | §4.6 ll. 891–896 | Interpretation / Simulation | **PARTLY VALID** | l. 891–893: "Had within-woman ordering produced an excess of roughly five Haflaga events…the joint test would have detected it about half the time" — the following sentence already limits to "whole-day autoregression or exact repetition", but the headline sentence is mechanism-free. | MINOR | No | No | Yes (mild) | NO (G2) | Tie to the DGP: "Under the AR(1) alternative at ρ = 0.25, which added about five Haflaga events on average, the joint test rejected in 51% of datasets." | Trivial wording | High |
| C9 | §4.2 "The reason is…" too definite | R-Stat D3 | Conditioning on the multiset explains the shift; which feature of the multiset does is not decomposed | §4.2 ll. 676–680 | Interpretation | **PARTLY VALID** | l. 678: "The reason is that women with narrow cycle-length multisets generate many exact repetitions under any ordering". Discussion (ll. 924–929) is appropriately hedged; this one sentence is not. | MINOR | No | No | Yes (mild) | NO (G2) | "A principal contributor is that…" | Trivial wording | High |
| C10 | Re-centre title/abstract on the statistical lesson | R-Lit R1 (Must fix), R-Lit R2 (Must fix) | Title/abstract emphasise Halacha; BJ will desk-reject a domain analysis; "BJ has no Case Study section" | Title; Abstract ll. 101–123; §1 | BJ fit / Structure | **PARTLY VALID; one premise INCORRECT** | Title already leads with "Null-Model Choice for Rule-Defined Pattern Counts in Clustered Longitudinal Sequences"; §1 ¶3 states the general problem and §5.6 the transferable lesson. Abstract, however, opens with two sentences on Halacha before the statistical question. R-Lit's claim that BJ lacks a Case Study article type is false (BJ author guidelines list Case Studies; R-Stat quotes them). | MINOR–MODERATE (positioning) | Yes (positioning) | No | No | NO (G2) | Reorder the abstract so the first sentence poses the statistical question (null-model choice for dependent rule-defined counts in clustered sequences) and the Halachic rules enter as the application. No change to title. | Trivial wording | Medium |
| C11 | Missing citations on exchangeability under constraints / randomisation tests | R-Lit R2 (Good 2005); R-Lit R1 ("permutation tests for dependent data [arxiv]", unspecified) | Ground the conditional nulls in established resampling theory | §1 ¶3; §3.3 | Literature | **PARTLY VALID** | `good2005` and `edgington2007` are in the .bib but are **no longer cited** in the current text (dropped in the rewrite); Berrett et al., Hemerik & Goeman, Friedrich et al., Lee & Braun, Winkler, Nichols are cited. The arXiv item is unidentifiable — cannot act. | MINOR | No | No | No | NO (G2) | Add `\citep{edgington2007,good2005}` at the first mention of randomisation tests (§1 ¶3 or §3.3). | Trivial | High |
| C12 | Ecochard/Fehring cohort overlap should be framed as a data-origin caveat, not speculation | R-Lit R1, R2 | Provide a "Data Origin" contrast; acknowledge as potential correlation without asserting identity | §2.3 ll. 250–265 | Data / Literature | **ALREADY ADEQUATELY ADDRESSED** | ll. 255–258 give both cohort sizes, say the samples are "probably drawn from the same data-collection infrastructure and should not be regarded as independent samples", and draw the methodological consequence for H_W. That is the requested framing. | NONE | No | No | No | NO | None (optionally retitle the paragraph "Data origin and cohort overlap"). | No action | High |
| C13 | "Algebraic proximity" phrasing borders on a generative claim | R-Fit #2 | Rephrase to stress it is a consequence of the definitions, not a biological expectation | §2.2 ll. 229–234; §4.4 | Interpretation | **ALREADY ADEQUATELY ADDRESSED** | l. 234: "makes positive co-occurrence with Haflaga plausible but does not entail it"; §3.5 last paragraph and §4.4 repeat that attribution would need a generative model. | NONE–MINOR | No | No | No | NO | Optional: add "by the definitions alone" to l. 233. | Trivial | High |
| C14 | Case-Study novelty / desk-rejection risk | R-Stat C4–5 ("moderately convincing, not decisive"; fit *adequate*), R-Fit C ("5/5 unusually strong"), R-Lit ("high risk") | Standard methods on one old niche dataset; lesson already familiar to statisticians | Whole paper | Case-study contribution / BJ fit | **OPTIONAL — editorial judgment, not a defect** | Manuscript does not claim methodological novelty (§1 ¶4, §5.4); contribution stated as the null-model contrast + calibration. R-Stat: "not a standard application of a standard method to a standard endpoint". Cannot be reduced by analysis; only by positioning (C10) and cover letter. | — | Yes (irreducible residual) | No | No | NO | C10 + ensure cover letter leads with the lesson (it already does). Accept residual risk. | No action beyond C10 | Medium |
| C15 | Interpretation of H_W non-rejection and of heterogeneity | (none — R-Stat D and R-Fit G explicitly *clear* the manuscript) | — | Abstract l. 121–123; §3.3 l. 374; §5.1 ll. 924–929; §6 ll. 1066–1068 | Interpretation | **NO ISSUE (confirmed)** | Abstract: "consistent with…provides no evidence that…does not establish the absence of serial dependence"; §5.1: "simplest account consistent with both…supported by the data; it is not uniquely identified by them"; §6 same. Rule A/B satisfied in every section. | NONE | No | No | No | NO | Retain wording. | No action | High |

---

## Stage 5. Cross-reviewer synthesis

### A. Issues raised independently by two or more reviewers
- **The pairwise section (C1 + C2).** R-Stat objects on validity grounds, R-Fit on focus/length grounds. These are *different* arguments about the *same* section; agreement raises confidence that the section is the weakest part of the paper, not that R-Stat's specific validity claim is decisive. After checking: the validity concern is real in principle and negligible in practice (partial and stratified results above); the focus concern is legitimate. One action (demote to SI with a one-sentence stratified-permutation note) resolves both.
- **Positioning for desk review (C10/C14).** R-Lit (both passes) and, more mildly, R-Stat raise it; R-Fit rates fit 5/5. Agreement is partial and the strongest form (R-Lit) rests on a false premise about BJ's article types. Valid only to the extent the abstract's opening two sentences lead with religion rather than statistics.
- No two reviewers identified the same **valid MAJOR** issue.

### B. Single-reviewer issues worth taking seriously
- **C5 age inconsistency** and **C6 "fourfold"** (R-Stat) — both verified against the source; small, but exactly the kind of error that costs credibility at desk review.
- **C3 unverified degeneracy claim** (R-Stat) — the manuscript asserts a technical property (degenerate covariance) it has not demonstrated; that sentence should not survive to submission.
- **C7 "nominal size"** (R-Stat) — the text contradicts its own Table 4.
- **C4 provenance** (R-Stat) — already disclosed, but a checksum and rule list are cheap and directly address BJ's reproducibility emphasis.

### C. False positives / already solved / unnecessary
- **C12** Ecochard overlap framing — §2.3 already does exactly what R-Lit asks.
- **C13** proximity phrasing — already "plausible but does not entail".
- **C15** H_W / heterogeneity over-interpretation — both R-Stat and R-Fit explicitly confirm the current wording is appropriate; there is nothing to fix.
- **R-Lit: "BJ has no Case Study section"** — incorrect; BJ's author guidelines list Case Studies (R-Stat quotes them).
- **R-Lit: "first quantitative analysis" claim weakened by Ecochard** — the manuscript claims only "the first formal randomisation-based test of the halachic rules on empirical data of which we are aware" (§1, §6); Ecochard et al. do not test halachic rules.
- **R-Lit second pass** — mostly a work plan and clarifying questions; it contains no additional verifiable criticism.
- **R-Stat B1 remedy "make the section descriptive only / drop p-values"** — disproportionate given the stratified check; a one-sentence qualification suffices.
- **R-Stat B2 remedy "implement both omitted rules under all three nulls"** — reasonable as a supplement but not required for the validity of any reported result; it is an author decision needing exact halachic definitions.

### D. Reviewer disagreements
- **Submit now?** R-Fit: yes. R-Stat: no, after four surgical fixes. R-Lit: no, after reframing. R-Stat's position is best supported: it identifies concrete, verified defects (C3, C5, C6, C7) that R-Fit missed, while R-Lit's reframing demand overshoots what the current title and §1 already do and rests partly on a wrong premise.
- **Pairwise section.** R-Stat: repair the inference; R-Fit: move to SI. Given the stratified check, R-Fit's structural remedy plus a one-line qualification is the better-supported and cheaper path.
- **Case-Study fit.** R-Fit 5/5 vs R-Lit "high desk-rejection risk" vs R-Stat "adequate, 60% to external review". R-Stat's middle position is the most credible: the fit risk is real, irreducible by analysis, and only marginally movable by abstract positioning.

---

## Stage 6. Pre-submission priority list

### GROUP 1 — MUST FIX BEFORE SUBMISSION
1. **C6** "fourfold excess" → factual error in the Conclusion's headline number.
2. **C5** Age range 21–43 (Methods) vs 18–42 (Limitations) → factual inconsistency in the sample description.
3. **C3 (wording only)** Remove the unverified "null covariance is degenerate" claim and the unverified "would not alter any conclusion" claim; state plainly what was and was not done with the two omitted rules.
4. **C7** "nominal size for all tests" (Abstract, §4.6, Conclusion) contradicts Table 4 → "at or below nominal (marginal/Holm); approximately nominal (joint)".

### GROUP 2 — WORTH FIXING IF LOW-COST
5. **C1 + C2** Add one sentence to the pairwise methods/results noting the marginal-independence null, Spearman(nᵢ, Haflaga rate) = .25, and the follow-up-stratified permutation p = .005; move §3.5/§4.4/Table 4/Fig. 9 to Supporting Information with a short pointer.
6. **C8** Tie the "five extra Haflaga events → 51%" sentence to the AR(1) ρ = .25 mechanism.
7. **C9** §4.2 "The reason is" → "A principal contributor is".
8. **C10** Reorder the abstract's opening so the statistical question comes first; title unchanged.
9. **C11** Cite Edgington & Onghena (2007) and Good (2005) at the first mention of randomisation tests (entries already in .bib).
10. **C4 (minimum form)** Add SHA-256 of `data/FilteredData.csv`, a numbered list of inclusion rules, and an explicit "cannot be regenerated from released code" sentence to the Data Availability Statement; step counts once the author confirms raw totals.

### GROUP 3 — DO NOT DELAY SUBMISSION FOR THIS
- Implementing the two Chozer Chalila rules under all three nulls (C3, full form) — author decision; needs definitions.
- Reconstructing the raw→filtered script (C4, full form) — author task; archive not reachable here.
- Alternative statistics (max-T, parametric/hierarchical null), further simulation families — R-Stat explicitly says no further simulation is needed.
- Retitling §2.3 "Data origin" (C12), adding "by the definitions alone" (C13).
- Generalisability to observant Jewish women — disclosed; replication is the answer.
- Perplexity's unspecified "arXiv permutation tests for dependent data" — unidentifiable.

---

## Stage 7. Stopping-rule assessment
1. Unresolved FATAL issue? **NO.**
2. Unresolved MAJOR statistical-validity issue? **NO.** The only validity objection (C1) is confirmed immaterial (partial r_s .254 vs .269; stratified p .0047 vs .0036) and concerns an exploratory result.
3. Unresolved MAJOR interpretation issue? **NO.** Two reviewers independently confirm the H_W/heterogeneity wording is bounded; remaining items (C8, C9) are single sentences.
4. Unresolved MAJOR literature/novelty issue? **NO.** Two classic references dropped in the rewrite should be restored; no prior work duplicates the contribution.
5. Unresolved MAJOR BJ Case-Study-fit issue? **NO** as a defect; a residual editorial risk remains that analysis cannot remove (R-Stat: ~60% to external review).
6. Two or more reviewers, same valid MAJOR issue? **NO.**
7. Remaining criticisms are predominantly **minor and optional**, with four verified reporting errors (Group 1) that are trivial to fix.
8. **ONE TARGETED REVISION THEN SUBMIT.** A further blind LLM round would most likely re-litigate positioning and request additional robustness analyses; it is unlikely to find another verifiable defect of the C5/C6 kind now that three independent reports have been checked line by line.

---

## Stage 8. Final triage verdict

**B. NEARLY READY.** Four specific, verified reporting corrections (Group 1) should be made; the Group 2 items are cheap and worth doing in the same pass; no broad review cycle is warranted.

### Minimal pre-submission action plan (do not implement yet)
1. **"fourfold excess"** — Conclusion, l. 1059 — replace with "at 2.3 times their global-null expectation (26 vs 11.3; z = 4.50)". *Manuscript text.*
2. **Age range** — §3.1 l. 280 and §5.7 l. 1023 — state enrolment eligibility 18–42 (Fehring et al.) and observed 21–43 in the analysed sample in §3.1; make l. 1023 refer to eligibility. *Manuscript text; author confirms Age variable.*
3. **Degeneracy claim** — §3.1 ll. 299–304 — delete "(its null covariance is degenerate…)" and "would not alter any conclusion below"; replace with: the two rules were not implemented; zero occurrences established by inspection; excluded before any null analysis; the five-pattern tests are unaffected by the exclusion. *Manuscript text.*
4. **"Nominal size"** — Abstract l. 118, §4.6 l. 809, Conclusion l. 1062 — "type-I error at or below nominal for marginal and Holm procedures, approximately nominal for the joint test". *Manuscript text.*
5. **Pairwise section** — §3.5/§4.4 — add the marginal-independence/stratified-permutation sentence (numbers above); move section, Table 4 and Fig. 9 to SI with a pointer. *Manuscript text (structural); one number from an existing check.*
6. **Power sentence** — §4.6 l. 891 — tie to AR(1) ρ = .25. *Text.*
7. **"The reason is"** — §4.2 l. 678 — "A principal contributor is". *Text.*
8. **Abstract opening** — ll. 101–107 — lead with the statistical question. *Text.*
9. **Citations** — §1 ¶3 or §3.3 — add `edgington2007`, `good2005`. *Literature (entries exist).*
10. **Provenance** — Data Availability Statement — SHA-256 of the filtered file, numbered inclusion rules, explicit non-regenerability statement; step counts pending author's raw totals. *Text + one checksum.*
