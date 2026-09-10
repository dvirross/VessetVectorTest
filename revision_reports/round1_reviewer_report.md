# Round 1 Reviewer Report — SMMR Reviewer Simulation
**Paper**: "Are Halachic Vesset Patterns Overrepresented in Empirical Menstrual Cycle Data? A Permutation Test"  
**Date**: 2026-09-10  
**Status**: Post-revision review (acting as SMMR Reviewer 2)

---

## Summary Decision

**Major Revision** — The revised manuscript represents substantial improvement over the original submission. The authors have addressed most reviewer concerns, added a power analysis, implemented permutation-based dependence tests, and provided the CONSORT flow diagram. However, several issues remain that require attention before the manuscript can be accepted.

---

## Point-by-Point Assessment

### Previously Raised Concerns: Status

**SMMR 2.1 (Terminology "stratified randomisation")**: ✅ RESOLVED  
The text now consistently uses "within-woman permutation null" throughout, correctly describing a shuffle without replacement. The previous revision incorrectly changed "stratified randomisation" to "bootstrap (sampling with replacement)," which was factually wrong. The current revision has been corrected and the implementation (shuffle without replacement) now matches the description.

**SMMR 2.2 (Markov limitation)**: ✅ RESOLVED  
The permutation null's conservativeness relative to a Markov chain model is now clearly articulated.

**SMMR 2.3 (Dependence structure)**: ✅ RESOLVED  
The polynomial hierarchy and mathematical dependencies are well-explained.

**SMMR 2.4 (Mahalanobis formula)**: ✅ RESOLVED  
The paper now correctly describes $D_M^2 = (\mathbf{O}-\boldsymbol{\mu})^\top \hat{\Sigma}^{-1} (\mathbf{O}-\boldsymbol{\mu})$ and its geometric interpretation.

**SMMR 2.5 (Covariance matrix estimation)**: ✅ RESOLVED  
Condition numbers (4.4–5.4) are reported; the paragraph in §3.4 describes empirical estimation from simulation.

**SMMR 2.6 (Pseudocode/algorithm)**: ⚠️ PARTIALLY ADDRESSED  
A prose description of Σ̂ estimation is present, but no numbered algorithm block has been provided. For reproducibility, a formal algorithm box (with pseudocode for the full simulation pipeline: permutation, count, covariance estimation, test statistic computation) would be preferable.

**SMMR 2.7 (CONSORT flow)**: ✅ RESOLVED  
Figure 2 (CONSORT-style participant flow) is now included at §3.1.

**SMMR 2.8 (Power analysis)**: ✅ RESOLVED  
Table 5 (power analysis under $H_W$) is added with MDC computed for each pattern. The formula MDC = μ_W + 2.487σ_W is clearly stated.

**SMMR 2.9 (Chi-square: small cells + confound)**: ✅ RESOLVED  
Table 3 now includes a $p_{\text{perm}}$ column from 50,000 permutations per pair, which:  
- Eliminates the small expected-cell problem (no asymptotic approximation)  
- Controls the cycle-count confound by preserving column marginals

---

### New Issues in This Revision

**R-NEW-1 (Minor): MDC formula missing factor explanation**  
Table 5 caption reports MDC but doesn't explain where 2.487 comes from in the table or caption itself. The formula is in the text, but adding "(z_0.05 + z_0.80) = 1.645 + 0.842 = 2.487" in the caption footnote or table note would make Table 5 self-contained for readers who read tables independently of the text.

**R-NEW-2 (Minor): §4.3 "Results" paragraph lacks z-score table cross-reference**  
The §4.3 results paragraph mentions Dilug-in-Dilug z=−2.07 but this value appears in Table 2, not in the text for the other patterns. Readers comparing Table 2 and Table 5 will note that z values differ slightly (e.g., Haflaga: −1.25 in Table 2 vs. −1.24 in Table 5; Week: −0.47 in both). A sentence acknowledging that Tables 2 and 5 use slightly different simulation runs (different seeds/runs: Table 2 was the original simulation, Table 5 is a fresh 50K run) would prevent confusion.

**R-NEW-3 (Minor): Missing F3 element**  
The manuscript mentions (§3.4/discussion) that the Mahalanobis test is more powerful than applying Hotelling's T² or Bonferroni-corrected marginal tests. However, no quantitative comparison with Hotelling T² is given. A brief paragraph or footnote noting the relationship (under Gaussian approximation, D_M² ~ χ² with df=5; Hotelling T² is an exact multivariate test requiring the distributional assumption; the permutation-based D_M is distribution-free) would sharpen the methodological claim.

**R-NEW-4 (Moderate): Sentence about Week pattern in §5.3 is unclear**  
The paragraph on Week in §5.3 discusses weekly biological rhythms. With the corrected weekly anchor set {22, 29, 36, 43, 50} (haflaga values), corresponding to cycle lengths of 21, 28, 35, 42, 49 days, the "anchoring" to multiples of 7 is clear. However, the text says "anchored to multiples of 7" without specifying which value of the anchor is biologically most important. Since cycle lengths ≥35 days would be unusually long (beyond the normal range), nearly all Week patterns in this dataset come from L=21 or L=28 cycles. A sentence clarifying that the Week pattern predominantly captures 28-day cycles in this dataset would help.

**R-NEW-5 (Moderate): Relabelling invariance paragraph inconsistency**  
The §5.6 relabelling invariance paragraph says "we verified that all pattern counts are identical whether computed on H_k or directly on L_k with the anchor adjusted accordingly (Week-Dilug at L=29 instead of H=30)." This claim is correct: using H_k with Week-Dilug anchor=30 gives 18, and using L_k with anchor=29 gives 18 (verified computationally). However, the paragraph does NOT mention that using L_k with anchor=30 (unadjusted) gives only 12, which would be misleading. The phrase "with the anchor adjusted accordingly" is important and should be emphasised more strongly.

**R-NEW-6 (Minor): Table 3 column widths may overflow**  
The updated Table 3 now has 7 columns with specified widths. With the added $p_{\text{perm}}$ column (C{1.2cm}), the total table width may exceed text width in the SMMR format. Authors should verify that the table compiles without overflow warnings.

**R-NEW-7 (Editorial): "Within-woman null" in abstract still says "bootstrap"**  
The keywords section (line ~118) still includes "multinomial sampling" but the abstract now says "permutation null" for H_W. The keywords should include "within-woman permutation" as a distinct keyword to help indexing.

---

### Remaining Statistical Concerns

**S-1 (Moderate): Power analysis — MDC exceeds observed for ALL patterns**  
Table 5 shows that all five observed counts fall below their MDC, confirming individual underpowerment. However, the interpretation in the text states "the composite Mahalanobis test achieves significance by pooling correlated evidence." This is slightly misleading because the Mahalanobis test is significant under the GLOBAL null (p<.001) but NOT under the within-woman null (p=.121). The power discussion should clarify: MDC is computed under H_W (the within-woman permutation null), and under H_W even the joint test is non-significant. The claim that the Mahalanobis test "achieves significance" in the power analysis paragraph appears to refer to the global null result — this should be made explicit.

**S-2 (Minor): p=.004 for Haflaga×DiD in Table 3 is very close to the usual threshold for claiming robustness after permutation testing**  
The paper claims that the permutation p-values "confirm" the Fisher's exact findings. While true numerically, a sentence explicitly noting that the permutation test p-value (0.004) is below 0.01, providing strong evidence even after accounting for distributional assumptions, would strengthen the claim.

---

## Recommendation

**Major Revision required** with focus on:
1. (**Required**) Clarify the power analysis paragraph: distinguish global-null significance from within-woman-null non-significance (S-1)
2. (**Required**) Add algorithm pseudocode for reproducibility (SMMR 2.6)
3. (**Recommended**) Add F3: brief Hotelling T² comparison paragraph (R-NEW-3)
4. (**Minor**) Address R-NEW-1, R-NEW-2, R-NEW-4, R-NEW-5, R-NEW-7

The fundamental statistical analyses are sound, and the methodological contributions (polynomial hierarchy, distance-based joint test, within-woman permutation null, permutation chi-square) are well-motivated and correctly implemented.
