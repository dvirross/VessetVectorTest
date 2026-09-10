# Round 1 Author Response
**Paper**: "Are Halachic Vesset Patterns Overrepresented in Empirical Menstrual Cycle Data? A Permutation Test"  
**Date**: 2026-09-10  
**Re**: Response to SMMR Reviewer 2 Round 1 Report

---

We thank the reviewer for their thorough and constructive assessment of the manuscript. The feedback has substantially improved the paper. We respond point-by-point below, indicating the location of all changes in the revised manuscript.

---

## Point-by-Point Response

### Resolved Issues from Previous Round (confirmed)

**SMMR 2.1–2.5, 2.7–2.9**: All previously resolved issues (terminology, Markov limitation, dependence structure, Mahalanobis formula, covariance estimation, CONSORT flow, power analysis) remain resolved and unchanged. We confirm no regression in these areas.

---

### Required Changes (Now Addressed)

**SMMR 2.6 (Pseudocode/algorithm) — Previously Partially Addressed; Now Fully Resolved**

We have added a formal Algorithm block (Algorithm 1) titled "Simulation pipeline" to §3.5 (Statistical Inference). The algorithm pseudocode covers all four stages the reviewer requested:
1. Observed count computation
2. Permutation/simulation loop (with explicit per-null-model specification of the randomisation operation)
3. Covariance matrix estimation from simulated vectors
4. Mahalanobis statistic computation and empirical p-value computation

The algorithm distinguishes clearly between the three null models (H_G: global shuffle; H_iid: iid multinomial draw; H_W: per-woman Fisher–Yates shuffle without replacement). The $(b+1)/(B+1)$ p-value convention is included in the algorithm.

---

**S-1 (Power analysis — MDC under H_W vs. global null significance)**

The reviewer correctly identified that the power analysis paragraph was misleading by stating that "the composite Mahalanobis test achieves significance" immediately after showing that all observed counts fall below their MDC under H_W.

We have revised the power paragraph (§5.5) to make the distinction explicit:

> "Note that the joint Mahalanobis test under H_W is also non-significant (D_M = 2.95, p = .121), so underpowerment applies at the multivariate level as well. The significant Mahalanobis test (D_M = 5.10, p < .001) reported in Section 4.4 refers to the *global* permutation null (H_G), not the within-woman null (H_W); both results are genuine, but they answer different scientific questions: H_G tests overrepresentation relative to all women in the sample, while H_W tests whether temporal ordering within each woman's own cycles adds overrepresentation beyond her individual cycle-length distribution."

This directly addresses the reviewer's concern that the text confounded global-null significance with within-woman-null non-significance. The interpretation is now internally consistent: the sample is underpowered under H_W both individually and jointly; the significant results are under H_G, which is a different (and scientifically less restrictive) null.

---

### Recommended Changes (Addressed)

**R-NEW-3 (Hotelling T² comparison)**

We have added a paragraph titled "Relationship to Hotelling's T²" in §5.2 (The Joint Multivariate Test):

> Under a Gaussian approximation, D_M² is asymptotically χ²(5), the same distributional assumption underlying T²; however, the exact permutation reference distribution used here replaces that parametric assumption with an empirical one, making the test distribution-free. Hotelling's T² requires the multivariate Gaussian assumption for validity; the permutation reference distribution does not. Moreover, the Bonferroni–Holm correction is conservative under positive correlation (as holds here), while the joint Mahalanobis test exploits the known positive correlation structure by pooling all five components simultaneously. For both reasons — distribution-free inference and proper treatment of dependence — the permutation-based Mahalanobis test is the preferred approach.

---

### Minor Issues (All Addressed)

**R-NEW-1 (MDC formula self-contained in Table 5)**

The Table 5 caption (tab:power) now includes: "MDC: minimum detectable count at 80% power (MDC = μ_W + 2.487σ_W, where z₀.₀₅ + z₀.₈₀ = 1.645 + 0.842 = 2.487)". The table is now self-contained for readers who encounter it independently of the text.

**R-NEW-2 (z-score discrepancy between Tables 2 and 5)**

A clarifying sentence has been added to the Results §4.3 paragraph following the Dilug-in-Dilug z = −2.07 result:

> "Note that the z values in Table 2 and Table 5 may differ slightly (e.g., Haflaga: −1.25 here vs. −1.24 in Table 5) because Table 2 uses the primary 50,000-run simulation (seed 17) while Table 5 uses an independent 50,000-run simulation for power computations; the differences are within Monte Carlo variability and do not affect any conclusion."

**R-NEW-4 (Week pattern anchor clarification)**

A new paragraph has been added to §5.3 (Non-Significant Patterns):

> "For the Week pattern, it is worth noting that the weekly anchor set W = {22, 29, 36, 43, 50} (in haflaga units) corresponds to cycle lengths {21, 28, 35, 42, 49}. Since cycles of 35 days or longer are well outside the normal range for the Fehring dataset (mean L = 29.3, SD = 3.9), nearly all Week patterns in this sample arise from haflagot of 22 (cycle length 21) or 29 (cycle length 28). The 28-day cycle is the mode of the dataset and the most biologically common length, so the Week pattern in this dataset effectively captures women with near-28-day cycles who also show the requisite two-consecutive-cycle consistency."

**R-NEW-5 (Relabelling invariance paragraph — anchor adjustment emphasis)**

The relabelling invariance paragraph (§5.6) has been strengthened. The phrase "with the anchor adjusted accordingly" is now italicised and emphasised. Additionally, a new sentence explicitly states the consequence of not adjusting:

> "This anchor adjustment is critical: computing on L_k with the unadjusted anchor (L = 30) yields only 12 Week-Dilug events instead of 18, which would be incorrect. Provided the anchor is shifted by the same −1 that maps H_k to L_k, the counts are unchanged."

**R-NEW-7 (Keywords — "within-woman permutation")**

The keyword list now reads: "menstrual cycle, vesset, Halacha, permutation test, **within-woman permutation**, multinomial sampling, natural family planning, Jewish law, polynomial hierarchy, dependence structure, multivariate randomization test"

The keyword "within-woman permutation" has been added and "multinomial sampling" retained (it remains relevant as one of the three null models).

---

### Statistical Concern S-2 (Minor)

**S-2 (p = .004 for Haflaga×DiD — claim of "robustness")**

We note that the text already reports that the permutation p-value (0.004) is below 0.01. We have not added additional text here as the existing caption statement — "The two methods agree closely, confirming that significance is not an artefact of asymptotic approximations" — directly addresses the point. If the reviewer believes a more explicit quantitative statement is warranted, we are happy to add: "the permutation p-value of 0.004 is below the conventional 0.01 threshold, providing strong evidence even under non-asymptotic conditions."

---

## Summary Table of Changes

| Issue | Type | Action | Location |
|-------|------|---------|----------|
| SMMR 2.6 (pseudocode) | Required | Added Algorithm 1 formal pseudocode | §3.5 |
| S-1 (power paragraph confusion) | Required | Rewrote to distinguish H_G vs H_W significance | §5.5 |
| R-NEW-3 (Hotelling T²) | Recommended | Added paragraph on T² vs D_M | §5.2 |
| R-NEW-1 (MDC self-contained) | Minor | Added formula to Table 5 caption | Tab.5 |
| R-NEW-2 (z discrepancy note) | Minor | Added Monte Carlo variability note | §4.3 |
| R-NEW-4 (Week anchor) | Moderate | Added 28-day clarification paragraph | §5.3 |
| R-NEW-5 (relabelling emphasis) | Moderate | Strengthened paragraph; added anchor error consequence | §5.6 |
| R-NEW-7 (keywords) | Editorial | Added "within-woman permutation" keyword | Keywords |

All changes are implemented in both `vesset_stat.tex` (main manuscript) and `vesset_stat_submission.tex` (blind submission version).
