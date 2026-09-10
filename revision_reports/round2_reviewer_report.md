# Round 2 Reviewer Report — SMMR Reviewer Simulation
**Paper**: "Are Halachic Vesset Patterns Overrepresented in Empirical Menstrual Cycle Data? A Permutation Test"  
**Date**: 2026-09-10  
**Status**: Second revision review (acting as SMMR Reviewer 2)

---

## Summary Decision

**Minor Revision** — The second revision has comprehensively addressed all required and recommended issues from Round 1. The addition of Algorithm 1 (simulation pseudocode), the corrected power analysis paragraph distinguishing H_G from H_W significance, the Hotelling T² comparison, and the minor textual clarifications are all satisfactory. The manuscript is substantially improved and the statistical reasoning is now clearly articulated throughout.

A small number of minor issues remain, none of which affect the validity of the analysis or the conclusions.

---

## Point-by-Point Assessment

### Previously Raised Concerns: Status

**S-1 (Power analysis paragraph)**: ✅ RESOLVED  
The revised paragraph now explicitly states: "the joint Mahalanobis test under H_W is also non-significant (D_M = 2.95, p = .121), so underpowerment applies at the multivariate level as well. The significant Mahalanobis test (D_M = 5.10, p < .001) refers to the global permutation null (H_G), not the within-woman null (H_W)." This is the correct and scientifically appropriate framing. The distinction between the two nulls is now consistently maintained throughout §5.5.

**SMMR 2.6 (Algorithm pseudocode)**: ✅ RESOLVED  
Algorithm 1 is clear and covers all required steps: observation, simulation loop with per-null-model specification, covariance estimation, and p-value computation. The $(b+1)/(B+1)$ convention is explicitly stated in the algorithm. The H_W step correctly specifies "per-woman Fisher–Yates shuffle without replacement," matching the implementation.

**R-NEW-3 (Hotelling T²)**: ✅ RESOLVED  
The new paragraph in §5.2 correctly distinguishes the permutation-based D_M² from Hotelling's T² on two axes: (a) distributional assumption (parametric Gaussian vs. empirical permutation reference) and (b) treatment of inter-component correlation (Bonferroni conservative vs. joint pooling). The statement "D_M² is asymptotically χ²(5)" is correct under Gaussian approximation; this connection to the χ² distribution is helpful for readers familiar with classical multivariate testing.

**R-NEW-1 (MDC formula self-contained)**: ✅ RESOLVED  
The Table 5 caption now reads "MDC = μ_W + 2.487σ_W, where z₀.₀₅ + z₀.₈₀ = 1.645 + 0.842 = 2.487." This is correct and makes the table fully self-contained.

**R-NEW-2 (z-score discrepancy)**: ✅ RESOLVED  
The added note in §4.3 correctly attributes the small differences (e.g., −1.25 vs. −1.24 for Haflaga) to Monte Carlo variability across independent simulation runs. This prevents reader confusion.

**R-NEW-4 (Week anchor clarification)**: ✅ RESOLVED  
The new paragraph in §5.3 correctly identifies that the weekly anchor set {22, 29, 36, 43, 50} in haflaga units corresponds to cycle lengths {21, 28, 35, 42, 49}, and that cycles ≥35 days are uncommon in this sample. The statement that "the Week pattern in this dataset effectively captures women with near-28-day cycles" is factually correct and addresses the biological interpretation question.

**R-NEW-5 (Relabelling invariance)**: ✅ RESOLVED  
The paragraph now uses italics to emphasise "with the anchor adjusted accordingly" and explicitly warns: "computing on L_k with the unadjusted anchor (L = 30) yields only 12 Week-Dilug events instead of 18, which would be incorrect." This is the correct and important warning, preventing replication errors.

**R-NEW-7 (Keywords)**: ✅ RESOLVED  
"within-woman permutation" has been added to the keyword list.

---

### Remaining Issues

**R2-1 (Minor/Editorial): Algorithm 1 — notation inconsistency in line 8**

In Algorithm 1 (main tex), line 8 reads:
```
p_{D_M}^{(m)} ← (|{r : D_M^{(m,r)} ≥ D_M^{(m)}| + 1) / (B+1)
```
The curly brace in `|{r : D_M^{(m,r)} ≥ D_M^{(m)}|` appears to be missing a closing brace: it should be `|\{r : D_M^{(m,r)} \geq D_M^{(m)}\}|`. This is likely a minor LaTeX rendering issue in the pseudocode; please verify the compiled output renders correctly and add the missing closing brace if needed.

**R2-2 (Minor): S-2 — permutation p-value sentence (recommended addition)**

In the Round 1 author response, the authors acknowledged S-2 but declined to add an explicit sentence about the 0.004 permutation p-value being below 0.01. On reflection, a single sentence in the Table 3 caption or in §4.5 text — "The Haflaga×Dilug-in-Dilug permutation p-value (0.004) falls below 0.01, providing strong evidence against the null that co-occurrence arises by chance even under non-asymptotic conditions" — would strengthen the paper modestly and is recommended (though not required for acceptance).

**R2-3 (Minor): §5.2 Hotelling T² paragraph — clarify "under Gaussian approximation"**

The new paragraph states "under a Gaussian approximation, D_M² is asymptotically χ²(5)." This is correct, but it would be more precise to note that the approximation requires the count vector to be approximately multivariate normally distributed, which holds when pattern counts are large (central limit theorem). With observed counts as small as 18 (Week-Dilug), this approximation may not be fully reliable — which is one of the reasons why the permutation reference distribution is preferred. Adding a parenthetical "(an approximation that may not hold well with small observed counts, motivating the permutation approach)" would strengthen the methodological justification for the permutation-based test.

**R2-4 (Very minor): §5.6 relabelling invariance — sentence structure is slightly awkward**

The phrase "This anchor adjustment is critical: computing on L_k with the unadjusted anchor (L = 30) yields only 12 Week-Dilug events instead of 18, which would be incorrect" reads well, but "which would be incorrect" is slightly ambiguous (incorrect relative to what?). Consider: "...yields only 12 Week-Dilug events instead of the correct 18." This is purely cosmetic.

---

### Editorial Note

**R-NEW-6 (Table 3 column widths)**: The reviewer notes that the previous concern about Table 3 column widths was raised in the Round 1 report but is not explicitly addressed in the author response. The authors are asked to confirm that the table compiles without overflow warnings in the SMMR-formatted PDF. (This was listed as a minor issue requiring only verification, not text changes, and may already be fine.)

---

## Recommendation

**Minor Revision** with the following hierarchy:

1. (**Required**) Fix the notation in Algorithm 1 line 8 — missing closing brace (R2-1)
2. (**Recommended**) Add the S-2 sentence about p = .004 (R2-2)
3. (**Recommended**) Clarify the Gaussian approximation caveat in §5.2 (R2-3)
4. (**Optional**) Polish §5.6 phrasing (R2-4)
5. (**Verification only**) Confirm Table 3 does not overflow in compiled PDF (R-NEW-6)

The fundamental statistical analysis, methodological contributions, and scientific conclusions are all well-supported. The manuscript is approaching publication quality. After the minor revision, the paper should be suitable for acceptance.
