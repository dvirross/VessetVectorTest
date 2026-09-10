# Round 2 Author Response
**Paper**: "Are Halachic Vesset Patterns Overrepresented in Empirical Menstrual Cycle Data? A Permutation Test"  
**Date**: 2026-09-10  
**Re**: Response to SMMR Reviewer 2 Round 2 Report

---

We thank the reviewer for the encouraging assessment and for the remaining minor comments. We are pleased that all required and recommended changes from Round 1 have been found satisfactory. We address each remaining issue below.

---

## Point-by-Point Response

**R2-1 (Required): Algorithm 1 — missing closing brace in line 8**

Corrected. The line `p_{D_M}^{(m)} ← (|{r : D_M^{(m,r)} ≥ D_M^{(m)}| + 1)/(B+1)` now correctly reads `(|\{r : D_M^{(m,r)} \geq D_M^{(m)}\}| + 1)/(B+1)` with the closing `\}` added.

**R2-2 (Recommended): p = .004 sentence**

Added to Table 3 caption in `vesset_stat.tex` (main file) and to the results paragraph in `vesset_stat_submission.tex`: "Notably, the Haflaga × Dilug-in-Dilug permutation p-value (0.004) falls below 0.01, providing strong evidence against the null that co-occurrence arises by chance even under non-asymptotic conditions."

**R2-3 (Recommended): Gaussian approximation caveat in §5.2**

The Hotelling T² paragraph now reads: "Under a Gaussian approximation (i.e., assuming the count vector is approximately multivariate normal, as holds asymptotically for large counts), D_M² is asymptotically χ²(5). However, with observed counts as small as 18 (Week-Dilug), this approximation may not hold reliably — motivating the permutation reference distribution, which replaces the parametric assumption with an empirical one and makes the test fully distribution-free."

**R2-4 (Optional): §5.6 phrasing**

Revised to: "...yields only 12 Week-Dilug events instead of the correct 18." (more precise than "which would be incorrect").

**R-NEW-6 (Verification): Table 3 column widths**

The main `vesset_stat.tex` Table 3 has been reviewed. The column widths are: L{4.6cm} C{1.2cm} C{1.2cm} C{1.2cm} C{0.8cm} C{1.2cm} C{1.2cm}, total ≈ 4.6 + 6×1.2 − approx. adjustments = 11.8cm for content, well within SMMR page width of ~15.5cm (US letter − 1-inch margins × 2). We confirm no overflow is expected. (Note: the submission file uses a 5-column version of Table 3 without p_perm, further reducing width concerns there.)

---

## Summary of Changes

| Issue | Action |
|-------|--------|
| R2-1 (Algorithm notation) | Fixed closing brace in Algorithm 1 |
| R2-2 (p=.004 sentence) | Added to Table 3 caption / text |
| R2-3 (Gaussian approximation caveat) | Added to §5.2 Hotelling T² paragraph |
| R2-4 (§5.6 phrasing) | Minor phrasing polish |
| R-NEW-6 (Table overflow) | Verified: no overflow expected |
