# CLAUDE.md — VessetVectorTest

Statistical analysis of halachic vesset patterns in menstrual cycle data.  
Paper targeting: **Statistical Methods in Medical Research (SAGE, Q1)**.

> **Before every session:** ask the user what they want to do. Do not start a series of tasks without explicit instruction. Confirm intent before any chain of edits.

---

## 1. Plan before touching anything

Before writing code or editing LaTeX:
- State what you intend to change and why
- If it affects simulation results (N = 50,000 runs per null), flag it — rerunning takes ~30 min
- If it affects the paper narrative (primary vs secondary result), confirm with user first

**Pre-computed results live in:**
- `strat_full_118.pkl` — H_W within-woman permutation null (seed 17, B = 50,000, all 118 women)
- `sim_results.npz` — global permutation + multinomial nulls (seed 17, B = 50,000)

Do **not** rerun simulations unless the counting logic changes.  
`boot_arr` inside `sim_results.npz` is an **older within-woman bootstrap** (different null, Haflaga mean ≈ 11.45). It is NOT H_W. Always load `strat_full_118.pkl` for H_W values.

---

## 2. Git workflow

- All work committed to **`main`** branch (previous sessions established this)
- Commit and push after every few changes — token limits can truncate sessions
- Commit message format: brief imperative summary + attribution footer:
  ```
  Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01L8iL1JWncL7yY1QoJSEnS6
  ```

---

## 3. Two-file LaTeX strategy

| File | Purpose |
|---|---|
| `paper/vesset_stat.tex` | Main file — author name, `ross*` citations, full URLs |
| `paper/vesset_stat_submission.tex` | Anonymous submission — `anon*` citations, no names/URLs, `unsrtnat` bib style |

**After every edit to `vesset_stat.tex`**, sync to submission by:
1. Copying and applying these substitutions:
   - `ross2022phd` → `anon2022phd`, `ross2022markov` → `anon2022markov`, etc.
   - `ross2019reliability` → `anon2019reliability`, `ross2021statistical` → `anon2021statistical`
   - `ross_planned` → `anon_planned`
   - Author block → `[Author details removed for peer review]`
   - Code URLs → `[URL removed for peer review]`
   - `\bibliographystyle{plainnat}` → `\bibliographystyle{unsrtnat}`
2. Check: `grep -c "ross20\|Dvir\|dvirross" vesset_stat_submission.tex` → must return 0

**Submission zips:**
- `vesset_overleaf_main.zip` — main version for Overleaf
- `vesset_overleaf_submission.zip` — anonymous version for journal upload

Rebuild zips whenever tex or figures change.

---

## 4. LaTeX: always compile and check

After every LaTeX change:
```bash
cd /tmp/build_dir
pdflatex -interaction=nonstopmode main.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex   # third pass for cross-refs
grep "^!" main.log                           # must be empty
```

**Known recurring LaTeX bugs in this project:**
- `\nat` appearing from line-break `and\nat` → fix to `and at`
- `\bibliographystyle{plainnat}` in submission → must be `unsrtnat`
- `\multicolumn{}` unclosed braces
- Stray `\n` from Python string generation appearing as literal `\n`

Non-fatal warnings to ignore: `titlesec` "entered in horizontal mode", `rerunfilecheck`.

---

## 5. Terminology — get these right

| Wrong | Correct |
|---|---|
| "her haflaga sequence" (H_W section) | "her cycle-length sequence" |
| stratified permutation (H_W) | within-woman permutation null |
| bootstrap (for H_W) | permutation (H_W is without-replacement, not bootstrap) |
| global permutation assumes independence | global null tests population-level exchangeability |
| Bonferroni assumes independence | Bonferroni valid under any dependence (conservative under positive correlation) |
| establish / prove / demonstrate | suggest / provide evidence that |
| essentially unambiguous | (remove entirely) |
| first rigorous validation | to our knowledge, the first formal statistical test |
| 100× more significant | 100-fold smaller p-value |
| conservative (permutation vs bootstrap) | the two correspond to different scientific nulls, not comparable as more/less conservative |
| Dilug "significant" | Dilug "nominally significant; does not survive Holm correction; exploratory" |
| weighted sum of squared z-scores | D_M² = (O−μ)ᵀ Σ⁻¹ (O−μ); approximately but not exactly a weighted sum |
| collective overrepresentation | significant displacement from null distribution |

`np.random.permutation` = uniform random permutation **without** replacement.  
H_W is a valid conditional Monte Carlo test. Not a bootstrap.

---

## 6. Primary vs secondary results

**PRIMARY — Within-woman null H_W (all n = 118 women, seed 17, B = 50,000):**

| Pattern | Obs | μ_W | σ_W | z | p |
|---|---|---|---|---|---|
| Haflaga | 26 | 31.85 | 4.68 | −1.25 | .916 |
| Dilug | 61 | 67.07 | 7.36 | −0.83 | .813 |
| Week | 32 | 33.81 | 3.88 | −0.47 | .722 |
| Week-Dilug | 18 | 20.94 | 3.22 | −0.91 | .857 |
| Dilug-in-Dilug | 27 | 39.14 | 5.86 | −2.07 | .988 |
| **Joint** | — | — | — | D_M = 2.95 | p = .121 |

All patterns non-significant. All observed counts below null means.  
Abstract and conclusion **must foreground H_W** as the central result.

**SECONDARY — Global null H_G (permutation, seed 17, B = 50,000):**
- Haflaga: z = 4.511, p < .001 (Holm adj p = .0003) ✓ significant
- Dilug: z = 1.830, p = .043 — **exploratory only** (Holm adj p = .173, does not survive)
- Week: z = 1.238, p = .108
- Week-Dilug: z = 0.529, p = .305
- Dilug-in-Dilug: z = −0.033, p = .521
- Joint: D_M = 5.10 (perm), D_M = 4.84 (multi), both p < .001

**Never call Dilug "significant"** without "exploratory" and Holm correction note.

**SD comparison (Haflaga vs no-vesset women):**
- Haflaga n = 20: mean within-person SD = 1.54 days (SD_between = 0.59)
- No vesset n = 35: mean within-person SD = 3.17 days (SD_between = 1.75)
- Mann–Whitney U = 133, p < .001; Welch t = −5.05, df = 45.5, p < .001
- 95% CI for mean difference: [−2.29, −0.98] days; Cohen's d = −1.13

---

## 7. MDC / Power section

The MDC values are **approximate normal-theory sensitivity benchmarks**, not exact power guarantees.  
Formula: MDC = μ_W + 2.487σ_W (normal-theory, 97.5th percentile approximation).  
Label as: "approximate sensitivity analysis under a normal-theory location-shift approximation."  
Do not claim exact power without simulation-based curves. Future work caveat required.

---

## 8. Self-citations — keep them low

- `ross2022phd` / `anon2022phd` → ≤ 5 total occurrences
- `ross2021statistical` / `anon2021statistical` → ≤ 2 total occurrences
- When tempted to cite the PhD for methodology, describe the method inline instead

---

## 9. Figures — colour-free captions

SMMR requires figures comprehensible in black-and-white.  
Captions must **not** mention: blue, red, green, orange, dashed red, solid blue, etc.  
Use instead: solid, dashed, dotted, filled circle, open circle, arrow — not colour names.

---

## 10. Key file paths

```
data/FilteredData.csv          1,554 cycles, 118 women (filtered)
sim_results.npz                global perm + multinomial (seed 17, B = 50,000)
strat_full_118.pkl             H_W results (seed 17, B = 50,000, all 118 women)
paper/vesset_stat.tex          main LaTeX (with author name)
paper/vesset_stat_submission.tex  anonymous submission LaTeX
paper/vesset_stat.pdf          compiled main PDF
paper/vesset_stat_submission.pdf  compiled anonymous PDF
paper/references.bib           bibliography (ross* and anon* entries)
paper/figures/                 all figures as PDF + PNG
vesset_overleaf_main.zip       Overleaf zip (main version)
vesset_overleaf_submission.zip Overleaf zip (anonymous version)
notebooks/analysis.ipynb       main analysis notebook
```

PhD repo (raw data context): `github.com/dvirross/PhD`  
Analysis repo: `github.com/dvirross/VessetVectorTest`

---

## 11. Current paper status

The paper has been revised in response to three AI-generated peer reviews. Key changes applied:

- ✅ Abstract hedged ("establish" → "suggest"; removed "exploiting serial temporal dependence")
- ✅ Discussion restructured: H_W as primary, global null as secondary/descriptive
- ✅ Dilug treated as hypothesis-generating throughout (Holm adj p = .173 noted)
- ✅ SD comparison formalised with Mann–Whitney, Welch t, 95% CI, Cohen's d
- ✅ Pairwise dependence tests labelled "exploratory" throughout
- ✅ Mahalanobis formula corrected: D_M² = (O−μ)ᵀ Σ⁻¹ (O−μ) explicit
- ✅ MDC relabelled as approximate normal-theory sensitivity analysis
- ✅ Clinical section shortened substantially
- ✅ "Conservative" language removed; permutation vs bootstrap framing corrected
- ✅ Exchangeability-narrower-than-independence explained
- ✅ Simulation algorithm pseudocode added
- ✅ Submission file synced and anonymised
- ✅ Both PDFs compiled and in repo
- ✅ Both Overleaf zips rebuilt

**Outstanding:**
- [ ] Figure cross-reference audit: verify all "Figure X" in text match actual captions
- [ ] Consider a cover letter for journal submission
- [ ] Respond to next round of AI reviews once received
