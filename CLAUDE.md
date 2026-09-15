# CLAUDE.md — VessetVectorTest

Statistical analysis of halachic vesset patterns in menstrual cycle data.  
Paper targeting: **Biometrical Journal (Wiley, Q1)** — shifted from SMMR 2026-09-15.

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

## 3. Biometrical Journal requirements (target journal)

**Submission portal**: https://authors.wiley.com/journal/BIMJ  
**Article type**: Research Article (4500 words recommended) or Case Study (4500 words)  
**Peer review**: Single-blind — reviewers anonymous, authors not. **Submission does NOT need to be anonymous.** The anonymous submission file (`vesset_stat_submission.tex`) is no longer needed for this journal.

### Mandatory formatting constraints

| Requirement | Current status | Action needed |
|---|---|---|
| **No footnotes** | ❌ Paper has footnotes | Convert to parenthetical text or appendix |
| **≤5 keywords, alphabetical** | ❌ Currently 11 keywords | Reduce and sort |
| **Unstructured abstract** | ✅ Already unstructured | No change |
| **Data Availability Statement** | ❌ Missing | Add to manuscript |
| **Conflict of Interest statement** | ❌ Missing | Add (state none) |
| **IRB/ethics statement** | ❌ Missing | Add (secondary data; IRB per Fehring 2013) |
| **APA reference style** | ❌ Currently plainnat | Change `\bibliographystyle` to `apalike` |
| **ORCID** | ❌ Not in manuscript | Add to title page |
| **Reproducible Research ZIP** | ✅ Code on GitHub; data public | Prepare code+data ZIP at revision |
| **Figures as separate files** | ✅ Already PDF/PNG | Required at revision stage |

### Scope fit

BJ requires methodological development motivated by a real problem — our paper fits as a **Case Study** (novel application of permutation-based Mahalanobis joint test to sequentially-defined count vectors) or Research Article (the H_W within-woman permutation null + polynomial hierarchy are methodological contributions). Framing should emphasise: (1) the permutation-based Mahalanobis test as a distribution-free method for correlated count vectors, (2) the within-woman permutation null as a design-adapted exchangeability block, (3) the real-world menstrual cycle application as motivation.

### Cover letter required statements (BJ-specific)
- All authors concur with the submission
- All funding listed in Acknowledgements (include any funding or state none)
- The manuscript has not been submitted elsewhere
- Ethics: secondary analysis of publicly archived dataset (Fehring et al. 2013, Marquette University); IRB approval obtained for original data collection; no direct patient contact in this study
- Conflicts of interest: none (or disclose)
- Closely related work (PhD thesis, prior publications): reference in cover letter; submit as "Supplementary material for review only"

### Things NOT required for BJ (vs SMMR)
- No anonymous submission file needed
- No 250-word abstract limit (BJ has no stated limit)
- No colour-free figure requirement (BJ does not mandate B&W figures)
- Software names in Courier New typeface (NumPy → `NumPy`)
- Sentences must not start with symbols or numbers

---

## 4. Two-file LaTeX strategy

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

## 5. Notebook: keep in sync with the paper

The single source of truth for all analyses and figures is `notebooks/analysis.ipynb`.
**After every paper change that affects an analysis or figure, update the notebook too.**

### Figure generators

| Paper figure | Source |
|---|---|
| `fig1_cycle_distribution` | notebook cell (cycle length histogram + per-woman scatter) |
| `fig2_consort_flow` | `paper/figures/gen_consort_flow.py` — run directly |
| `fig3_null_distributions` | notebook cell (global perm + multinomial PMF bars) |
| `fig4_zscores` | notebook cell (z-score grouped bar chart) |
| `fig5_cdfs` | notebook cell (empirical CDFs) |
| `fig6_stratified_comparison` | notebook cell (H_G / H_iid / H_W KDE overlay) |
| `fig7_dm_nulldist` | notebook cell (3-panel Mahalanobis D_M null histograms) |
| `fig8_chisq_heatmap` | `paper/figures/regen_heatmap.py` — run directly (§ markers, hardcoded values) |
| `fig9_joint_3panels` | notebook cell (3-panel 3D scatter, focused 3-pattern joint test) |
| `figA1_heaping_diagnostic` | notebook cell (haflaga bar chart, heaping sensitivity check) |

### Notebook invariants

- All `plt.savefig` calls must use the **current** figure filename (fig1…figA1); never old names.
- New analyses added to the paper must have a corresponding notebook cell.
- Stale cells (superseded analyses/figures) should be clearly marked `# STALE` or removed.

---

## 6. LaTeX: always compile and check


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

## 7. Terminology — get these right

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

## 8. Primary vs secondary results

**PRIMARY — Within-woman null H_W (all n = 118 women, seed 17, B = 50,000):**

| Pattern | Obs | μ_W | σ_W | z | p |
|---|---|---|---|---|---|
| Haflaga | 26 | 31.85 | 4.68 | −1.25 | .249 |
| Dilug | 61 | 67.07 | 7.36 | −0.83 | .455 |
| Week | 32 | 33.81 | 3.88 | −0.47 | .741 |
| Week-Dilug | 18 | 20.94 | 3.22 | −0.91 | .454 |
| Dilug-in-Dilug | 27 | 39.14 | 5.86 | −2.07 | .040 (Holm adj .198) |
| **Joint** | — | — | — | D_M = 2.95 | p = .121 |

All patterns non-significant. All observed counts below null means.  
Abstract and conclusion **must foreground H_W** as the central result.

**SECONDARY — Global null H_G (permutation, seed 17, B = 50,000):**
- Haflaga: z = 4.511, p < .001 (Holm adj p = .0003) ✓ significant
- Dilug: z = 1.830, p = .087 — **exploratory only** (Holm adj p = .346, does not survive; two-sided)
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

## 9. MDC / Power section

The MDC values are **approximate normal-theory sensitivity benchmarks**, not exact power guarantees.  
Formula: MDC = μ_W + 2.487σ_W (normal-theory, 97.5th percentile approximation).  
Label as: "approximate sensitivity analysis under a normal-theory location-shift approximation."  
Do not claim exact power without simulation-based curves. Future work caveat required.

---

## 10. Self-citations — keep them low

- `ross2022phd` / `anon2022phd` → ≤ 5 total occurrences
- `ross2021statistical` / `anon2021statistical` → ≤ 2 total occurrences
- When tempted to cite the PhD for methodology, describe the method inline instead

---

## 11. Figures

Biometrical Journal does **not** require colour-free figures (unlike SMMR).  
Figures submitted as separate PDF/EPS/TIFF files at revision stage.  
Figure legends must appear both beneath each image and as a complete list in the manuscript text.  
Do not create figures using LaTeX code — use Python (matplotlib) and export as PDF.

---

## 12. Key file paths

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

## 13. Current paper status

The paper has been revised in response to three rounds of AI-generated peer reviews (SMMR). All critiques resolved.

### Round 1–2 revisions (earlier sessions)
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
- ✅ §4.1 absorbed into §4.2 (section merged)
- ✅ §5.1+§5.2 merged into single "Pattern-Specific Findings" section
- ✅ §5.3 trimmed from ~1,000 → ~300 words
- ✅ Limitations condensed from 9 paragraphs to 4
- ✅ R1 #4: Added condition-number caveat in §3.4 (Mahalanobis power note)
- ✅ R1 #8: "selects" → "is associated with"; added cycle-count precision caveat in SD section
- ✅ R1 minor (tied values): permutation tie-handling note added
- ✅ R2 3.3: Marginal results framed as complementary summaries
- ✅ R1 minor (cycles ≥ 35 days): replaced "well outside normal range" with data-supported frequency (9.6%)
- ✅ R3 #3: Conclusion novelty claim softened; Mahalanobis test framed as adapted approach
- ✅ CONSORT figure: fixed ≤ 60 → ≤ 54 days; SD 8.7 → SD 6.1; range 1–45 → range 5–45
- ✅ fig:chisq caption: updated to reflect only Haflaga×DiD surviving Holm (adj p = .040)

### Round 3 revisions (commit d06d161, 2026-09-15)
- ✅ R1-3.3: Joint Mahalanobis test explicitly labelled **primary**; individual Holm-corrected marginals labelled **secondary descriptive summaries**
- ✅ R1-3.3: Acknowledged similar quadratic-form tests (energy-distance); reframed as adaptation, not first-ever
- ✅ R1-Minor: Onset-bias × H_W interaction: H_W conditions on each woman's cycle-length multiset, so consistent within-woman onset shifts are preserved identically in permuted sequences → H_W unaffected by that misclassification form
- ✅ R2-2: "addresses cycle-count confound" → "mitigates…reduces but does not eliminate"
- ✅ R3-2: Added underpowered caveat to abstract (all observed counts below MDC benchmark)
- ✅ R3-2: Added non-detection sentence in conclusion: underpowered; not proof of absence
- ✅ R3-5: Added confirmatory/exploratory/robustness paragraph at start of Discussion
- ✅ R3-Minor: Added Week z-score difference note in fig4 caption
- ✅ CLAUDE.md §9: Removed all specific colour names from all figure captions (blue/green/red/orange → solid/dashed/dotted/light fill/hatched)
- ✅ Cover letter written and pushed (commit a2c1644)
- ✅ Submission file regenerated and re-verified (0 non-anon occurrences)
- ✅ Both PDFs compiled (36 pages each) and in repo
- ✅ Both Overleaf zips rebuilt

**Note:** R1 Minor figure legend/title overlap was raised against the old two-panel 3D scatter figures (fig6_5d_permutation, fig7_5d_multinomial). Those were replaced by the combined 3-panel figure (fig9_joint_3panels) and the D_M null distribution figure (fig7_dm_nulldist), which were redesigned to fix those layout issues. The reviewer never saw the new figures — this concern is resolved.

### Round 4 revisions (2026-09-15) — responding to next reviewer batch
- ✅ **CRITICAL MATH FIX**: Set inclusion direction corrected in 4 locations (lines 158, 476, 478, 1172): "Week-Dilug ⊆ Haflaga-30" → "Haflaga-30 ⊆ Week-Dilug". Methods section (lines 361, 751-752) was already correct.
- ✅ CONSORT figure caption renamed: "CONSORT-style" → "Participant and data flow diagram"
- ✅ fig3 null-distributions caption: KDE described as "smoothed density estimates for visual comparison only; primary inference from empirical PMF"
- ✅ Hotelling's T² paragraph rewritten: asymptotic/parametric framing removed; Holm vs Mahalanobis correctly framed as answering different questions (FWER vs global null), not competing
- ✅ Conclusion "avoids conservatism" → "Holm and Mahalanobis are complementary, not competing"
- ✅ Dataset section: added 5-cycle minimum explanation, parity-not-recorded note, two-cohort disclosure
- ✅ Heaping diagnostic rewritten as sensitivity diagnostic; added "formal test not performed" caveat
- ✅ Global null Haflaga result tempered: "biologically interpretable" → "consistent with expectation"; generalisability caveat added
- ✅ Pairwise section: added deterministic vs stochastic distinction paragraph (set-containment is arithmetic, not stochastic evidence)
- ✅ fig:chisq caption: "consistent with polynomial-hierarchy" → "compatible with proximity argument (all pairwise tests exploratory)"
- ✅ Exact exceedance counts added for joint tests (e.g., "1 of 50,000 null replicates exceeded")
- ✅ Σ̂ formula added explicitly: sample covariance of null replicates, not observations
- ✅ "First formal test" qualifier added to Introduction and Conclusion (with "no systematic literature search" caveat)
- ✅ Broken cross-reference labels added: \label{sec:methods}, \label{sec:limitations}
- ✅ Submission file anonymised (0 non-anon occurrences)
- ✅ Both PDFs recompiled (37 pages each)
- ✅ Both Overleaf zips rebuilt

### Round 4 continued (same session)
- ✅ R2-2.2: Explicit H_0/H_1 added for 3-pattern primary test; structural basis for inclusion (containment + polynomial proximity) and exclusion (Dilug hypothesis-generating; Week no structural link) stated formally
- ✅ R2-2.5: Rate-based Spearman analysis already in paper; code added to notebook (cell 33, before chi-square section); seed 17, B=50,000 permutations
- ✅ R3 Major 4 (derivations): sec:mathdep prose shortened; Discussion "Joint Multivariate Test" subsection condensed (~40% shorter); Taylor expansion analogy compressed
- ✅ R2 Minor 3.1: Min run-length examples added to Methods (before Algorithm 1): lengths 2, 3, 4 with concrete sequences
- ✅ Notebook: rate-based Spearman analysis code added as cell 33
- ✅ Both PDFs recompiled (38 pages each, no errors)
- ✅ Submission re-anonymised (0 non-anon occurrences)
- ✅ Both Overleaf zips rebuilt

### Round 5 revisions (2026-09-15) — two-sided switch + audit
- ✅ Abstract trimmed to ~217 words (was ~283; SMMR limit 250)
- ✅ Hotelling's T² paragraph added: normality assumption violated by small integer counts; parametric form inappropriate
- ✅ max-T paragraph added: less sensitive to multi-pattern displacement; same conclusion in our data
- ✅ **All marginal tests switched to two-sided p-values** (scientific rationale: patterns confirm regularity in either direction; DiD underrepresentation is informative)
- ✅ Two-sided p-value formula added to Methods and Algorithm pseudocode
- ✅ H_W table p-values updated: .916/.813/.722/.857/.988 → .249/.455/.741/.454/.040
- ✅ H_G Dilug: p=.043 (bold) → p=.087 (not significant); Holm adj .173 → .346; throughout text
- ✅ DiD H_W note rewritten: underrepresentation consistent with DiD detecting genuine second-order cycle regularity (was: one-sided dismissal)
- ✅ "all p ≥ .722" → "all Holm-corrected p ≥ .198" in abstract, results, conclusion
- ✅ MDC formula: z_{0.025}+z_{0.80}=2.802 (was z_{0.05}+z_{0.80}=2.487, one-sided)
- ✅ MDC table values updated: Haflaga/Week 44→45, Dilug 86→88, WD 29→30, DiD 54→56
- ✅ Figure cross-reference audit completed: fig:consort had no in-text \ref — added pointer; H_W table caption p-value formula updated to two-sided
- ✅ \dag (†, cross shape) replaced with $^*$ (asterisk) in pairwise table and caption
- ✅ Both PDFs recompiled, anonymisation verified (0 non-anon occurrences), both Overleaf zips rebuilt

### Round 6 additions (2026-09-15) — notebook audit + sync
- ✅ Notebook audit: identified all 10 paper figures and their generators
- ✅ Fixed all savefig filenames in notebook to match current paper figure names (fig0→fig1, fig_stratified→fig6, fig1_null→fig3, fig2_zscores→fig4, fig4_cdfs→fig5, fig5_chisq→fig8)
- ✅ Cell 29 (stale Euclidean joint figure) replaced with fig7_dm_nulldist (Mahalanobis D_M, 3-panel)
- ✅ Cell 38 (stale old 3D scatters fig6/fig7) replaced with fig9_joint_3panels generator
- ✅ Added figA1_heaping_diagnostic notebook cell (haflaga distribution bar chart)
- ✅ Added note cell listing standalone scripts (fig2_consort_flow, fig8_chisq_heatmap)
- ✅ regen_heatmap.py added to paper/figures/ (permanent location, previously only in scratchpad)
- ✅ CLAUDE.md §4 added: Notebook requirements and figure generator table

### Round 7 revisions (2026-09-15) — strategic reframe + citation correction
- ✅ **Primary/secondary hierarchy clarified throughout**: 5-pattern Mahalanobis is now unambiguously primary; 3-pattern analysis relabelled "exploratory focused subspace analysis" serving as geometric complement
- ✅ Abstract: Updated to foreground 5-pattern primary (D_M=5.10, p<.001); removed "underpowered" → "limited sensitivity" framing
- ✅ Methods §3.4: Rewritten so primary joint test = 5-pattern; 3-pattern = exploratory with rationale for exclusion of Dilug and Week
- ✅ Results "Joint Multivariate Test": Replaced 3-pattern primary with 5-pattern primary result block; 3-pattern demoted to exploratory
- ✅ "Confirmatory vs exploratory" paragraph: Fixed item (iii) to 5-pattern D_M=5.10; 3-pattern moved to exploratory list
- ✅ Conclusion: Updated to 5-pattern primary; 3-pattern described as "exploratory focused analysis"
- ✅ Figure 9 suptitle + caption: Reframed as exploratory geometric complement to primary 5-pattern analysis
- ✅ Fig 4 caption: Two-sided z-thresholds (|z|=1.96, |z|=2.576); Dilug no longer described as exceeding threshold
- ✅ Fig 6: Regenerated with per-panel two-sided H_W p-value annotations
- ✅ New intro paragraph: Correctly names Ecochard et al. 2024 (not Cortet); notes Fehring/Schneider overlap; discloses probable dataset overlap (Marquette 1,649/159 vs present 1,554/118); explains H_W robustness to autocorrelation
- ✅ New intro paragraph: Added Winkler 2014 and Nichols 2002 citations for exchangeability-block permutation precedent
- ✅ Novelty claims softened in 3 locations: "first formal hypothesis test" → "applies a formal permutation-based hypothesis test… to our knowledge the first such test applied to this specific domain"
- ✅ MDC section: Relabelled as "approximate sensitivity benchmarks"; removed power claim language
- ✅ references.bib: Added ecochard2024, winkler2014, nichols2002 (replacing incorrect cortet2024)
- ✅ Notebook cells patched: fig4 two-sided thresholds; fig6 H_W p-value annotations; fig9 exploratory suptitle
- ✅ Both PDFs recompiled clean (0 errors)
- ✅ Submission re-anonymised (0 non-anon occurrences)
- ✅ Both Overleaf zips rebuilt

**Outstanding (require author decisions):**
- [ ] Journal targeting: Submit to application journal (Statistics in Medicine, JRSS-C, BMC MRM) rather than SMMR; cover letter revision needed
- [ ] R1 Major 1: New simulation study — power comparison vs Hotelling's T², max-T, Bonferroni-Holm (text response added; simulation deferred)
- [ ] R1 Major 2 + R2-2.8: Simulation-based power analysis injecting synthetic vesset events at known prevalence rates (MDC table is normal-theory approximation; deferred to future work)
- [ ] R3 Major 4 (figure): Move one figure to supplement (deferred)
