# External replication on the Utah Creighton Model cohort — results (2026-09-17)

> **Superseded numbers (2026-09-22).** After this report was written the Week-Dilug rule was generalised from H = 30 to all one-weekday-progression anchors H ≡ 2 (mod 7) = {23, 30, 37, …} and the whole pipeline was rerun. Utah Week-Dilug count 40 → 42 (2 events at H = 37); joint D_M 9.99/9.53/1.20 → 10.04/9.56/1.22 (H_W p .921 → .916). Current numbers: `external_replication/results/comparison.md` and `summary.json`; the manuscript and SI Section S4 are up to date.


All numbers from `external_replication/results/` (B = 50,000, seed 17; same engine, seeds and
constants as the Fehring analysis). Preprocessing and deviations:
`external_replication/reports/utah_preprocessing_and_deviations.md`. Reproduce with

```bash
python external_replication/scripts/prepare_utah.py           # raw deposit -> external_replication/data/sequences.csv
python external_replication/scripts/run_external.py --nproc 3   # ~30 min on 3 cores -> external_replication/results/
```

## 1. Samples

| | Fehring (Marquette) | Utah (Creighton Model) |
|---|---|---|
| women / cycles | 118 / 1,554 | 270 / 2,139 |
| cycles per woman: mean (SD), median, range | 13.2 (6.1), 12, 5–45 | 7.9 (2.7), 7, 5–15 |
| cycle length L: mean (SD), median, range | 29.3 (3.9), 28.5, 18–54 | 30.0 (6.2), 29, 17–98 |
| share L ≥ 35 | 9.6 % | 12.0 % |
| age (analysed women) | 21–43, M 31.8 (n = 99) | 19–40, M 27.2, SD 4.4 (n = 270) |
| collection | NFP trial 2008–2011, Marquette method | CrM charting, CMFS 1990–97 (162 women), TTP 2003–06 (9), CEIBA 2009–13 (99) |
| Week anchors (H) | 22, 29, 36, 43, 50 | 22, 29, …, 99 |
| observed counts (Haf, Dil, Wk, W-D, DiD) | 26, 61, 32, 18, 27 | **41, 76, 35, 40, 41** |
| women establishing each pattern | 20, 48, 26, 18, 23 | 38, 66, 33, 35, 39 |
| Haflaga runs at H = 30; DiD runs with |d| = 1 | 3; 16 | 8; 18 |

## 2. Marginal tests (two-sided Monte Carlo p; Holm over five patterns within each null)

| Pattern | H_G μ (σ) | z | p | Holm | H_iid μ (σ) | z | p | Holm | H_W μ (σ) | z | p | Holm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Haflaga (41) | 12.42 (3.43) | +8.33 | <.001 (0/50,000 ≥ 41) | .0002 | 12.56 (3.47) | +8.20 | <.001 | .0002 | 35.60 (4.78) | +1.13 | .313 | 1 |
| Dilug (76) | 54.95 (6.99) | +3.01 | .0056 | .021 | 54.84 (7.05) | +3.00 | .0052 | .021 | 75.71 (7.75) | +0.04 | 1.00 | 1 |
| Week (35) | 27.09 (4.18) | +1.89 | .083 | .083 | 27.17 (5.05) | +1.55 | .156 | .156 | 35.61 (3.97) | −0.15 | .986 | 1 |
| Week-Dilug (40) | 27.40 (4.17) | +3.02 | .0052 | .021 | 27.52 (5.08) | +2.46 | .025 | .050 | 38.56 (4.10) | +0.35 | .812 | 1 |
| Dilug-in-Dilug (41) | 27.08 (5.04) | +2.76 | .012 | .024 | 27.09 (5.09) | +2.73 | .014 | .042 | 41.48 (5.95) | −0.08 | 1.00 | 1 |

## 3. Joint five-pattern Mahalanobis test (primary) and checks

| Null | D_M | p (replicates ≥ D_M) | split-batch D_M, p | 3-pattern D_M, p (exploratory) | cond(Σ̂) | max |r| |
|---|---|---|---|---|---|---|
| H_G | 9.99 | <.001 (0 of 50,000) | 10.03, <.001 | 9.10, <.001 | 4.52 | .11 |
| H_iid | 9.53 | <.001 (0) | 9.52, <.001 | 8.81, <.001 | 4.45 | .14 |
| H_W | 1.20 | .921 (46,057) | 1.20, .922 | 1.17, .715 | 3.97 | .12 |

Fehring: 5.08 / 4.87 / 2.93 (p .126); split 5.07 / 4.87 / 2.95 (.119); 3-pattern 4.50 / 4.39 / 2.68 (.067).
H_W covariance and correlation matrices are in `summary.json` (`W.covariance`, `W.correlation`).

## 4. Sensitivity

- **Leave-one-woman-out (H_W, exact from per-woman counts):** joint p .861–.974 (Fehring .055–.224);
  no marginal p below .05 for any omission (Haflaga .248–.467).
- **First-12 truncation:** 2,109 cycles; counts (40, 75, 34, 36, 41); H_W D_M = 1.34, p = .878
  (Fehring 1,289 cycles, 2.49, .289).
- **Harmonised support 18 ≤ L ≤ 54:** removes 28 cycles; 258 women, 2,050 cycles, no woman split
  (every out-of-range cycle sits at a sequence end); counts unchanged (41, 76, 35, 40, 41);
  H_G D_M 9.27 (p <.001), H_iid 8.85 (<.001), H_W 1.15 (.933). Not the Fehring inclusion criterion —
  a support restriction only.
- **Follow-up influence:** no woman has ≥ 20 cycles; the top follow-up quintile (n_i ≥ 11: 56 women)
  holds 31.7 % of cycles and 56 % of Haflaga events (Fehring: n_i ≥ 14, 28 women, 36.9 % of cycles,
  50 % of Haflaga events). Spearman(n_i, events) = .30, .18, .16, .17, .18.
- **Tracking-quality:** the deposit carries no skipped-tracking flag; the 186 date gaps (each an
  excluded cycle in the source) are handled by splitting, never bridging; 26 cycles > 54 days remain
  in the primary analysis and are removed in the harmonised sensitivity with no change in conclusion.

## 5. Simulation validation (same sizes as the paper)

| | H_G | H_iid | H_W |
|---|---|---|---|
| marginal size (5 patterns) | .039 .042 .030 .032 .046 | .028 .039 .038 .036 .048 | .045 .042 .044 .036 .035 |
| Holm family | .039 | .038 | .034 |
| joint plug-in / split, α = .05 | .051 / .051 | .052 / .046 | .049 / .051 |
| joint plug-in / split, α = .01 | .010 / .008 | .009 / .009 | .011 / .011 |

MC SE at .05 with 10,000 datasets ≈ .002. Power of H_W (400 datasets × 2,000 permutations):

| Alternative | Haf | Dil | Wk | W-D | DiD | Holm-any | Joint |
|---|---|---|---|---|---|---|---|
| AR(1) ρ = 0 | .030 | .048 | .035 | .032 | .035 | .030 | .045 |
| ρ = .25 | .135 | .278 | .110 | .060 | .132 | .223 | .443 |
| ρ = .50 | .682 | .818 | .358 | .200 | .465 | .877 | .993 |
| ρ = .75 | 1 | .992 | .778 | .595 | .585 | 1 | 1 |
| persistence q = 0 | .045 | .025 | .030 | .058 | .035 | .028 | .058 |
| q = .10 | .960 | .180 | .540 | .360 | .058 | .938 | .970 |
| q = .20 | 1 | .422 | .942 | .782 | .070 | 1 | 1 |
| q = .30 | 1 | .718 | .988 | .908 | .362 | 1 | 1 |

Fehring joint power: .507 / .995 / 1 (AR ρ .25/.5/.75); .983 / 1 / 1 (q .1/.2/.3). Comparable.

## 6. Exploratory pairwise (SI-type analyses)

Per-woman rates, Spearman, permutation p, Holm over ten pairs: only **Haflaga × Week** is non-null
(r_s = .27, p = .0001, Holm .001; binary version χ² = 19.9, Holm .001). A Haflaga run at a Week anchor
value necessarily produces a Week event (containment analogous to Haflaga-30 ⊆ Week-Dilug), so this
is partly structural. Fehring's Haflaga × Dilug-in-Dilug association (r_s = .27, Holm .038) **does not
replicate** (r_s = −.09, p = .14; partial r given n_i = −.14; follow-up-stratified p = .21).
Haflaga × Week-Dilug: r_s = .09 (Fehring .23). SD comparison: 38 Haflaga women, within-woman SD of H
1.89 (1.37) vs 124 no-pattern women 4.61 (4.47); Welch t = −5.91, d = −0.68 — expected by construction.

## 7. Answer to the replication question

Yes. In an independent cohort (different method, sites, decades, investigators), the count vector is
displaced even further from the pooled nulls than in Fehring (D_M ≈ 10 vs 5; all five counts above
their global means, Haflaga z = 8.3), and the displacement vanishes completely once each woman's own
multiset is conditioned on (D_M = 1.20, p = .92; every |z| ≤ 1.13; LOO .86–.97; truncation .88;
harmonised support .93). Fehring's one borderline H_W feature (Dilug-in-Dilug deficit, Holm .19) is
absent here (z = −0.08), and the exploratory Haflaga × Dilug-in-Dilug between-woman association does
not replicate. The result is consistent with between-woman heterogeneity in cycle-length regularity
as the source of the global signal; as before, this explanation is supported, not uniquely identified.

## 8. Files

Everything lives in `external_replication/`: `scripts/` (analysis_engine.py — self-contained copy of the
analysis as functions of a Dataset, verified bit-identical to the cached Fehring results; run_external.py;
external_figures.py; prepare_utah.py), `data/` (raw deposit, sequences.csv, women.csv, preprocessing.json,
README.md), `results/` (null_replicates.npz, summary.json, sensitivity_L18-54.json, sim_validation.npz,
comparison.md/json, provenance.json, run.log, figures/figE1–E6), `reports/` (this file, the preprocessing
and deviation log, the draft manuscript subsection).
The only file changed outside the folder is `scripts/patterns.py` (Week anchors defined by the calendar
rule; Fehring anchors, counts and all cached results unchanged). Environment: Python 3.11.15,
numpy 2.4.6, scipy 1.17.1, pandas 3.0.5, matplotlib 3.11.2 (`results/provenance.json` records seeds, B,
hashes and commit).
