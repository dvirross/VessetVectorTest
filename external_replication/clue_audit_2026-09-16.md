# External replication — stage 1 audit (Clue), 2026-09-16

Status: **STOPPED at the agreed checkpoint** (repository audit, Fehring reproduction, Clue
download, Clue file audit). No Clue Monte Carlo has been run. No analysis code or cached
Fehring result has been modified.

## 1. Repository audit (authoritative implementation)

| Component | Location | Notes |
|---|---|---|
| Raw → analysed Fehring data | `scripts/filter_raw.py` (`filter_raw`) | 1,665/159 → drop 3 duplicated-record women (16 rows) → ≥5 cycles (−38 women, −87 rows) → drop first copy of nfp8107 (8 rows) → **1,554 cycles / 118 women**; asserts contiguous CycleNumber |
| Five counting rules, vectorised | `scripts/patterns.py` `Counter.run_starts/count/count_by_woman` | H = L + 1; maximal contiguous runs; never cross woman boundary; Week anchors = {min+3, min+10, …}; Week-Dilug anchor H = 30 |
| Five counting rules, reference loop | `patterns.count_all_loop` | mirrors the notebook definitions |
| Equivalence test | `patterns.verify_equivalence` | observed data + 3,000 × 3 nulls = 9,000 sequences |
| H_G generator | `patterns.perm_global` | pooled permutation without replacement, n_i preserved by position |
| H_iid generator | `patterns.multinomial_iid` | i.i.d. draws from pooled empirical distribution |
| H_W generator | `patterns.perm_within` | independent permutation inside each woman's block |
| Two-sided Monte Carlo p (+1 convention) | `patterns.two_sided_p` | min{1, 2·min((b⁺+1)/(B+1), (b⁻+1)/(B+1))} |
| Holm | `patterns.holm_adjust` (and a local copy `holm` in `postprocess.py`; identical algorithm) | |
| Five-pattern Mahalanobis, plug-in and split-batch | `patterns.mahalanobis_test(null, obs, ref=None)` | split = `ref=arr[:B/2]`, `null=arr[B/2:]` |
| Three-pattern exploratory | `postprocess.py` (`arr[:, [0,3,4]]`) and notebook fig 9 | |
| Null replicate run | `scripts/rerun_nulls.py` → `results/null_replicates.npz` | B = 50,000, seed 17, one RNG stream shared by the three nulls; stores per-woman H_W counts (B×118×5 int8) |
| All reported statistics | `scripts/postprocess.py` → `results/summary.json` | marginal, joint, split, 3-pattern, run details, SD comparison, influence, LOO, pairwise |
| Size (Part A), power (Part B), truncation (Part C) | `scripts/sim_validation.py` → `results/sim_validation.npz` | B_REF 20,000; N_FRESH 10,000; R_ALT 400; B_IN 2,000; AR(1) ρ ∈ {0,.25,.5,.75}; persistence q ∈ {0,.1,.2,.3}; first-12 truncation with 20,000 H_W replicates |
| Leave-one-woman-out | `postprocess.py` (uses additive `within_by_woman` counts — no re-simulation) | |
| Pairwise SI analyses | `postprocess.py` | rate-based Spearman + permutation p (50,000, seed 17) + Holm over 10 pairs; binary ever/never (Fisher + permutation χ²); follow-up-quintile-stratified permutation and partial Spearman for the two structurally expected pairs |
| Tables/figures | `scripts/build_notebook.py` → `notebooks/analysis.ipynb` (fig1, 3–10, A1); `paper/figures/gen_consort_flow.py` (fig2) | |

Unequal-follow-up sensitivity in the repository is **first-12 truncation + LOO + influence
share of women with ≥ 20 cycles**; there is no ≥8/≥10-cycle analysis (consistent with the
manuscript).

**Code vs manuscript disagreements found: none.** Every number in `paper/vesset_stat.tex`
that was checked (preprocessing counts; Tables 1–2; joint, split-batch and 3-pattern D_M and
p; LOO ranges; Part A/B/C values; SI pairwise values) matches `results/summary.json` and
`results/sim_validation.npz`. Minor notes only: `sim_validation.py` uses `Pool(3)` while
CLAUDE.md says "4 cores"; `postprocess.py` re-implements Holm locally instead of importing
`holm_adjust` (same algorithm).

## 2. Fehring reproduction (all passed, this session)

- `python scripts/filter_raw.py data/RawData.csv` → step counts exactly as documented;
  output `equals()` the released `data/FilteredData.csv`
  (SHA-256 `508ee88e…7811` as stated in README/DAS; raw `9fa76583…2a89`).
- `python scripts/patterns.py` → 118 women, 1,554 cycles, contiguous; vectorised = loop on
  observed data and 9,000 null sequences; observed counts (26, 61, 32, 18, 27).
- Fresh `rerun_nulls.py 50000 17` in a scratch copy (35 s): the `observed`, `perm`, `multi`,
  `within`, `within_by_woman` arrays are **bit-identical** to `results/null_replicates.npz`.
  Hence all benchmark values follow: H_G Haflaga 11.31 / 3.27 / z 4.50 / p .0002 / Holm .001;
  Dilug 48.92 / 6.57 / z 1.84 / p .085 / Holm .34; Week 27.02; W-D 16.25; DiD 27.17;
  D_M 5.08 (p .0001, 6 exceed), H_iid 4.87 (17), H_W 2.93 (p .1257; split 2.95 / .1186);
  3-pattern 4.50 / 4.39 / 2.68 (p .0666); LOO joint .055–.224.
- Cached validation: size marginal 2.9–4.9 %, Holm 3.2–3.7 %, joint 4.8–5.3 %; Part C
  1,289 cycles, (19, 51, 21, 16, 21), D_M 2.49, p .289.
- Environment: Python 3.11.15, numpy 2.4.6, scipy 1.17.1, pandas 3.0.5, matplotlib 3.11.2.

## 3. Clue source, checksum, retrieval

- Repository: https://github.com/iurteaga/menstrual_cycle_analysis (branch `master`,
  no tags). HEAD at retrieval: `5f21f094dbdb9a6c302424b66a52b7c38e99f683` (2023-07-07).
- File: `data/cycle_length_data/cycle_lengths.npz`, added in commit
  `c042ee4463d55148048a1adea5c55b7fb2714b69` (2021-07-29, "Updated with prediction work"),
  unchanged since.
- Raw URL: https://raw.githubusercontent.com/iurteaga/menstrual_cycle_analysis/5f21f094dbdb9a6c302424b66a52b7c38e99f683/data/cycle_length_data/cycle_lengths.npz
- Retrieved 2026-09-16T14:43Z by `git clone`. Size 175,030 bytes.
- **SHA-256 `ccf46330ad721e024ee5c25506c649dc3b37744563102040d99d01df1c32b83f`**.

## 4. Clue file audit — what is directly verifiable

Keys: `data_model`, `I`, `C`, `hyperparameters`, `cycle_lengths`, `cycle_skipped`,
`true_params`.

| Item | Value |
|---|---|
| `cycle_lengths` | int64, shape (5000, 11); no NaN, no padding, no zeros; every row has exactly 11 values |
| range | 10–142 days; mean 24.49, median 22; 7.4 % of cycles < 18, 1.6 % > 54; 2,370 of 5,000 rows entirely within 18–54 |
| IDs / age / contraception / pregnancy / procedure flags | **none present** |
| duplicates | 0 identical rows |
| ordering | columns are positional cycles 1..11; column means flat (24.3–24.6) |
| `data_model` | `'generalized_poisson'` |
| `I`, `C` | 5000, 11 |
| `hyperparameters` | [160, 4, 2, 20, 1, inf, 2, 20] = (κ, γ, α_ξ, β_ξ, ξ_max, x_max, α, β) |
| `cycle_skipped` | int64 (5000, 11), values 0–6; 9.3 % of cells > 0 |
| `true_params` | dict of per-row `lambda` (30.0–53.9), `xi` (−0.999 to −0.221), `pi` (0.0001–0.47) |

### The public file is a simulated dataset, not Clue user data

1. Its own metadata (`data_model`, `hyperparameters`, `true_params`, `cycle_skipped`) is
   exactly what `src/prediction/data_functions.py::generate_sim_model` writes with
   `np.savez_compressed` for a **generalized-Poisson simulation**; the loader in `get_data`
   comments that `true_params` is "Only possible if simulated data".
2. The stored hyperparameters and I = 5000, C = 11 coincide with the *simulated-data*
   example in `scripts/README.md` (`-data_model generalized_poisson -hyperparameters
   160 4 2 20 2 20 -I 5000 -C 11`; the code pads the 6-vector with ξ_max = 1, x_max = inf,
   giving the 8-vector in the file).
3. The rows are generated from the stored parameters: correlation of each row's mean
   non-skipped length with λ/(1−ξ) is 0.944; mean of x/(s+1) − λ/(1−ξ) is 0.002 (SD 0.77);
   mean within-row variance of non-skipped cycles 6.86 vs theoretical GP variance 6.88;
   the maximum non-skipped length is 43 and every cycle flagged `cycle_skipped ≥ 1` is ≥ 26;
   share of skipped cells 9.26 % vs mean π 9.11 %.
4. Its distribution is implausible for real menstruators (median 22 days, mean of
   non-skipped cycles 22.1) and does not match the published Clue cohort (mean 30.7, SD 7.9,
   median 29; MLHC 2021 Table 1).
5. Every associated paper's Data Availability Statement (npj Digital Medicine 2020, arXiv
   2102.12439, JAMIA 2022) says the Clue database "cannot be made directly available to the
   reader"; MLHC 2021 footnote 5 says "The synthetic dataset can be generated with the Python
   codebase". The scripts README's phrase "real data (contained in
   ../cycle_length_data/cycle_lengths.npz)" is contradicted by the file's contents.

### Separation requested in the task

A. **Documented about the published Clue cohort** (not verifiable in the file): 186,106
   users, 2,047,166 cycles, cycle length 30.7 ± 7.9 (median 29), users aged 21–33, no
   hormonal birth control, cycles "excluded by the user" removed, ≥ 11 cycles per user with
   first 10 used for training.
B. **Verifiable in the public file**: none of the above. The file contains no user
   identifiers, no age, no contraception, pregnancy or procedure information, and its
   generative origin is recorded in its own metadata. It cannot be linked to the cohort.

## 5. Preprocessing / sample-size proposal — conditional

An "independent external replication in the Clue cohort" **cannot be performed with this
file**: it would be a replication on synthetic data drawn from an exchangeable hierarchical
model (each user's cycles are i.i.d. given (λ_i, ξ_i, π_i)), so H_W is true by construction
and H_G/H_iid are false by construction (between-user heterogeneity in λ_i). Any result
would be a known-answer check of the engine, not evidence about menstrual data.

If the file were nevertheless analysed under the current Fehring rules, the sample would be:
all 5,000 rows retained (fixed 11 complete cycles each, ≥ 5 satisfied, no gaps, no
duplicates) → 55,000 cycles; first-12 truncation inapplicable (C = 11 < 12); harmonised
18–54 sensitivity: 2,370 users with all 11 cycles in range (26,070 cycles), or sequence
splitting at out-of-range cycles for the others; `cycle_skipped` gives a ground-truth
"skipped-period" flag for a tracking-artefact sensitivity.

Runtime with the existing engine (single core): H_W 12 ms/replicate → B = 50,000 in
≈ 10 min; H_G ≈ 1.5 min; H_iid ≈ 2.5 min; per-woman H_W counts for LOO 1.25 GB (int8);
Part A ≈ 0.5 core-h; Part B (8 × 400 × 2,000) ≈ 21 core-h (≈ 5–7 h on 4 cores).

## 6. Options for the author

1. **Request the real Clue data** under a data-use agreement (as the papers instruct) —
   the only route to an actual Clue replication.
2. **Replace the external dataset** with another public longitudinal cycle-length source,
   selected *before* any pattern counts are seen.
3. **Re-label the exercise** as a synthetic known-answer validation of the three nulls on
   the Li/Urteaga generative model (5,000 × 11), reported as such — never as a Clue
   replication.
4. Stop here.

## 7. Engine verification addendum (2026-09-17)

`scripts/analysis_engine.py` reproduces every cached Fehring result: null replicates under
all three nulls and the per-woman H_W counts are bit-identical (B = 50,000, seed 17); every
value in `results/summary.json` is unchanged; validation Parts A (three nulls), B (all eight
power scenarios, 400 x 2,000 each) and C (first-12 truncation) match the cached
`sim_validation.npz` exactly. `scripts/run_external.py` was smoke-tested end to end on the
Fehring sequences.
