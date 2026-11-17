# Utah replication: preprocessing verification and sensitivity (2026-11-17)

**Purpose.** Verify that the Utah replication was preprocessed correctly, with emphasis on the
assumptions behind the order of operations, and compare each decision with the one taken on the
Marquette data, where the equivalent decisions shaped the whole analysis. Two independent read-only
audits were run (one on the Marquette pipeline, one re-deriving the Utah pipeline from the raw deposit
without importing any project code), followed by a twelve-specification sensitivity analysis
(`scripts/preprocessing_sensitivity.py` → `results/preprocessing_sensitivity.json`; B = 50,000,
seed 17, same engine as the published analysis). This note extends the 2026-11-01 audit and the
run-selection sensitivity of 2026-11-15 (SI Table S8). Nothing in the published analysis was changed.

## 1. What each pipeline actually does

| Decision | Marquette (`scripts/filter_raw.py`, verified against `Filtering_original.ipynb`; output identical) | Utah (`scripts/prepare_utah.py`, re-derived row for row) |
|---|---|---|
| Duplicate records | 3 women with duplicated (ID, CycleNumber) rows dropped (16 rows); nfp8107's first copy of cycles 1–8 dropped (8 rows), second copy kept as a correction (decision taken with the supervisor, dissertation p. 17) | none found (0 duplicate rows, 0 duplicate (id, start)) |
| Minimum cycles per woman | ≥ 5 recorded rows (38 women, 87 rows removed) | ≥ 5 **consecutive** cycles in a single run |
| Gap handling / run selection | **none possible**: the archive has no cycle dates (80 columns, none a cycle start or end); "consecutive" rests on the archive's own `CycleNumber` (1..n in every woman; 0 skipped numbers in the raw file) | split at every date gap (next start ≠ previous end + 1); keep the longest run ≥ 5, earliest on ties |
| Length definition | `LengthofCycle` as archived, taken on trust as onset-to-day-before-next-onset; H = L + 1 | recorded `cycle_length` = end − start + 1 for all 3,144 recorded cycles; in all 2,557 abutting pairs next start = end + 1, so L is the onset-to-onset interval; H = L + 1 is the same quantity |
| Length-range filter | none; archive already 18–54 (0 cycles outside) | none; 17–98 kept (28 cycles outside 18–54 in 18 women; sensitivity in SI S4) |
| Age filter | none (99/118 with age, 21–43; one woman aged 43 retained) | none (deposit 18–40; analysed 19–40) |
| Conception / pregnancy cycles | no identifiable flag; nothing excluded (`ReproductiveCategory` codes undocumented) | 180 conception cycles (no length) excluded; all are the woman's last record; 7 'Missing' flags kept (all terminal; 6 analysed) |
| Missing values | none in `LengthofCycle` | blank length only in conception cycles |

The two pipelines agree on the minimum of five cycles, the absence of length and age filters, and
H = L + 1. The Utah gap-split rule has no Marquette counterpart and could not have had one; the Utah
conception exclusion has no Marquette counterpart because the archive carries no such flag.

## 2. Findings of the Utah re-derivation (independent code, pandas only)

Confirmed: 3,324 rows / 581 women; 180 conception cycles, all terminal; 7 'Missing' flags, all
terminal; 2,743 successive pairs = 2,557 abutting + 186 gaps + 0 overlaps (170 gaps between two
recorded cycles, 16 whose second member is the conception cycle); 732 runs; 270 women / 2,139 cycles,
identical to `data/sequences.csv`; cycles per woman 5–15; ages 19–40.

New facts that matter for interpretation:

- **Gaps are whole missing cycles.** 185 of 186 gaps coincide with a jump in `cycle_number`; gap
  length per skipped number has median 29 days and lies in 18–54 in 177 of 185. The one exception
  is a 1-day gap (woman 42, 4 recorded cycles, not analysed). One `cycle_number` jump has abutting
  dates (woman 517, 5 → 7; analysed; treated as consecutive).
- **The deposit is already a filtered subset of the parent numbering.** 416 of 581 women start at
  `cycle_number` > 1 (859 leading numbers absent) and 285 internal numbers are absent; the README does
  not say why. Gap rates differ by cohort era: CEIBA 10.6 % of recorded-to-recorded pairs, CMFS 3.9 %,
  TTP 3.4 %.
- **What the gap rule alone costs.** 28 women with ≥ 5 recorded cycles but no single run ≥ 5 (169
  cycles; longest run 4 in 16 of them) are excluded only because a gap split their record. A further
  16 women have 5 records but only 4 recorded cycles and are lost only because the conception cycle has
  no length. The 7 women with two eligible runs (37 cycles) were already covered by SI Table S8.
- **Terminal cycles.** The end date of a woman's last recorded cycle is not fixed by a following onset.
  Last recorded non-conception cycles are longer and more dispersed (n = 401; mean 32.1, SD 9.9, max 98)
  than other recorded cycles (mean 30.1, SD 6.1), consistent with some being end-of-observation spans.
  181 of the 270 analysed runs end at the woman's last record. Unverifiable from the file.
- **The 7 'Missing' flags** are all terminal, like conception cycles, and plausibly denote an unknown
  pregnancy outcome; keeping them is defensible but not neutral.

Findings on the Marquette side that bear on the comparison: the archive's per-woman
`MeanCycleLength` equals the mean of the released rows for 129 of 137 women with the field; for 5 women
no contiguous window or leave-one-out explains the archived mean, so the archive's summary may have been
computed on a cycle set that differs from the released rows. Whether unrecorded cycles sit between
consecutive `CycleNumber`s in Marquette is therefore **not verifiable**, which is exactly the situation
the `no_split` specifications below emulate for Utah.

## 3. Sensitivity of the Utah results to the preprocessing decisions

All specifications: three nulls with B = 50,000, seed 17; H_W permutes within each analysed run; the
five-pattern D_M with its Monte Carlo p. H_G and H_iid give p < .001 in every row.

| specification | runs / women / cycles | observed (Haf, Dil, Wk, W-D, DiD) | D_M H_G / H_iid / H_W | p (H_W) | z under H_W | min Holm p (H_W) |
|---|---|---|---|---|---|---|
| baseline (published) | 270 / 270 / 2,139 | 41, 76, 35, 42, 41 | 10.04 / 9.56 / 1.22 | .916 | 1.13, 0.04, −0.15, 0.41, −0.08 | 1 |
| ties → latest run | 270 / 270 / 2,139 | 40, 76, 34, 42, 41 | 9.89 / 9.40 / 1.18 | .926 | 1.05, 0.08, −0.08, 0.50, −0.01 | 1 |
| first run ≥ 5 instead of longest | 270 / 270 / 2,139 | 41, 76, 35, 42, 41 | 10.04 / 9.56 / 1.22 | .916 | identical to baseline | 1 |
| gaps ≤ 7 days not split | 270 / 270 / 2,139 | 41, 76, 35, 42, 41 | 10.04 / 9.56 / 1.22 | .916 | identical to baseline | 1 |
| gaps ≤ 45 days not split (≈ one missing cycle bridged) | 285 / 285 / 2,403 | 47, 92, 39, 50, 50 | 11.04 / 10.41 / 1.70 | .716 | 1.35, 0.44, −0.17, 0.72, 0.29 | 1 |
| no split at all; recorded lengths (Marquette-equivalent treatment) | 298 / 298 / 2,538 | 51, 102, 43, 51, 52 | 11.81 / 11.21 / 2.11 | .485 | 1.61, 0.97, −0.10, 0.64, 0.26 | .68 |
| no split; L = onset-to-onset across gaps | 298 / 298 / 2,538 | 41, 83, 42, 48, 41 | 10.23 / 9.52 / 1.15 | .934 | 0.65, 0.48, 0.54, 0.55, −0.01 | 1 |
| 'Missing'-flag cycles dropped | 269 / 269 / 2,129 | 41, 76, 35, 42, 41 | 10.17 / 9.55 / 1.31 | .888 | 1.21, 0.06, −0.09, 0.48, −0.03 | 1 |
| all runs ≥ 4 (window of Dilug-in-Dilug) | 351 / 325 / 2,472 | 45, 85, 42, 51, 43 | 10.88 / 10.24 / 1.39 | .859 | 1.12, 0.08, 0.40, 0.76, −0.16 | 1 |
| last recorded cycle of every woman dropped | 234 / 234 / 1,814 | 35, 67, 28, 34, 35 | 8.95 / 8.57 / 1.16 | .931 | 0.83, 0.22, −0.66, 0.31, −0.08 | 1 |
| last record dropped and no split | 264 / 264 / 2,158 | 41, 91, 35, 37, 45 | 10.33 / 9.85 / 1.85 | .639 | 0.89, 1.34, −0.66, −0.17, 0.29 | 1 |
| also split at the abutting number jump (woman 517) | 270 / 270 / 2,135 | 41, 76, 35, 42, 41 | 10.06 / 9.53 / 1.19 | .922 | 1.11, 0.03, −0.16, 0.36, −0.10 | 1 |

Together with SI Table S8 (all runs ≥ 5 within-run 1.11 / .943; pooled within-woman 1.19 / .923; all
runs ≥ 3 1.43 / .843), the within-woman D_M ranges from 1.11 to 2.11 and its p from .49 to .94 across
fifteen specifications; no marginal H_W test reaches Holm p < .68 in any of them. The specification that
moves the result most is the Marquette-equivalent one (ignore gaps, recorded lengths), which creates
artificial adjacency across 185 missing cycles and still leaves the count vector inside its conditional
null; replacing the recorded length by the onset-to-onset span, which is the correct quantity when a
gap is bridged, returns the result to the baseline (1.15, p = .934).

## 4. Assessment

- The implemented pipeline does what the deviation log says, in the stated order (duplicates →
  conception cycles → gap split → run selection), and every claimed number re-derives.
- The decisions are more conservative than Marquette's, not less: Marquette could not detect gaps and
  so implicitly counts across any unrecorded cycle; Utah never does. The one direction in which the Utah
  rule is stricter (28 women lost to the gap rule alone, 16 to the conception exclusion) does not change
  any conclusion when relaxed.
- The two unverifiable items are the same in both cohorts and cannot be settled from the files: whether
  `LengthofCycle` / `cycle_length` follow the same onset convention (both are consistent with it), and
  whether unrecorded cycles sit between consecutive Marquette rows.
- Open author decisions (unchanged from the 2026-11-01 audit): the halachic semantics of counting a
  pattern across a gap in the released file; whether to pursue the separately released CEIBA data
  (https://hive.utah.edu/records/kfcm2-ynw54); whether the 7 'Missing' flags should be treated like
  conception cycles (numerically immaterial).
- Suggested for the paper, if desired: one sentence in SI Section S4 pointing to this table, or an
  extension of Table S8 by the two rows that matter most (no split with recorded lengths; last record
  dropped). No change to the main text is required by these results.

## Files

`scripts/preprocessing_sensitivity.py`, `results/preprocessing_sensitivity.json` (this folder);
`scripts/all_runs_sensitivity.py`, `results/all_runs_sensitivity.json` (2026-11-15);
`reports/2026-11-01_preprocessing_audit_and_open_questions.md`;
`reports/utah_preprocessing_and_deviations.md`; Marquette: `scripts/filter_raw.py`,
`scripts/Filtering_original.ipynb`, `data/RawData.csv` (no date columns).
