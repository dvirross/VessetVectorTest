# Utah external replication — dataset selection, preprocessing record and deviation log (2026-09-17)

## Dataset selection
Requirements: individual-level, per-woman consecutive cycle lengths in whole days with verifiable
order and gap information; at least 5 cycles per woman for a useful share of women; reproductive-age
women; provenance independent of the Marquette trial; open or documented access. The University of
Utah deposit below was the only openly downloadable source meeting all of them (the Colombo
databases at Padua meet them but require an application; the Tremin/BIMORA calendars are
perimenopausal; app cohorts are under data-use agreements; mcPHASES has too few cycles per
participant).

## Source
Stanford, J. B., & Najmabadi, S. (2023). *Menstrual Cycles Length of Women in the USA and Canada,
1990–2013*. The Hive, University of Utah. https://doi.org/10.7278/S50d-4gxs-s4hj (CC BY-NC).
Downloaded by the author from Hive and uploaded to the session on 2026-09-17 (the host is blocked
from the sandbox). Files in `external_replication/data/raw/`:

| file | SHA-256 |
|---|---|
| hive_download.zip | ae0c63c3584c248bc10b141f49dbd9df34f550d56c8aff30dbb15471d455d71c |
| CrMcyclelength_share.csv | ad67a3eed731a038d6502a1eb42a9dcddf8f6c69a607a61b85f3798198043dac |
| Najmabadi_README20230824.txt | 4bbb2c93efada3c84fe841bc5ef82820c30f8c1dce7e7ebc07d98558138a4678 |

## What is verifiable in the file
3,324 rows, 581 women; columns new_id, age (18–40, at first cycle), cycle_number, cycle_start_date,
cycle_end_date, cycle_length (blank for the 180 conception cycles), conception_cycle (No 3,137;
Yes 180; Missing 7). Recorded cycle_length equals end − start + 1 for every recorded cycle
(asserted). Start dates 1990-01-21 to 2012-10-13; year distribution identifies the three cohorts
(CMFS 1990–97, TTP 2003–06, CEIBA 2009–13). No duplicate rows. 2,743 within-woman successive
pairs: 2,557 abut (next start = end + 1), 186 have a gap, none overlap; gaps coincide with jumps in
cycle_number (an excluded cycle), so the deposit already removes some cycles. All 180 conception
cycles are a woman's last record. No contraception, pregnancy-history or procedure variables beyond
the README eligibility statement (regular bleeding, not pregnant, not exclusively breastfeeding,
recent OC users after 1–2 bleeds).

## Preprocessing applied (`external_replication/scripts/prepare_utah.py`)
1. L = recorded cycle_length (identical to next start − start for consecutive cycles; keeps each
   woman's last recorded cycle). H = L + 1.
2. Conception cycles (no recorded length) excluded (180); the 7 cycles with flag "Missing" but a
   recorded length that matches the dates kept as ordinary cycles.
3. Split at every gap (170 gap breaks among usable cycles; 169 with a cycle-number jump). One
   cycle-number jump without a date gap treated as consecutive (dates abut).
4. Longest run of ≥ 5 consecutive cycles per woman (earliest on ties); 7 women had more than one
   eligible run.
5. No cycle-length cutoff.

Result: **270 women, 2,139 cycles**; 5–15 cycles per woman (mean 7.92, SD 2.68, median 7);
L 17–98 (mean 29.98, SD 6.16, median 29; 2 cycles < 18, 26 cycles > 54; 12.0% ≥ 35);
age 19–40 (mean 27.2, SD 4.4); cohort era: CMFS 162, CEIBA 99, TTP 9 women.
`sequences.csv` SHA-256 08dddb0c04c6a107b95eefc6ef081fabedfe2e22bd8c9832ce6e1af45a9ee37f.

## Deviations from the 2026-09-16 pre-specification (all fixed before any null was run)
- Recorded length instead of start-date difference (identical where both exist; retains the last
  cycle of each woman, as in Fehring's file).
- "Missing" conception flag with recorded length → kept (7 cycles).
- **Week anchor definition corrected in `scripts/patterns.py`.** The manuscript defines Week as a
  repeated haflaga at a same-weekday value (L a multiple of 7: H = 22, 29, 36, …). The code
  derived the anchors from the data range, `range(H_min + 3, H_max + 1, 7)`, which equals the
  manuscript set on Fehring (H_min = 19) but gave {21, 28, …} on Utah (H_min = 18). The function now
  returns {H : (H − 1) mod 7 = 0} within the observed range. Fehring anchors, counts and all cached
  results are unchanged. Utah observed Week count: 44 under the old range-derived set, **35** under
  the manuscript rule; the latter is used. The manuscript's shorthand
  "$\mathcal{W} = \{H_{\min}+3, H_{\min}+10, \ldots\}$" (Section 2, Week item) should be reworded at
  revision to the calendar definition; no number in the paper changes.

## Applicability notes
- First-12 truncation applies (women have up to 15 cycles).
- Harmonised 18–54 sensitivity removes 28 cycles and splits sequences at those points.
- Validation at full size (10,000 fresh datasets vs 20,000 replicates; 8 × 400 × 2,000 power)
  is feasible (≈ 40 min on 3 cores) and is run unreduced.
