"""Build the analysed external sequence file from the University of Utah deposit
"Menstrual Cycles Length of Women in the USA and Canada, 1990-2013"
(Stanford & Najmabadi 2023; The Hive, DOI 10.7278/S50d-4gxs-s4hj; CC BY-NC).
Raw file: data/external/utah/raw/CrMcyclelength_share.csv (columns new_id, age, cycle_number,
cycle_start_date, cycle_end_date, cycle_length, conception_cycle).

Rules (pre-specified 2026-09-16; two refinements fixed on 2026-09-17 after inspecting the file
structure and before any pattern count was computed):
  1. Cycle length L = recorded ``cycle_length`` (days from the first day of menses to the day
     before the next menses). The script asserts that L equals end date - start date + 1 for
     every recorded cycle, so L is identical to the pre-specified "next start - start" for
     consecutive cycles and additionally retains each woman's last recorded cycle, as in the
     Fehring file. Haflaga value H = L + 1 (scripts/patterns.py).
  2. A cycle is usable when ``cycle_length`` is recorded. The 180 conception cycles have no
     recorded length (pregnancy, no menstrual end) and are excluded; they are always a woman's
     last record. Seven cycles have conception_cycle = "Missing" but a recorded length that
     matches the dates; they are kept as ordinary cycles (refinement).
  3. Two successive usable cycles are consecutive when the next start equals this end + 1.
     A woman's record is split at every gap (a gap always corresponds to at least one
     excluded cycle); gaps are never bridged.
  4. Each run of >= MIN_CYCLES (5) consecutive cycles is eligible; a woman contributes her
     single longest eligible run (earliest on ties), so each woman is one contiguous series.
  5. No cycle-length cutoff. Exact duplicate (woman, start date) rows would be removed and
     counted (there are none).
  6. Output schema = analysed columns of data/FilteredData.csv (ClientID, CycleNumber,
     LengthofCycle), CycleNumber contiguous within every woman; women.csv gives age and cohort
     era for the retained women.

Usage: python scripts/external/prepare_utah.py [RAW_CSV]
Writes data/external/utah/sequences.csv, women.csv and preprocessing.json.
"""
import hashlib, json, os, sys
from collections import Counter
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT_DIR = os.path.join(ROOT, "data", "external", "utah")
RAW = os.path.join(OUT_DIR, "raw", "CrMcyclelength_share.csv")
MIN_CYCLES = 5


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(raw_path=RAW):
    raw = pd.read_csv(raw_path)
    raw.columns = [c.lstrip("﻿") for c in raw.columns]
    d = raw.copy()
    d["start"] = pd.to_datetime(d["cycle_start_date"], format="%m/%d/%y")
    d["end"] = pd.to_datetime(d["cycle_end_date"], format="%m/%d/%y")
    dups = int(d.duplicated(["new_id", "start"]).sum())
    d = d.drop_duplicates(["new_id", "start"]).sort_values(["new_id", "start"]).reset_index(drop=True)
    span = (d["end"] - d["start"]).dt.days + 1
    rec = d["cycle_length"].notna()
    assert (span[rec] == d.loc[rec, "cycle_length"]).all(), "recorded length != end - start + 1"
    assert (d.loc[~rec, "conception_cycle"] == "Yes").all(), "blank length outside conception cycles"

    log = dict(source_file=os.path.basename(raw_path), source_sha256=sha256(raw_path),
               n_rows_raw=int(len(raw)), n_women_raw=int(raw["new_id"].nunique()),
               duplicate_rows_removed=dups, min_cycles=MIN_CYCLES,
               conception_cycles_excluded=int((~rec).sum()),
               conception_cycles_all_terminal=bool(int((d.groupby("new_id")["conception_cycle"].last() == "Yes").sum()) == int((~rec).sum())),
               cycles_with_missing_conception_flag_kept=int(((d["conception_cycle"] == "Missing") & rec).sum()),
               recorded_length_equals_dates="asserted",
               rules="see module docstring")

    seqs, women_rows = [], []
    n_gap_breaks = n_gap_with_cn_jump = n_cn_jump_without_gap = 0
    run_hist = Counter()
    for wid_, g in d.groupby("new_id", sort=True):
        g = g.reset_index(drop=True)
        runs, cur = [], []
        for k in range(len(g)):
            if pd.isna(g.loc[k, "cycle_length"]):          # conception cycle: unusable, ends the run
                if cur: runs.append(cur); cur = []
                continue
            if cur:
                gap = (g.loc[k, "start"] - g.loc[k - 1, "end"]).days - 1
                cn_jump = int(g.loc[k, "cycle_number"] - g.loc[k - 1, "cycle_number"])
                if gap != 0:
                    n_gap_breaks += 1
                    n_gap_with_cn_jump += int(cn_jump > 1)
                    runs.append(cur); cur = []
                elif cn_jump > 1:
                    n_cn_jump_without_gap += 1     # dates abut: treated as consecutive
            cur.append(int(g.loc[k, "cycle_length"]))
        if cur: runs.append(cur)
        for r in runs: run_hist[len(r)] += 1
        elig = [r for r in runs if len(r) >= MIN_CYCLES]
        if elig:
            best = max(elig, key=len)          # first maximum = earliest run on ties
            seqs.append((int(wid_), best))
            first_year = int(g.loc[0, "start"].year)
            era = "CMFS 1990-1997" if first_year <= 1997 else ("TTP 2003-2006" if first_year <= 2006 else "CEIBA 2009-2013")
            women_rows.append(dict(ClientID=int(wid_), age=int(g.loc[0, "age"]), first_cycle_year=first_year, cohort_era=era,
                                   n_records=int(len(g)), n_runs=len(runs), n_eligible_runs=len(elig), n_cycles_analysed=len(best)))
    log.update(gap_breaks=n_gap_breaks, gap_breaks_with_cycle_number_jump=n_gap_with_cn_jump,
               cycle_number_jumps_without_date_gap=n_cn_jump_without_gap,
               run_length_histogram={int(k): v for k, v in sorted(run_hist.items())},
               women_with_eligible_run=len(seqs), women_without_eligible_run=int(d["new_id"].nunique() - len(seqs)),
               women_with_several_eligible_runs=int(sum(w["n_eligible_runs"] > 1 for w in women_rows)))
    rows = [(wid_, k, L_) for wid_, s in seqs for k, L_ in enumerate(s, 1)]
    out = pd.DataFrame(rows, columns=["ClientID", "CycleNumber", "LengthofCycle"])
    women = pd.DataFrame(women_rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_csv(os.path.join(OUT_DIR, "sequences.csv"), index=False)
    women.to_csv(os.path.join(OUT_DIR, "women.csv"), index=False)
    log.update(n_cycles_analysed=int(len(out)), n_women_analysed=int(out["ClientID"].nunique()),
               cycles_per_woman=dict(mean=round(float(women.n_cycles_analysed.mean()), 2), sd=round(float(women.n_cycles_analysed.std()), 2),
                                     min=int(women.n_cycles_analysed.min()), max=int(women.n_cycles_analysed.max())),
               L_min=int(out["LengthofCycle"].min()), L_max=int(out["LengthofCycle"].max()),
               n_cycles_L_below_18=int((out["LengthofCycle"] < 18).sum()), n_cycles_L_above_54=int((out["LengthofCycle"] > 54).sum()),
               age=dict(min=int(women.age.min()), max=int(women.age.max()), mean=round(float(women.age.mean()), 1), sd=round(float(women.age.std()), 1)),
               cohort_era_counts=women.cohort_era.value_counts().to_dict(),
               output_sha256=sha256(os.path.join(OUT_DIR, "sequences.csv")))
    with open(os.path.join(OUT_DIR, "preprocessing.json"), "w") as fh:
        json.dump(log, fh, indent=1)
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else RAW)
