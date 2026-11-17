"""Utah replication: sensitivity of the results to the preprocessing decisions (added 2026-11-17).

Complements all_runs_sensitivity.py (run-selection rule). Each specification rebuilds the analysed
sequences from the raw deposit with one decision changed, then runs the engine's three nulls
(B = 50,000, seed 17; H_W within each analysed run) and the five-pattern Mahalanobis test.

  baseline                   longest run >= 5 per woman, split at every date gap (the published analysis)
  tie_latest                 ties between equally long runs broken towards the latest run (baseline: earliest)
  first_run_ge5              the woman's FIRST run >= 5 instead of her longest
  gap_tol_7d                 a date gap of <= 7 days does not split a record (recording-artefact tolerance)
  gap_tol_45d                a gap of <= 45 days (about one missing cycle) does not split; recorded lengths kept
  no_split_recorded_L        gaps ignored entirely: each woman's recorded cycles in date order form one
                             sequence (the treatment the Marquette archive implicitly receives, since it
                             carries no dates); >= 5 recorded cycles per woman
  no_split_start_to_start    as above, but L = next onset - this onset, so a gap inflates the interval
                             instead of creating artificial adjacency (last cycle: recorded length)
  drop_missing_flag          the 7 cycles whose conception_cycle flag is 'Missing' end a run and are dropped
  all5_run_min4              every run >= 4 cycles (the window Dilug-in-Dilug needs), H_W within each run
  drop_last_record           each woman's last recorded cycle removed (its end date is not fixed by a
                             following onset and terminal lengths are longer and more dispersed)
  drop_last_record_no_split  the two previous ideas combined: no gap split and last record removed
  split_cn_jump_too          also split at the one cycle_number jump with abutting dates (woman 517)

Usage: python external_replication/scripts/preprocessing_sensitivity.py
-> external_replication/results/preprocessing_sensitivity.json (about 9 min)
"""
import json, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from all_runs_sensitivity import analyse, segment, baseline_runs, RAW, B, SEED  # noqa: E402

OUT = os.path.join(EXT, "results", "preprocessing_sensitivity.json")


def load_raw():
    raw = pd.read_csv(RAW)
    raw.columns = [c.lstrip("﻿") for c in raw.columns]
    d = raw.copy()
    d["start"] = pd.to_datetime(d["cycle_start_date"], format="%m/%d/%y")
    d["end"] = pd.to_datetime(d["cycle_end_date"], format="%m/%d/%y")
    return d.drop_duplicates(["new_id", "start"]).sort_values(["new_id", "start"]).reset_index(drop=True)


def runs_with(d, gap_tol=0, drop_missing_flag=False, start_to_start=False, drop_last_record=False,
              split_cn_jump=False):
    """Segment like prepare_utah.py, with one decision changed. gap_tol=None: never split at a gap."""
    runs = []
    for w, g in d.groupby("new_id", sort=True):
        g = g.reset_index(drop=True)
        cur = []
        for k in range(len(g)):
            if pd.isna(g.loc[k, "cycle_length"]):              # conception cycle: no length, ends the run
                if cur:
                    runs.append((int(w), cur)); cur = []
                continue
            if drop_missing_flag and str(g.loc[k, "conception_cycle"]).strip().lower() == "missing":
                if cur:
                    runs.append((int(w), cur)); cur = []
                continue
            if drop_last_record and k == len(g) - 1:
                continue
            gap = (g.loc[k, "start"] - g.loc[k - 1, "end"]).days - 1 if cur else 0
            cnj = int(g.loc[k, "cycle_number"] - g.loc[k - 1, "cycle_number"]) if cur else 1
            if cur and ((gap_tol is not None and gap > gap_tol) or (split_cn_jump and cnj != 1)):
                runs.append((int(w), cur)); cur = []
            if start_to_start and k + 1 < len(g):
                L = int((g.loc[k + 1, "start"] - g.loc[k, "start"]).days)
            else:
                L = int(g.loc[k, "cycle_length"])
            cur.append(L)
        if cur:
            runs.append((int(w), cur))
    return runs


def longest(runs, tie="earliest", minlen=5):
    best = {}
    for w, r in runs:                       # runs are in date order
        if len(r) < minlen:
            continue
        if w not in best or len(r) > len(best[w]) or (tie == "latest" and len(r) == len(best[w])):
            best[w] = r
    return [(w, best[w]) for w in sorted(best)]


def first_run(runs, minlen=5):
    best = {}
    for w, r in runs:
        if len(r) >= minlen and w not in best:
            best[w] = r
    return [(w, best[w]) for w in sorted(best)]


def main():
    d = load_raw()
    base_runs = segment()
    specs = {
        "baseline": baseline_runs(base_runs),
        "tie_latest": longest(base_runs, tie="latest"),
        "first_run_ge5": first_run(base_runs),
        "gap_tol_7d": longest(runs_with(d, gap_tol=7)),
        "gap_tol_45d": longest(runs_with(d, gap_tol=45)),
        "no_split_recorded_L": longest(runs_with(d, gap_tol=None)),
        "no_split_start_to_start": longest(runs_with(d, gap_tol=None, start_to_start=True)),
        "drop_missing_flag": longest(runs_with(d, drop_missing_flag=True)),
        "all5_run_min4": [(w, r) for w, r in base_runs if len(r) >= 4],
        "drop_last_record": longest(runs_with(d, drop_last_record=True)),
        "drop_last_record_no_split": longest(runs_with(d, gap_tol=None, drop_last_record=True)),
        "split_cn_jump_too": longest(runs_with(d, split_cn_jump=True)),
    }
    assert [L for _, r in specs["baseline"] for L in r] == \
        pd.read_csv(os.path.join(EXT, "data", "sequences.csv"))["LengthofCycle"].tolist()
    out = dict(B=B, seed=SEED, W_note="within-run permutation for every specification")
    for name, segs in specs.items():
        res, _ = analyse(name, segs)
        Ls = [L for _, r in segs for L in r]
        res["L_range"] = [int(min(Ls)), int(max(Ls))]
        res["n_L_outside_18_54"] = int(sum(1 for L in Ls if L < 18 or L > 54))
        out[name] = res
        with open(OUT, "w") as fh:
            json.dump(out, fh, indent=1)
    print(f"wrote {os.path.relpath(OUT, EXT)}")


if __name__ == "__main__":
    main()
