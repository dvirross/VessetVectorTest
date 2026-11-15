"""Utah replication: sensitivity of the results to the run-selection rule (added 2026-11-15).

Responds to the 2026-11-01 preprocessing audit (branch audit/utah-preprocessing-2026-11-01),
which noted that the pre-specified rule keeps only each woman's single longest run of
>= 5 consecutive cycles, discarding 37 cycles in secondary eligible runs of 7 women.

The raw deposit is segmented exactly as in prepare_utah.py (duplicates dropped; a run ends at
a conception cycle or at any date gap; gaps are never bridged). Four specifications are
analysed with the engine's own null generators (seed 17, B = 50,000):

  baseline        longest run >= 5 per woman (the published analysis; must reproduce the
                  cached counts and sequences exactly)
  all5_run        every run >= 5; each run is a separate counting unit; H_W permutes within
                  each run
  all5_woman      every run >= 5; H_W pools each woman's values over all her runs and
                  permutes them into the run slots (run lengths and boundaries preserved;
                  events are still never counted across a gap)
  all3_run        every run >= 3 (the shortest window any pattern needs except Dilug-in-Dilug,
                  which needs 4); H_W within each run

H_G and H_iid always pool every analysed cycle. Output: results/all_runs_sensitivity.json.
"""
import json, os, sys, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from analysis_engine import Dataset, run_nulls  # noqa: E402
from patterns import (two_sided_p, holm_adjust, mahalanobis_test, perm_global,  # noqa: E402
                      PATTERN_LABELS)

RAW = os.path.join(EXT, "data", "raw", "CrMcyclelength_share.csv")
SEQ = os.path.join(EXT, "data", "sequences.csv")
OUT = os.path.join(EXT, "results", "all_runs_sensitivity.json")
B, SEED = 50_000, 17


def segment(raw_path=RAW):
    """List of (woman id, list of cycle lengths L) for every run, as in prepare_utah.py."""
    raw = pd.read_csv(raw_path)
    raw.columns = [c.lstrip("﻿") for c in raw.columns]
    d = raw.copy()
    d["start"] = pd.to_datetime(d["cycle_start_date"], format="%m/%d/%y")
    d["end"] = pd.to_datetime(d["cycle_end_date"], format="%m/%d/%y")
    d = d.drop_duplicates(["new_id", "start"]).sort_values(["new_id", "start"]).reset_index(drop=True)
    runs = []
    for w, g in d.groupby("new_id", sort=True):
        g = g.reset_index(drop=True)
        cur = []
        for k in range(len(g)):
            if pd.isna(g.loc[k, "cycle_length"]):
                if cur:
                    runs.append((int(w), cur)); cur = []
                continue
            if cur and (g.loc[k, "start"] - g.loc[k - 1, "end"]).days - 1 != 0:
                runs.append((int(w), cur)); cur = []
            cur.append(int(g.loc[k, "cycle_length"]))
        if cur:
            runs.append((int(w), cur))
    return runs


def baseline_runs(runs):
    best = {}
    for w, r in runs:                       # runs are in date order: '>' keeps the earliest on ties
        if len(r) >= 5 and (w not in best or len(r) > len(best[w])):
            best[w] = r
    return [(w, best[w]) for w in sorted(best)]


def pooled_within(ds, woman_of_run, rng, B):
    """H_W variant: permute each woman's values across all her runs (run slots fixed)."""
    idx = {}
    for i, w in enumerate(woman_of_run):
        s, f = ds.bounds[i]
        idx.setdefault(w, []).append(np.arange(s, f))
    groups = [np.concatenate(v) for v in idx.values()]
    out = np.zeros((B, 5), np.int16)
    for r in range(B):
        X = ds.H.copy()
        for ix in groups:
            X[ix] = ds.H[ix][rng.permutation(len(ix))]
        out[r] = ds.C.count(X)
    return out


def summarise(obs, null):
    mu, sd = null.mean(0), null.std(0, ddof=1)
    p = np.array([two_sided_p(null[:, j], obs[j]) for j in range(5)])
    DM, pDM, _, S = mahalanobis_test(null, obs)
    return dict(null_mean=mu.round(2).tolist(), null_sd=sd.round(2).tolist(),
                z=((obs - mu) / sd).round(2).tolist(), p_two=p.round(5).tolist(),
                holm=holm_adjust(p).round(4).tolist(), DM=round(DM, 3), p_DM=round(pDM, 5),
                cond=round(float(np.linalg.cond(S)), 2))


def analyse(name, segs, pooled=False, log=print):
    ds = Dataset.from_sequences([np.array(r) + 1 for _, r in segs], name)
    women = [w for w, _ in segs]
    t0 = time.time()
    R = run_nulls(ds, B=B, seed=SEED, store_by_woman=False, log=None)
    obs = R["observed"]
    res = dict(n_runs=len(segs), n_women=len(set(women)), n_cycles=int(ds.n_cycles),
               run_length_range=[int(min(len(r) for _, r in segs)), int(max(len(r) for _, r in segs))],
               observed=obs.tolist(), G=summarise(obs, R["perm"]), iid=summarise(obs, R["multi"]))
    if pooled:
        res["W"] = summarise(obs, pooled_within(ds, women, np.random.default_rng(SEED), B))
        res["W_note"] = "within-woman permutation pooled over each woman's runs"
    else:
        res["W"] = summarise(obs, R["within"])
        res["W_note"] = "within-run permutation"
    log(f"{name:12s} runs {res['n_runs']:3d} women {res['n_women']} cycles {res['n_cycles']} O={res['observed']} | "
        f"G D_M {res['G']['DM']:.2f} p {res['G']['p_DM']:.5f} | W D_M {res['W']['DM']:.2f} p {res['W']['p_DM']:.4f} "
        f"min Holm {min(res['W']['holm']):.3f} ({time.time()-t0:.0f}s)")
    return res, R


def main():
    runs = segment()
    base = baseline_runs(runs)
    # check: the baseline reproduces the released sequences file row for row
    seq = pd.read_csv(SEQ)
    rebuilt = [L for _, r in base for L in r]
    assert len(base) == seq["ClientID"].nunique() and rebuilt == seq["LengthofCycle"].tolist(), \
        "baseline segmentation does not reproduce data/sequences.csv"
    all5 = [(w, r) for w, r in runs if len(r) >= 5]
    all3 = [(w, r) for w, r in runs if len(r) >= 3]
    multi = sorted({w for w, _ in all5 if sum(1 for v, _ in all5 if v == w) > 1})
    out = dict(B=B, seed=SEED, labels=PATTERN_LABELS,
               segmentation=dict(n_runs_total=len(runs),
                                 women_with_more_than_one_run_ge5=len(multi),
                                 cycles_in_secondary_runs_ge5=int(sum(len(r) for _, r in all5) - sum(len(r) for _, r in base)),
                                 note="runs end at a conception cycle or a date gap, as in prepare_utah.py"))
    res, R = analyse("baseline", base)
    cached = np.load(os.path.join(EXT, "results", "null_replicates.npz"))
    res["reproduces_cached_replicates"] = bool(all(np.array_equal(R[k], cached[k]) for k in ("perm", "multi", "within")))
    out["baseline"] = res
    out["all5_run"] = analyse("all5_run", all5)[0]
    out["all5_woman"] = analyse("all5_woman", all5, pooled=True)[0]
    out["all3_run"] = analyse("all3_run", all3)[0]
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"baseline reproduces cached replicates: {res['reproduces_cached_replicates']}")
    print(f"wrote {os.path.relpath(OUT, EXT)}")


if __name__ == "__main__":
    main()
