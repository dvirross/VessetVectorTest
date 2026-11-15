"""Independent re-implementation of the counting rules, the H_G and H_W nulls and the joint
Mahalanobis test (added 2026-11-15; SI Section S8).

Deliberately imports nothing from this repository: the five rules are coded in plain Python
directly from their definitions in Section 2.1 of the paper (a maximal run of qualifying
positions counts once), the nulls use Python's own `random` module with a different seed,
and the Mahalanobis statistic is computed afresh. Agreement with the main pipeline guards
against an error shared by the main analysis and the replication driver, which use the same
counting module.

Usage: python -I scripts/independent_check.py [B]   (default B = 20,000)
-> results/independent_check.json
"""
import csv, json, math, os, random, sys
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
COHORTS = {"Marquette": "data/FilteredData.csv", "Utah": "external_replication/data/sequences.csv"}
SEED = 2026


def read(path):
    seqs = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            seqs.setdefault(r["ClientID"], []).append(int(r["LengthofCycle"]) + 1)   # interval H = L + 1
    return list(seqs.values())


def n_runs(flags):
    return sum(1 for i, f in enumerate(flags) if f and (i == 0 or not flags[i - 1]))


def make_counter(all_h):
    lo, hi = min(all_h), max(all_h)
    week = {h for h in range(lo, hi + 1) if h % 7 == 1}       # cycle length a multiple of 7
    week_dilug = {h for h in range(lo, hi + 1) if h % 7 == 2} # one weekday later each cycle

    def count(s):
        n = len(s)
        d = {i: s[i] - s[i - 1] for i in range(1, n)}
        haf = [i >= 2 and s[i] == s[i - 1] == s[i - 2] for i in range(n)]
        dil = [i >= 2 and d[i] != 0 and d[i] == d[i - 1] for i in range(n)]
        eq = [i >= 1 and s[i] == s[i - 1] for i in range(n)]
        wk = [eq[i] and s[i] in week for i in range(n)]
        wd = [eq[i] and s[i] in week_dilug for i in range(n)]
        dd = {i: d[i] - d[i - 1] for i in range(2, n)}
        did = [i >= 3 and dd[i] != 0 and dd[i] == dd[i - 1] for i in range(n)]
        return [n_runs(haf), n_runs(dil), n_runs(wk), n_runs(wd), n_runs(did)]
    return count


def analyse(seqs, B, rnd):
    all_h = [h for s in seqs for h in s]
    count = make_counter(all_h)
    total = lambda ss: np.sum([count(s) for s in ss], axis=0)
    O = total(seqs)
    lens = [len(s) for s in seqs]
    G, W = [], []
    for _ in range(B):
        W.append(total([rnd.sample(s, len(s)) for s in seqs]))
        pool = all_h[:]
        rnd.shuffle(pool)
        out, k = [], 0
        for L in lens:
            out.append(pool[k:k + L]); k += L
        G.append(total(out))
    res = dict(observed=O.tolist())
    for lab, N in (("G", np.array(G, float)), ("W", np.array(W, float))):
        mu = N.mean(0)
        Ci = np.linalg.inv(np.cov(N.T))
        dm = lambda v: math.sqrt(float((v - mu) @ Ci @ (v - mu)))
        DM = dm(O)
        null = np.array([dm(v) for v in N])
        res[lab] = dict(null_mean=mu.round(2).tolist(), z=((O - mu) / N.std(0, ddof=1)).round(2).tolist(),
                        DM=round(DM, 3), p_DM=round(float((np.sum(null >= DM) + 1) / (B + 1)), 5))
    return res


def main(B=20_000):
    rnd = random.Random(SEED)
    out = dict(B=B, seed=SEED, note="no code shared with scripts/patterns.py")
    for name, path in COHORTS.items():
        out[name] = r = analyse(read(os.path.join(ROOT, path)), B, rnd)
        print(f"{name}: observed {r['observed']} | H_G D_M {r['G']['DM']:.2f} (p {r['G']['p_DM']:.5f}) "
              f"| H_W D_M {r['W']['DM']:.2f} (p {r['W']['p_DM']:.3f}) | H_W z {r['W']['z']}")
    with open(os.path.join(ROOT, "results", "independent_check.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote results/independent_check.json")


if __name__ == "__main__":
    main(int(sys.argv[1])) if len(sys.argv) > 1 else main()
