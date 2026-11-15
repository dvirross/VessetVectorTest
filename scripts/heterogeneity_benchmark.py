"""Heterogeneity-only benchmark (added 2026-11-15; SI Section S8).

Question: what does the three-null comparison return for a cohort in which women differ in
level and regularity but the ORDER of each woman's cycles carries no information?

For each cohort, R synthetic cohorts are generated: woman i keeps her n_i, and her intervals
are drawn i.i.d. from a normal distribution with her observed mean and SD (SD floored at
0.3 days), rounded to whole days and clipped to the cohort's observed range. Each synthetic
cohort is analysed with the counting code and test statistics of the main analysis, with B
replicates per null. Within-woman exchangeability holds by construction, so H_W is true;
H_G and H_iid are false whenever women differ.

Reported: distribution over synthetic cohorts of the H_G Haflaga z and joint D_M, the
H_G / H_W joint rejection rates at alpha = .05, and where the observed cohort falls.

Usage: python scripts/heterogeneity_benchmark.py [R B]   (default R = 500, B = 2,000)
-> results/heterogeneity_benchmark.json
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from patterns import load_data, Counter, weekly_anchor_set, perm_global, perm_within, mahalanobis_test

ROOT = os.path.join(os.path.dirname(__file__), "..")
COHORTS = {"Marquette": ("data/FilteredData.csv", "results/summary.json"),
           "Utah": ("external_replication/data/sequences.csv", "external_replication/results/summary.json")}
SEED = 17


def synthetic(H, bounds, rng):
    lo, hi = H.min(), H.max()
    out = [np.clip(np.rint(rng.normal(H[s:f].mean(), max(H[s:f].std(ddof=1), 0.3), f - s)), lo, hi)
           for s, f in bounds]
    return np.concatenate(out).astype(np.int64)


def one(H, wid, bounds, B, rng):
    C = Counter(wid, weekly_anchor_set(H))
    O = C.count(H)
    G = np.array([C.count(perm_global(H, rng)) for _ in range(B)])
    W = np.array([C.count(perm_within(H, bounds, rng)) for _ in range(B)])
    zG = (O[0] - G[:, 0].mean()) / G[:, 0].std(ddof=1)
    dG, pG, _, _ = mahalanobis_test(G, O)
    dW, pW, _, _ = mahalanobis_test(W, O)
    return zG, dG, pG, dW, pW


def main(R=500, B=2000):
    out = dict(R=R, B=B, seed=SEED, generator="per-woman i.i.d. rounded normal (own mean, SD >= 0.3), own n_i")
    for name, (path, summ) in COHORTS.items():
        H, wid, bounds = load_data(os.path.join(ROOT, path))
        rng = np.random.default_rng(SEED)
        res = np.array([one(synthetic(H, bounds, rng), wid, bounds, B, rng) for _ in range(R)])
        zG, dG, pG, dW, pW = res.T
        q = lambda x: [round(float(v), 2) for v in np.percentile(x, [2.5, 50, 97.5])]
        out[name] = dict(
            n_women=len(bounds), n_cycles=int(len(H)),
            haflaga_zG_q025_50_975=q(zG), DM_G_q025_50_975=q(dG), DM_W_q025_50_975=q(dW),
            reject_G_05=round(float(np.mean(pG <= .05)), 3), reject_W_05=round(float(np.mean(pW <= .05)), 3),
            min_haflaga_zG=round(float(zG.min()), 2))
        r = out[name]
        print(f"{name}: Haflaga z_G median {r['haflaga_zG_q025_50_975'][1]} (95% {r['haflaga_zG_q025_50_975'][0]}"
              f"-{r['haflaga_zG_q025_50_975'][2]}, min {r['min_haflaga_zG']}); D_M(G) median {r['DM_G_q025_50_975'][1]}; "
              f"reject H_G {r['reject_G_05']:.3f}, reject H_W {r['reject_W_05']:.3f}")
    with open(os.path.join(ROOT, "results", "heterogeneity_benchmark.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote results/heterogeneity_benchmark.json")


if __name__ == "__main__":
    main(*map(int, sys.argv[1:3])) if len(sys.argv) == 3 else main()
