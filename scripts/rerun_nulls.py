"""Re-run the three null models (H_G, H_iid, H_W) with B replicates.

Writes results/null_replicates.npz containing the observed count vector, the
(B, 5) replicate count matrices for each null, and for H_W the per-woman
(B, n_women, 5) count array used for the leave-one-woman-out analysis.

Usage: python scripts/rerun_nulls.py [B] [seed]
"""
import sys, time, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from patterns import (load_data, weekly_anchor_set, Counter, perm_global,
                      perm_within, multinomial_iid, PATTERN_LABELS)

B = int(sys.argv[1]) if len(sys.argv) > 1 else 50_000
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 17
ROOT = os.path.join(os.path.dirname(__file__), "..")

H, wid, bounds = load_data(os.path.join(ROOT, "data", "FilteredData.csv"))
n_w = len(bounds)
C = Counter(wid, weekly_anchor_set(H))
obs = C.count(H)
vals, cnt = np.unique(H, return_counts=True)
probs = cnt / cnt.sum()

rng = np.random.default_rng(SEED)
perm = np.zeros((B, 5), np.int16)
multi = np.zeros((B, 5), np.int16)
within = np.zeros((B, 5), np.int16)
within_by_woman = np.zeros((B, n_w, 5), np.int8)
t0 = time.time()
for r in range(B):
    perm[r] = C.count(perm_global(H, rng))
    multi[r] = C.count(multinomial_iid(H, rng, vals, probs))
    bw = C.count_by_woman(perm_within(H, bounds, rng), n_w)
    within_by_woman[r] = bw
    within[r] = bw.sum(axis=0)
    if (r + 1) % 10_000 == 0:
        print(f"{r+1}/{B}  {time.time()-t0:.0f}s", flush=True)

os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
np.savez_compressed(os.path.join(ROOT, "results", "null_replicates.npz"),
                    observed=obs, perm=perm, multi=multi, within=within,
                    within_by_woman=within_by_woman, seed=SEED, B=B,
                    n_cycles=len(H), n_women=n_w, labels=np.array(PATTERN_LABELS))
print("observed", dict(zip(PATTERN_LABELS, obs.tolist())))
for name, arr in [("H_G", perm), ("H_iid", multi), ("H_W", within)]:
    print(name, "mean", arr.mean(0).round(2), "sd", arr.std(0, ddof=1).round(2))
