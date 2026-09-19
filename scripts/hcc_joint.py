"""Joint (Mahalanobis) test with Haflaga Chozer Chalila added as a sixth component.

Re-creates the three null replicate streams with the SAME seed and RNG consumption order
as rerun_nulls.py, so the first five columns reproduce results/null_replicates.npz
exactly, and counts HCC on the same permuted sequences. Then compares the 5-component
joint test (as in the paper) with the 6-component version under each null.

Usage: python scripts/hcc_joint.py [B] [seed]   -> results/hcc_joint.npz
"""
import os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from patterns import (load_data, weekly_anchor_set, Counter, perm_global, perm_within,
                      multinomial_iid, mahalanobis_test, PATTERN_LABELS)
from hcc_rule import HCCCounter

ROOT = os.path.join(os.path.dirname(__file__), "..")
B = int(sys.argv[1]) if len(sys.argv) > 1 else 50_000
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 17

H, wid, bounds = load_data(os.path.join(ROOT, "data", "FilteredData.csv"))
C5 = Counter(wid, weekly_anchor_set(H)); C6 = HCCCounter(wid, bounds)
obs = np.r_[C5.count(H), C6.count(H)]
vals, cnt = np.unique(H, return_counts=True); probs = cnt / cnt.sum()

rng = np.random.default_rng(SEED)
out = {k: np.zeros((B, 6), np.int16) for k in ("perm", "multi", "within")}
t0 = time.time()
for r in range(B):                       # identical call order to rerun_nulls.py
    x = perm_global(H, rng);                out["perm"][r, :5] = C5.count(x);   out["perm"][r, 5] = C6.count(x)
    x = multinomial_iid(H, rng, vals, probs); out["multi"][r, :5] = C5.count(x);  out["multi"][r, 5] = C6.count(x)
    x = perm_within(H, bounds, rng);        out["within"][r, :5] = C5.count(x); out["within"][r, 5] = C6.count(x)
    if (r + 1) % 10_000 == 0:
        print(f"{r+1}/{B}  {time.time()-t0:.0f}s", flush=True)

ref = np.load(os.path.join(ROOT, "results", "null_replicates.npz"))
for k in ("perm", "multi", "within"):
    paired = np.array_equal(out[k][:, :5], ref[k][:B])
    print(f"{k}: first five columns identical to results/null_replicates.npz -> {paired}")
    assert paired, "replicate streams not reproduced; comparison would be unpaired"

res = {}
print("\nobserved 6-vector:", obs.tolist(), "(", PATTERN_LABELS + ["HCC"], ")")
half = B // 2
for name, k in (("H_G", "perm"), ("H_iid", "multi"), ("H_W", "within")):
    a = out[k]
    DM5, p5, dm5, _ = mahalanobis_test(a[:, :5], obs[:5])
    DM6, p6, dm6, S6 = mahalanobis_test(a, obs)
    DM6s, p6s, _, _ = mahalanobis_test(a[half:], obs, ref=a[:half])
    n5 = int(np.count_nonzero(dm5 >= DM5)); n6 = int(np.count_nonzero(dm6 >= DM6))
    hcc_ge1 = a[:, 5] >= 1
    print(f"{name}: 5-pattern D_M={DM5:.3f} p={p5:.4f} ({n5} exceed) | "
          f"6-pattern D_M={DM6:.3f} p={p6:.4f} ({n6} exceed; split-batch D_M={DM6s:.3f} p={p6s:.4f}) | "
          f"HCC>=1 in {hcc_ge1.sum()} replicates, of which {np.count_nonzero(dm6[hcc_ge1] >= DM6)} exceed D_M^(6); "
          f"cond. number of 6x6 Sigma = {np.linalg.cond(S6):.0f}")
    res[k] = dict(DM5=DM5, p5=p5, n5=n5, DM6=DM6, p6=p6, n6=n6, DM6_split=DM6s, p6_split=p6s,
                  n_hcc_ge1=int(hcc_ge1.sum()), cond6=float(np.linalg.cond(S6)))
np.savez_compressed(os.path.join(ROOT, "results", "hcc_joint.npz"), observed=obs, seed=SEED, B=B,
                    **{f"{k}6": v for k, v in out.items()},
                    **{f"{k}_{s}": v for k, d in res.items() for s, v in d.items()})
print("saved results/hcc_joint.npz")
