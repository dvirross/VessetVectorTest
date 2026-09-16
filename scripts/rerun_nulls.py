"""Re-run the three null models (H_G, H_iid, H_W) with B replicates.

Writes results/null_replicates.npz containing the observed count vector, the
(B, 5) replicate count matrices for each null, and for H_W the per-woman
(B, n_women, 5) count array used for the leave-one-woman-out analysis.

Thin wrapper around ``analysis_engine.run_nulls`` (same RNG stream order as the
original script; output verified bit-identical to the cached file).

Usage: python scripts/rerun_nulls.py [B] [seed]
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from analysis_engine import Dataset, run_nulls
from patterns import PATTERN_LABELS

B = int(sys.argv[1]) if len(sys.argv) > 1 else 50_000
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 17
ROOT = os.path.join(os.path.dirname(__file__), "..")

ds = Dataset.from_csv(os.path.join(ROOT, "data", "FilteredData.csv"), "fehring")
R = run_nulls(ds, B, SEED)
os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
np.savez_compressed(os.path.join(ROOT, "results", "null_replicates.npz"), **R)
print("observed", dict(zip(PATTERN_LABELS, R["observed"].tolist())))
for name, arr in [("H_G", R["perm"]), ("H_iid", R["multi"]), ("H_W", R["within"])]:
    print(name, "mean", arr.mean(0).round(2), "sd", arr.std(0, ddof=1).round(2))
