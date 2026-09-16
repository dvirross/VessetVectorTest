"""Derive every reported statistic from results/null_replicates.npz.

Prints a human-readable summary and writes results/summary.json. Thin wrapper around
``analysis_engine.summarise`` (verified to reproduce the previous summary.json exactly;
the engine adds a few extra keys such as lower-tail counts and the null covariance).
"""
import sys, os, json
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from analysis_engine import Dataset, summarise

ROOT = os.path.join(os.path.dirname(__file__), "..")
ds = Dataset.from_csv(os.path.join(ROOT, "data", "FilteredData.csv"), "fehring")
R = np.load(os.path.join(ROOT, "results", "null_replicates.npz"))
out = summarise(ds, R)
with open(os.path.join(ROOT, "results", "summary.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print("\nWrote results/summary.json")
