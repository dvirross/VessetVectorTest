"""Simulation validation of the testing procedures used in the paper.

Part A  Size (type-I error) of the marginal two-sided Monte Carlo tests, the
        Holm-corrected family, and the joint Mahalanobis test under each null
        model, using fresh pseudo-observed datasets drawn from the same null.
        The joint test is evaluated in two ways: (i) plug-in, with the null
        mean/covariance and the D_M reference distribution taken from the
        same replicate batch (as in the paper); (ii) split-batch, with the
        mean/covariance from one batch and the D_M reference from another.
Part B  Power of the within-woman (H_W) tests against two families of
        within-woman temporal-dependence alternatives, each simulated dataset
        analysed with its own B_IN within-woman permutations:
          AR(1)       each woman's sequence is a Gaussian AR(1) process with
                      her observed mean and SD, rounded to whole days;
          persistence each cycle repeats the previous haflaga exactly with
                      probability q, otherwise is drawn from the woman's own
                      observed multiset (q = 0 is exchangeable: size check).
Part C  Sensitivity of the H_W analysis to unequal follow-up: each woman is
        truncated to her first 12 recorded cycles.

Thin wrapper around ``analysis_engine.run_validation`` (same seeds and RNG order as the
original script; Parts A, B and C verified against the cached file).
Writes results/sim_validation.npz and prints a summary.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from analysis_engine import Dataset, run_validation, print_validation
from patterns import PATTERN_LABELS

ROOT = os.path.join(os.path.dirname(__file__), "..")

if __name__ == "__main__":
    ds = Dataset.from_csv(os.path.join(ROOT, "data", "FilteredData.csv"), "fehring")
    res = run_validation(ds, n_proc=3)
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    np.savez_compressed(os.path.join(ROOT, "results", "sim_validation.npz"),
                        part_a=np.array(res["part_a"], dtype=object), part_b=np.array(res["part_b"], dtype=object),
                        part_c=np.array([res["part_c"]], dtype=object), labels=np.array(PATTERN_LABELS),
                        B_REF=res["B_REF"], N_FRESH=res["N_FRESH"], R_ALT=res["R_ALT"], B_IN=res["B_IN"])
    print_validation(res)
