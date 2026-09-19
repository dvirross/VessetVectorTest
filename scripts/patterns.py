"""Vesset pattern counting and null-model generators.

Single source of truth for the counting rules used in the paper. Two
implementations are provided: a reference loop version (mirrors the notebook
definitions) and a vectorised version used for the simulations. The test in
``verify_equivalence`` checks that they agree.

Counting convention: each function counts distinct maximal contiguous runs
(within a woman's sequence) that satisfy the rule; a run longer than the
minimum counts once; a resumed run after a discordant cycle counts again.
Runs never cross a woman boundary.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

PATTERN_LABELS = ["Haflaga", "Dilug", "Week", "Week-Dilug", "Dilug-in-Dilug"]
WEEK_DILUG_ANCHOR = 30


def load_data(path: str):
    """Return (interval values H_k = L_k + 1 (the halachic haflaga), woman-id integer array, boundaries).

    Also asserts that every woman's CycleNumber sequence is contiguous
    (1, 2, ..., n_i) so that no pattern window can bridge an excluded cycle.
    """
    data = pd.read_csv(path, low_memory=False)
    H = data["LengthofCycle"].to_numpy(dtype=np.int64) + 1
    ids = data["ClientID"].to_numpy()
    starts = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])
    bounds = list(zip(starts.tolist(), np.r_[starts[1:], len(H)].tolist()))
    wid = np.repeat(np.arange(len(bounds)), np.diff(np.r_[starts, len(H)]))
    for s, f in bounds:
        cn = data["CycleNumber"].to_numpy()[s:f].astype(int)
        if not np.array_equal(cn, np.arange(1, f - s + 1)):
            raise ValueError(f"Non-contiguous CycleNumber for woman index {s}")
    return H, wid, bounds


def weekly_anchor_set(H) -> set:
    """Haflaga values H = L + 1 whose cycle length L is a multiple of 7 (H = 7n + 1),
    i.e. successive onsets fall on the same weekday, restricted to the observed range."""
    mn, mx = int(H.min()), int(H.max())
    return {h for h in range(mn, mx + 1) if h % 7 == 1}


# ── Vectorised counter ─────────────────────────────────────────────────────────
class Counter:
    def __init__(self, wid: np.ndarray, weekly: set):
        self.wid = wid
        n = len(wid)
        self.same1 = np.zeros(n, bool)          # same woman as previous position
        self.same1[1:] = wid[1:] == wid[:-1]
        self.same2 = self.same1.copy()          # positions k-2, k-1, k all same woman
        self.same2[1:] &= self.same1[:-1]
        self.same3 = self.same2.copy()          # positions k-3..k all same woman
        self.same3[:2] = False
        self.same3[2:] &= self.same1[:-2]
        self.weekly = np.array(sorted(weekly), dtype=np.int64)

    @staticmethod
    def _n_runs(flag: np.ndarray) -> int:
        prev = np.zeros_like(flag)
        prev[1:] = flag[:-1]
        return int(np.count_nonzero(flag & ~prev))

    @staticmethod
    def _run_starts(flag: np.ndarray) -> np.ndarray:
        prev = np.zeros_like(flag)
        prev[1:] = flag[:-1]
        return flag & ~prev

    def run_starts(self, H: np.ndarray):
        """Boolean arrays marking the position at which each qualifying run is
        first detected, for the five patterns in PATTERN_LABELS order."""
        n = len(H)
        d = np.zeros(n, np.int64)
        d[1:] = H[1:] - H[:-1]
        eq = self.same1 & (d == 0)
        # Haflaga: >=3 consecutive equal
        triple = eq.copy()
        triple[1:] &= eq[:-1]
        # Dilug: >=3 cycles, constant non-zero first difference
        dil = self.same2 & (d != 0)
        dil[1:] &= d[1:] == d[:-1]
        dil[0] = False
        # Dilug-in-Dilug: >=4 cycles, constant non-zero second difference
        dd = np.zeros(n, np.int64)
        dd[1:] = d[1:] - d[:-1]
        did = self.same3 & (dd != 0)
        did[1:] &= dd[1:] == dd[:-1]
        did[0] = False
        return [self._run_starts(triple),
                self._run_starts(dil),
                self._run_starts(eq & np.isin(H, self.weekly)),       # Week
                self._run_starts(eq & (H == WEEK_DILUG_ANCHOR)),      # Week-Dilug
                self._run_starts(did)]

    def count(self, H: np.ndarray) -> np.ndarray:
        return np.array([int(np.count_nonzero(s)) for s in self.run_starts(H)])

    def count_by_woman(self, H: np.ndarray, n_women: int) -> np.ndarray:
        """(n_women, 5) matrix of run counts per woman."""
        return np.stack([np.bincount(self.wid[s], minlength=n_women)
                         for s in self.run_starts(H)], axis=1)


# ── Reference loop implementation (as in the notebook) ─────────────────────────
def count_all_loop(H, bounds, weekly) -> np.ndarray:
    def haflaga():
        n = 0
        for s, f in bounds:
            i = s + 2
            while i < f:
                if H[i] == H[i - 1] == H[i - 2]:
                    n += 1
                    while i + 1 < f and H[i + 1] == H[i]:
                        i += 1
                i += 1
        return n

    def dilug():
        n = 0
        for s, f in bounds:
            i = s + 2
            while i < f:
                d = H[i] - H[i - 1]
                if d != 0 and d == H[i - 1] - H[i - 2]:
                    n += 1
                    while i + 1 < f and H[i + 1] - H[i] == d:
                        i += 1
                i += 1
        return n

    def week():
        n = 0
        for s, f in bounds:
            i = s + 1
            while i < f:
                if H[i] in weekly and H[i] == H[i - 1]:
                    n += 1
                    while i + 1 < f and H[i + 1] == H[i]:
                        i += 1
                i += 1
        return n

    def week_dilug():
        n = 0
        for s, f in bounds:
            i = s + 1
            while i < f:
                if H[i] == WEEK_DILUG_ANCHOR == H[i - 1]:
                    n += 1
                    while i + 1 < f and H[i + 1] == H[i]:
                        i += 1
                i += 1
        return n

    def did():
        n = 0
        for s, f in bounds:
            i = s + 3
            while i < f:
                b1 = H[i - 2] - H[i - 3]
                b2 = H[i - 1] - H[i - 2]
                b3 = H[i] - H[i - 1]
                d1, d2 = b2 - b1, b3 - b2
                if d1 == d2 and d1 != 0:
                    n += 1
                    b, dd = b3, d2
                    while i + 1 < f and H[i + 1] == H[i] + b + dd:
                        i += 1
                        b += dd
                i += 1
        return n

    return np.array([haflaga(), dilug(), week(), week_dilug(), did()])


# ── Null-model generators ──────────────────────────────────────────────────────
def perm_global(H, rng):
    return rng.permutation(H)


def perm_within(H, bounds, rng):
    out = H.copy()
    for s, f in bounds:
        out[s:f] = rng.permutation(out[s:f])
    return out


def multinomial_iid(H, rng, vals=None, probs=None):
    if vals is None:
        vals, cnt = np.unique(H, return_counts=True)
        probs = cnt / cnt.sum()
    return rng.choice(vals, size=len(H), p=probs, replace=True)


# ── Inference helpers ──────────────────────────────────────────────────────────
def two_sided_p(null: np.ndarray, obs: float) -> float:
    B = len(null)
    up = (np.count_nonzero(null >= obs) + 1) / (B + 1)
    lo = (np.count_nonzero(null <= obs) + 1) / (B + 1)
    return min(1.0, 2 * min(up, lo))


def holm_adjust(p) -> np.ndarray:
    """Holm step-down adjusted p-values."""
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for k, i in enumerate(order):
        running = max(running, (m - k) * p[i])
        adj[i] = min(1.0, running)
    return adj


def mahalanobis_test(null: np.ndarray, obs: np.ndarray, ref: np.ndarray | None = None):
    """D_M of obs relative to null replicates.

    ``ref`` (default: ``null``) supplies the replicates used to estimate the
    null mean and covariance; the reference distribution of D_M is computed on
    ``null``. Passing disjoint batches gives the split-sample check.
    """
    if ref is None:
        ref = null
    mu = ref.mean(axis=0)
    Sigma = np.cov(ref.T)
    Sinv = np.linalg.inv(Sigma)
    diff = obs - mu
    DM = float(np.sqrt(diff @ Sinv @ diff))
    dn = null - mu
    dm_null = np.sqrt(np.einsum("ni,ij,nj->n", dn, Sinv, dn))
    p = (np.count_nonzero(dm_null >= DM) + 1) / (len(dm_null) + 1)
    return DM, p, dm_null, Sigma


def verify_equivalence(H, wid, bounds, n_perm=3000, seed=0):
    weekly = weekly_anchor_set(H)
    C = Counter(wid, weekly)
    rng = np.random.default_rng(seed)
    assert np.array_equal(C.count(H), count_all_loop(H, bounds, weekly))
    for _ in range(n_perm):
        for gen in (lambda: perm_global(H, rng), lambda: perm_within(H, bounds, rng),
                    lambda: multinomial_iid(H, rng)):
            x = gen()
            a, b = C.count(x), count_all_loop(x, bounds, weekly)
            if not np.array_equal(a, b):
                raise AssertionError(f"Mismatch: {a} vs {b}")
    return True


if __name__ == "__main__":
    import sys, time
    path = sys.argv[1] if len(sys.argv) > 1 else "data/FilteredData.csv"
    H, wid, bounds = load_data(path)
    print(f"{len(bounds)} women, {len(H)} cycles; CycleNumber contiguous for all women.")
    t0 = time.time()
    verify_equivalence(H, wid, bounds)
    print(f"Vectorised counter agrees with loop implementation on observed data and "
          f"9,000 random null sequences ({time.time()-t0:.1f}s).")
    C = Counter(wid, weekly_anchor_set(H))
    print("Observed counts:", dict(zip(PATTERN_LABELS, C.count(H).tolist())))
    t0 = time.time()
    rng = np.random.default_rng(1)
    for _ in range(2000):
        C.count(perm_within(H, bounds, rng))
    print(f"{(time.time()-t0)/2000*1e6:.0f} µs per within-woman replicate")
