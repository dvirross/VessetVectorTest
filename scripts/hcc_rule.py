"""Haflaga Chozer Chalila (HCC): the sixth rule computable from cycle lengths alone.

Definition (author, 2026-09-17): a block of m >= 3 consecutive haflaga values repeated
exactly, three times in a row (window 3m >= 9 cycles). (25,27,26,25,27,26,25,27,26) is an
event; (25,27,25,27,25,27) is not (block of length 2). Blocks may be longer than 3, e.g.
(25,27,26,28) x 3.

Counting convention (mirrors the five primary rules): distinct maximal m-periodic segments
of length >= 3m count once each; a segment whose block has a shorter period d | m (d < m)
is not counted at m (it is the same event at period d if d >= 3, a Haflaga if d = 1, and
not an HCC if d = 2). Rules that depend on calendar dates (e.g. Dilug Chozer Chalila) are
not computable from interval data and are not considered.

Usage: python scripts/hcc_rule.py [B] [seed]   -> results/hcc_null.npz
"""
import os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from patterns import load_data, perm_global, perm_within, multinomial_iid, two_sided_p

ROOT = os.path.join(os.path.dirname(__file__), "..")
MIN_BLOCK, REPS = 3, 3


def _minimal_period_is(block: np.ndarray, m: int) -> bool:
    for d in range(1, m):
        if m % d == 0 and np.array_equal(block[d:], block[:-d]):
            return False
    return True


def hcc_count_loop(H, bounds) -> int:
    """Reference implementation: per woman, per block length m, scan maximal runs."""
    n_events = 0
    for s, f in bounds:
        x = H[s:f]; n = len(x)
        for m in range(MIN_BLOCK, n // REPS + 1):
            k = m
            while k < n:
                if x[k] == x[k - m]:
                    j = k
                    while j + 1 < n and x[j + 1] == x[j + 1 - m]:
                        j += 1
                    # segment x[k-m : j+1], length j+1-(k-m)
                    if (j + 1 - (k - m)) >= REPS * m and _minimal_period_is(x[k - m:k], m):
                        n_events += 1
                    k = j + 1
                else:
                    k += 1
    return n_events


class HCCCounter:
    """Vectorised version: for each m, boolean 'equals value m back within the same woman',
    run detection with numpy, then a Python check on the (rare) long runs."""

    def __init__(self, wid: np.ndarray, bounds):
        self.wid = wid
        self.N = len(wid)
        self.m_max = max(f - s for s, f in bounds) // REPS
        self.same = {m: np.r_[np.zeros(m, bool), wid[m:] == wid[:-m]] for m in range(MIN_BLOCK, self.m_max + 1)}

    def count(self, H) -> int:
        n_events = 0
        N = self.N
        for m in range(MIN_BLOCK, self.m_max + 1):
            e = self.same[m].copy()
            e[m:] &= H[m:] == H[:-m]
            if not e.any():
                continue
            prev = np.r_[False, e[:-1]]
            nxt = np.r_[e[1:], False]
            starts = np.flatnonzero(e & ~prev)
            ends = np.flatnonzero(e & ~nxt)
            long = (ends - starts + 1) >= (REPS - 1) * m        # run of m-back equalities
            for k in starts[long]:
                if _minimal_period_is(H[k - m:k], m):
                    n_events += 1
        return n_events


def verify(H, wid, bounds, n_seq=3000, seed=0):
    C = HCCCounter(wid, bounds)
    rng = np.random.default_rng(seed)
    assert C.count(H) == hcc_count_loop(H, bounds)
    # sequences over a small alphabet (many periodic runs) plus planted HCC blocks
    for t in range(n_seq):
        x = rng.integers(25, 28, size=len(H))
        for s, f in bounds:
            if f - s >= 9 and rng.random() < 0.3:
                m = int(rng.integers(3, min(5, (f - s) // 3) + 1))
                blk = rng.integers(25, 32, size=m)
                x[s:s + REPS * m] = np.tile(blk, REPS)
        a, b = C.count(x), hcc_count_loop(x, bounds)
        assert a == b, (t, a, b)
    return True


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 50_000
    SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 17
    H, wid, bounds = load_data(os.path.join(ROOT, "data", "FilteredData.csv"))
    C = HCCCounter(wid, bounds)
    t0 = time.time(); verify(H, wid, bounds); print(f"equivalence OK ({time.time()-t0:.1f}s)")
    obs = C.count(H); print("observed HCC events:", obs)
    print("women with >= 9 cycles:", sum(f - s >= 9 for s, f in bounds), "of", len(bounds))
    # worked examples from the definition
    ex = lambda seq: hcc_count_loop(np.array(seq), [(0, len(seq))])
    print("examples:", ex([25,27,26,25,27,26,25,27,26]), ex([25,27,25,27,25,27]),
          ex([25,27,26,28,25,27,26,28,25,27,26,28]), ex([28]*9), ex([25,27]*6))
    vals, cnt = np.unique(H, return_counts=True); probs = cnt / cnt.sum()
    rng = np.random.default_rng(SEED)
    out = {k: np.zeros(B, np.int16) for k in ("perm", "multi", "within")}
    t0 = time.time()
    for r in range(B):
        out["perm"][r] = C.count(perm_global(H, rng))
        out["multi"][r] = C.count(multinomial_iid(H, rng, vals, probs))
        out["within"][r] = C.count(perm_within(H, bounds, rng))
        if (r + 1) % 10_000 == 0:
            print(f"{r+1}/{B}  {time.time()-t0:.0f}s", flush=True)
    np.savez_compressed(os.path.join(ROOT, "results", "hcc_null.npz"), observed=obs, seed=SEED, B=B, **out)
    for name, key in (("H_G", "perm"), ("H_iid", "multi"), ("H_W", "within")):
        a = out[key]
        print(f"{name}: mean {a.mean():.4f} sd {a.std(ddof=1):.4f} max {a.max()} "
              f"P(>=1) {np.mean(a >= 1):.4f} two-sided p(obs={obs}) {two_sided_p(a, obs):.3f}")
