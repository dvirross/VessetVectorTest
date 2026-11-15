"""Two exploratory complements to the primary H_W analysis (added 2026-11-15).

A. Directional statistics (from the cached replicates; no new simulation).
   The joint Mahalanobis statistic D_M is omnidirectional: it spreads its sensitivity over
   all five directions and cannot distinguish a coherent shift of every count in the same
   direction from scattered departures. Two one-dimensional summaries of the direction of
   the count vector are computed under each null:
     T_z   = sum_j (S_j - mu_j) / sigma_j   (equal weight per pattern; mu, sigma from the
             same replicates, as for D_M)
     T_tot = sum_j S_j                       (total rule-defined events)
   Two-sided Monte Carlo p-values use the convention of Eq. (1) of the paper. These
   statistics were chosen after the observed counts were seen to lie below their H_W means
   and are therefore exploratory.

B. Generic ordering statistics under H_W (positive control; fresh within-woman permutations).
   To show what within-woman ordering structure the data contain irrespective of the
   vesset rules, three standard serial statistics are tested against the same within-woman
   permutation null:
     r1     pooled within-woman lag-1 autocorrelation, sum x_k x_{k-1} / sum x_k^2 with x
            centred at each woman's own mean (the denominator and the centring are fixed
            under H_W because each woman's multiset is fixed)
     mad    mean absolute difference between consecutive intervals, mean |H_k - H_{k-1}|
     p_eq   proportion of consecutive pairs with H_k = H_{k-1} (the elementary event behind
            Haflaga, Week and Week-Dilug)
   Under H_W the expected r1 is negative (about -1/(n_i - 1) per woman), so the observed r1
   is compared with its null distribution, not with zero. Holm correction across the three.

Usage: python scripts/directional_ordering.py [DATA_CSV NULL_NPZ OUT_JSON]
  default: data/FilteredData.csv results/null_replicates.npz results/directional_ordering.json
  Utah:    external_replication/data/sequences.csv external_replication/results/null_replicates.npz
           results/directional_ordering_utah.json
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from patterns import load_data, weekly_anchor_set, Counter, two_sided_p, holm_adjust

ROOT = os.path.join(os.path.dirname(__file__), "..")
B_ORDER = 50_000
SEED_ORDER = 17
NULL_KEYS = {"G": "perm", "iid": "multi", "W": "within"}


def directional(rep):
    obs = rep["observed"].astype(float)
    out = {}
    for lab, key in NULL_KEYS.items():
        S = rep[key].astype(float)
        mu, sd = S.mean(0), S.std(0, ddof=1)
        tz_null = ((S - mu) / sd).sum(1)
        tz_obs = float(((obs - mu) / sd).sum())
        tt_null = S.sum(1)
        tt_obs = float(obs.sum())
        out[lab] = dict(
            T_z=round(tz_obs, 3), T_z_null_sd=round(float(tz_null.std(ddof=1)), 3),
            T_z_p_two=round(two_sided_p(tz_null, tz_obs), 5),
            T_z_p_lower=round((np.count_nonzero(tz_null <= tz_obs) + 1) / (len(S) + 1), 5),
            T_tot=int(tt_obs), T_tot_null_mean=round(float(tt_null.mean()), 2),
            T_tot_null_sd=round(float(tt_null.std(ddof=1)), 2),
            T_tot_z=round(float((tt_obs - tt_null.mean()) / tt_null.std(ddof=1)), 2),
            T_tot_p_two=round(two_sided_p(tt_null, tt_obs), 5),
            n_patterns_below_mean=int(np.sum(obs < mu)),
            P_all_five_at_or_below_mean=round(float(np.mean((S <= mu).all(1))), 4),
        )
    return out


class OrderingStats:
    """Vectorised r1, mad and p_eq for sequences stored contiguously by woman."""

    def __init__(self, H, wid):
        self.wid = wid
        self.same = np.r_[False, wid[1:] == wid[:-1]]            # position k has a predecessor in-woman
        n_w = wid.max() + 1
        means = np.bincount(wid, weights=H, minlength=n_w) / np.bincount(wid, minlength=n_w)
        self.mean_at = means[wid]

    def __call__(self, X):
        x = X - self.mean_at                                       # woman means are permutation-invariant
        s = self.same[1:]
        r1 = np.sum((x[1:] * x[:-1])[s]) / np.sum(x * x)
        d = (X[1:] - X[:-1])[s]
        return np.array([r1, np.mean(np.abs(d)), np.mean(d == 0)])


def within_perm_index(wid, rng):
    """Uniform within-woman permutation of positions (data stored contiguously by woman)."""
    return np.lexsort((rng.random(len(wid)), wid))


def ordering(H, wid, B=B_ORDER, seed=SEED_ORDER):
    f = OrderingStats(H.astype(float), wid)
    obs = f(H.astype(float))
    rng = np.random.default_rng(seed)
    Hf = H.astype(float)
    null = np.empty((B, 3))
    for r in range(B):
        null[r] = f(Hf[within_perm_index(wid, rng)])
    labels = ["lag1_autocorrelation", "mean_abs_consecutive_difference", "prop_equal_consecutive"]
    p = np.array([two_sided_p(null[:, k], obs[k]) for k in range(3)])
    holm = holm_adjust(p)
    return {lab: dict(observed=round(float(obs[k]), 4), null_mean=round(float(null[:, k].mean()), 4),
                      null_sd=round(float(null[:, k].std(ddof=1)), 4),
                      z=round(float((obs[k] - null[:, k].mean()) / null[:, k].std(ddof=1)), 2),
                      p_two=round(float(p[k]), 5), holm=round(float(holm[k]), 5))
            for k, lab in enumerate(labels)} | {"B": B, "seed": seed}


def run(DATA, NPZ, OUT=None, B=B_ORDER):
    H, wid, bounds = load_data(DATA)
    rep = np.load(NPZ)
    C = Counter(wid, weekly_anchor_set(H))
    assert np.array_equal(C.count(H), rep["observed"]), "observed counts differ from cached replicates file"
    res = dict(data=os.path.relpath(DATA, ROOT), n_cycles=int(len(H)), n_women=len(bounds),
               B_counts=int(rep[NULL_KEYS["W"]].shape[0]),
               directional=directional(rep), ordering_HW=ordering(H, wid, B=B))
    # robustness of B: women whose every interval corresponds to 18 <= L <= 54 (harmonised support)
    keep = np.array([H[s:f].min() >= 19 and H[s:f].max() <= 55 for s, f in bounds])
    if not keep.all():
        mask = keep[wid]
        wid_k = np.unique(wid[mask], return_inverse=True)[1]
        res["ordering_HW_women_all_L18_54"] = ordering(H[mask], wid_k, B=B) | {
            "n_women": int(keep.sum()), "n_cycles": int(mask.sum())}
    print(f"== {res['data']}: {res['n_women']} women, {res['n_cycles']} cycles")
    for lab, d in res["directional"].items():
        print(f"  H_{lab:3s} T_z = {d['T_z']:+.2f} (p2 {d['T_z_p_two']:.4f}); total {d['T_tot']} vs "
              f"{d['T_tot_null_mean']:.1f} (z {d['T_tot_z']:+.2f}, p2 {d['T_tot_p_two']:.4f}); "
              f"below mean {d['n_patterns_below_mean']}/5; P(all <= mean) {d['P_all_five_at_or_below_mean']:.3f}")
    for key in ("ordering_HW", "ordering_HW_women_all_L18_54"):
        if key not in res:
            continue
        print(f"  [{key}]")
        for lab, d in res[key].items():
            if isinstance(d, dict):
                print(f"  H_W {lab:33s} obs {d['observed']:.4f} null {d['null_mean']:.4f} ({d['null_sd']:.4f}) "
                      f"z {d['z']:+.2f} p2 {d['p_two']:.4f} Holm {d['holm']:.4f}")
    if OUT:
        with open(OUT, "w") as fh:
            json.dump(res, fh, indent=1)
        print(f"  wrote {os.path.relpath(OUT, ROOT)}")
    return res


if __name__ == "__main__":
    if len(sys.argv) == 4:
        run(*sys.argv[1:4])
    else:
        run(os.path.join(ROOT, "data", "FilteredData.csv"), os.path.join(ROOT, "results", "null_replicates.npz"),
            os.path.join(ROOT, "results", "directional_ordering.json"))
        run(os.path.join(ROOT, "external_replication", "data", "sequences.csv"),
            os.path.join(ROOT, "external_replication", "results", "null_replicates.npz"),
            os.path.join(ROOT, "results", "directional_ordering_utah.json"))
