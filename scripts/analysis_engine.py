"""Dataset-agnostic analysis engine for the vesset pattern study.

Every statistical step used for the paper (null replicates, marginal and joint tests,
leave-one-woman-out, pairwise exploratory analyses, simulation validation, truncation
sensitivity) is implemented here as a function of a ``Dataset`` object, so that the same
code runs on the Fehring sequences and on any external set of per-woman cycle-length
sequences. The counting rules themselves live in ``patterns.py``.

The original scripts (``rerun_nulls.py``, ``postprocess.py``, ``sim_validation.py``) are
thin wrappers around this module; their random-number streams are consumed in exactly the
same order as before, so the cached Fehring results are reproduced bit for bit.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from patterns import (load_data, weekly_anchor_set, Counter, perm_global, perm_within,  # noqa: E402
                      multinomial_iid, two_sided_p, holm_adjust, mahalanobis_test,
                      count_all_loop, PATTERN_LABELS)

SUBSPACE3 = [0, 3, 4]          # Haflaga, Week-Dilug, Dilug-in-Dilug (exploratory)


# ── Dataset container ──────────────────────────────────────────────────────────
@dataclass
class Dataset:
    """Per-woman contiguous sequences of haflaga values H = L + 1."""
    name: str
    H: np.ndarray
    wid: np.ndarray
    bounds: list
    weekly: set = field(init=False)
    C: Counter = field(init=False)

    def __post_init__(self):
        self.H = np.asarray(self.H, dtype=np.int64)
        self.wid = np.asarray(self.wid)
        self.weekly = weekly_anchor_set(self.H)
        self.C = Counter(self.wid, self.weekly)

    # constructors
    @classmethod
    def from_csv(cls, path, name):
        H, wid, bounds = load_data(path)
        return cls(name, H, wid, bounds)

    @classmethod
    def from_sequences(cls, seqs, name):
        seqs = [np.asarray(s, dtype=np.int64) for s in seqs if len(s) > 0]
        H = np.concatenate(seqs)
        lengths = np.array([len(s) for s in seqs])
        wid = np.repeat(np.arange(len(seqs)), lengths)
        starts = np.r_[0, np.cumsum(lengths)[:-1]]
        bounds = list(zip(starts.tolist(), (starts + lengths).tolist()))
        return cls(name, H, wid, bounds)

    # basic descriptors
    @property
    def n_women(self):
        return len(self.bounds)

    @property
    def n_cycles(self):
        return len(self.H)

    @property
    def n_i(self):
        return np.array([f - s for s, f in self.bounds])

    def sequences(self):
        return [self.H[s:f] for s, f in self.bounds]

    def observed(self):
        return self.C.count(self.H)

    def describe(self):
        L = self.H - 1
        n_i = self.n_i
        return dict(name=self.name, n_women=int(self.n_women), n_cycles=int(self.n_cycles),
                    cycles_per_woman_mean=round(float(n_i.mean()), 2),
                    cycles_per_woman_sd=round(float(n_i.std(ddof=1)), 2),
                    cycles_per_woman_min=int(n_i.min()), cycles_per_woman_max=int(n_i.max()),
                    cycles_per_woman_median=float(np.median(n_i)),
                    L_mean=round(float(L.mean()), 2), L_sd=round(float(L.std(ddof=1)), 2),
                    L_median=float(np.median(L)), L_min=int(L.min()), L_max=int(L.max()),
                    L_skew=round(float(stats.skew(L)), 2),
                    share_L_ge35=round(float((L >= 35).mean()), 3),
                    weekly_anchors=sorted(int(a) for a in self.weekly),
                    observed=dict(zip(PATTERN_LABELS, self.observed().tolist())))

    # derived datasets
    def truncated(self, max_cycles=12):
        return Dataset.from_sequences([s[:max_cycles] for s in self.sequences()],
                                      f"{self.name}_first{max_cycles}")

    def restricted_support(self, lo, hi, min_cycles=5, one_block_per_woman=False):
        """Sensitivity: keep only cycles with lo <= L <= hi. A woman's sequence is split
        at every out-of-range cycle (never bridged); each resulting run of at least
        ``min_cycles`` cycles becomes its own block. With ``one_block_per_woman`` only
        the longest run (earliest on ties) is kept. Returns (Dataset, info dict)."""
        blocks, n_women_kept, n_split_women, removed = [], 0, 0, 0
        for s in self.sequences():
            L = s - 1
            ok = (L >= lo) & (L <= hi)
            removed += int((~ok).sum())
            runs, cur = [], []
            for h, k in zip(s, ok):
                if k:
                    cur.append(h)
                else:
                    if cur:
                        runs.append(np.array(cur))
                    cur = []
            if cur:
                runs.append(np.array(cur))
            runs = [r for r in runs if len(r) >= min_cycles]
            if not runs:
                continue
            if one_block_per_woman:
                runs = [max(runs, key=len)]
            n_women_kept += 1
            n_split_women += int(len(runs) > 1)
            blocks.extend(runs)
        info = dict(lo=lo, hi=hi, min_cycles=min_cycles, cycles_removed=removed,
                    women_kept=n_women_kept, women_split_into_several_blocks=n_split_women,
                    n_blocks=len(blocks), n_cycles=int(sum(len(b) for b in blocks)))
        return Dataset.from_sequences(blocks, f"{self.name}_L{lo}-{hi}"), info

    def verify_counting(self, n_perm=3000, seed=0):
        """Vectorised counter == reference loop on observed data and 3*n_perm null sequences."""
        rng = np.random.default_rng(seed)
        assert np.array_equal(self.C.count(self.H), count_all_loop(self.H, self.bounds, self.weekly))
        for _ in range(n_perm):
            for x in (perm_global(self.H, rng), perm_within(self.H, self.bounds, rng),
                      multinomial_iid(self.H, rng)):
                a, b = self.C.count(x), count_all_loop(x, self.bounds, self.weekly)
                if not np.array_equal(a, b):
                    raise AssertionError(f"Mismatch: {a} vs {b}")
        return 3 * n_perm


# ── Null replicates ────────────────────────────────────────────────────────────
def run_nulls(ds: Dataset, B=50_000, seed=17, store_by_woman=True, log=print):
    """(B, 5) replicate count matrices under H_G, H_iid, H_W (one shared RNG stream,
    consumed in the order H_G, H_iid, H_W within each replicate) and, optionally, the
    per-woman H_W counts (B, n_women, 5) used for leave-one-woman-out."""
    import time
    n_w = ds.n_women
    vals, cnt = np.unique(ds.H, return_counts=True)
    probs = cnt / cnt.sum()
    rng = np.random.default_rng(seed)
    perm = np.zeros((B, 5), np.int16)
    multi = np.zeros((B, 5), np.int16)
    within = np.zeros((B, 5), np.int16)
    by_woman = np.zeros((B, n_w, 5), np.int8) if store_by_woman else None
    t0 = time.time()
    for r in range(B):
        perm[r] = ds.C.count(perm_global(ds.H, rng))
        multi[r] = ds.C.count(multinomial_iid(ds.H, rng, vals, probs))
        bw = ds.C.count_by_woman(perm_within(ds.H, ds.bounds, rng), n_w)
        if store_by_woman:
            by_woman[r] = bw
        within[r] = bw.sum(axis=0)
        if log and (r + 1) % 10_000 == 0:
            log(f"{r+1}/{B}  {time.time()-t0:.0f}s")
    out = dict(observed=ds.observed(), perm=perm, multi=multi, within=within, seed=seed, B=B,
               n_cycles=ds.n_cycles, n_women=n_w, labels=np.array(PATTERN_LABELS))
    if store_by_woman:
        out["within_by_woman"] = by_woman
    return out


# ── Marginal + joint summaries ─────────────────────────────────────────────────
def fmt_p(p):
    return "<.001" if p < .001 else f"{p:.3f}"


def summarise_null(arr, obs, B=None, label=None, log=None):
    """All marginal and joint statistics for one null's (B, 5) replicate matrix."""
    arr = np.asarray(arr, float)
    obs = np.asarray(obs, float)
    B = len(arr) if B is None else B
    mu, sd = arr.mean(0), arr.std(0, ddof=1)
    z = (obs - mu) / sd
    p = np.array([two_sided_p(arr[:, j], obs[j]) for j in range(5)])
    up = np.array([(np.count_nonzero(arr[:, j] >= obs[j]) + 1) / (B + 1) for j in range(5)])
    lo = np.array([(np.count_nonzero(arr[:, j] <= obs[j]) + 1) / (B + 1) for j in range(5)])
    n_ge = [int(np.count_nonzero(arr[:, j] >= obs[j])) for j in range(5)]
    n_le = [int(np.count_nonzero(arr[:, j] <= obs[j])) for j in range(5)]
    adj = holm_adjust(p)
    DM, pj, dm_null, Sig = mahalanobis_test(arr, obs)
    n_exc = int(np.count_nonzero(dm_null >= DM))
    half = B // 2
    DM_s, pj_s, _, _ = mahalanobis_test(arr[half:], obs, ref=arr[:half])
    DM3, pj3, dm3, _ = mahalanobis_test(arr[:, SUBSPACE3], obs[SUBSPACE3])
    corr = np.corrcoef(arr.T)
    corr_off = corr.copy(); np.fill_diagonal(corr_off, 0)
    out = dict(mean=mu.round(2).tolist(), sd=sd.round(2).tolist(), z=z.round(2).tolist(),
               p_two=p.round(4).tolist(), p_upper=up.round(5).tolist(), p_lower=lo.round(5).tolist(),
               n_ge_obs=n_ge, n_le_obs=n_le,
               holm=adj.round(4).tolist(), DM=round(DM, 3), p_DM=round(pj, 4), n_exceed_DM=n_exc,
               DM_split=round(DM_s, 3), p_DM_split=round(pj_s, 4),
               DM3=round(DM3, 3), p_DM3=round(pj3, 4), n_exceed_DM3=int(np.count_nonzero(dm3 >= DM3)),
               cond=round(float(np.linalg.cond(Sig)), 2), max_abs_corr=round(float(np.abs(corr_off).max()), 3),
               covariance=np.round(Sig, 3).tolist(), correlation=np.round(corr, 3).tolist(),
               pct95=np.percentile(arr, 95, axis=0).tolist(), pct975=np.percentile(arr, 97.5, axis=0).tolist())
    if log:
        log(f"\n── {label} ──")
        for j, l in enumerate(PATTERN_LABELS):
            log(f"{l:<15} obs {int(obs[j]):>3}  mean {mu[j]:6.2f}  sd {sd[j]:5.2f}  z {z[j]:+5.2f}  "
                f"p2 {p[j]:.4f} ({fmt_p(p[j])})  Holm {adj[j]:.4f}  #>=obs {n_ge[j]}")
        log(f"D_M = {DM:.3f}  p = {pj:.4f} ({n_exc} of {B} exceed)  | split-batch D_M = {DM_s:.3f} p = {pj_s:.4f}")
        log(f"3-pattern D_M = {DM3:.3f} p = {pj3:.4f} ({out['n_exceed_DM3']} exceed) | cond {out['cond']} | max|r| {out['max_abs_corr']}")
    return out


def run_details(ds: Dataset):
    """Haflaga runs at H = 30 and Dilug-in-Dilug runs with |second difference| = 1."""
    H, bounds = ds.H, ds.bounds
    haf30 = did1 = 0
    for s, f in bounds:
        x = H[s:f]; n = f - s
        i = 2
        while i < n:
            if x[i] == x[i-1] == x[i-2]:
                if x[i] == 30: haf30 += 1
                while i + 1 < n and x[i+1] == x[i]: i += 1
            i += 1
        i = 3
        while i < n:
            b1, b2, b3 = x[i-2]-x[i-3], x[i-1]-x[i-2], x[i]-x[i-1]
            d1, d2 = b2-b1, b3-b2
            if d1 == d2 and d1 != 0:
                if abs(d1) == 1: did1 += 1
                b, d = b3, d2
                while i + 1 < n and x[i+1] == x[i] + b + d: i += 1; b += d
            i += 1
    return dict(haflaga_runs_at_30=haf30, did_runs_abs_d_1=did1)


def sd_comparison(ds: Dataset, pw):
    """Within-woman SD of H: women with >=1 Haflaga run vs women with no pattern at all."""
    sd_i = np.array([ds.H[s:f].std(ddof=1) for s, f in ds.bounds])
    a, b_ = sd_i[pw[:, 0] > 0], sd_i[pw.sum(1) == 0]
    if len(a) < 2 or len(b_) < 2:
        return dict(n_haf=int(len(a)), n_none=int(len(b_)), note="too few women for a comparison")
    U, pU = stats.mannwhitneyu(a, b_, alternative="two-sided")
    t, pt = stats.ttest_ind(a, b_, equal_var=False)
    se = np.sqrt(a.var(ddof=1)/len(a) + b_.var(ddof=1)/len(b_))
    dfw = se**4 / ((a.var(ddof=1)/len(a))**2/(len(a)-1) + (b_.var(ddof=1)/len(b_))**2/(len(b_)-1))
    ci = (a.mean()-b_.mean()) + np.array([-1, 1]) * stats.t.ppf(.975, dfw) * se
    sp = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b_)-1)*b_.var(ddof=1)) / (len(a)+len(b_)-2))
    return dict(n_haf=int(len(a)), mean_sd_haf=round(a.mean(), 2), sd_between_haf=round(a.std(ddof=1), 2),
                n_none=int(len(b_)), mean_sd_none=round(b_.mean(), 2), sd_between_none=round(b_.std(ddof=1), 2),
                U=float(U), p_U=float(pU), t=round(float(t), 2), df=round(float(dfw), 1), p_t=float(pt),
                ci=ci.round(2).tolist(), cohen_d=round(float((a.mean()-b_.mean())/sp), 2))


def influence(ds: Dataset, pw, long_threshold=20):
    n_i = ds.n_i
    long = n_i >= long_threshold
    tot = pw.sum(0).astype(float)
    share = np.where(tot > 0, pw[long].sum(0) / np.where(tot > 0, tot, 1), np.nan)
    return dict(threshold=long_threshold, n_women_ge20=int(long.sum()),
                share_cycles_ge20=round(float(n_i[long].sum()/n_i.sum()), 3),
                share_events_ge20=np.round(share, 3).tolist(),
                spearman_n_vs_events=[round(float(stats.spearmanr(n_i, pw[:, j])[0]), 2) if pw[:, j].std() > 0 else None
                                      for j in range(5)])


def leave_one_out(obs, W, wbw, pw):
    """Exact LOO for H_W from additive per-woman replicate counts (no re-simulation)."""
    n_w = wbw.shape[1]
    loo_p, loo_pm = [], []
    for i in range(n_w):
        o_i = obs - pw[i]
        null_i = W - wbw[:, i, :]
        _, pj_i, _, _ = mahalanobis_test(null_i, o_i)
        loo_p.append(pj_i)
        loo_pm.append([two_sided_p(null_i[:, j], o_i[j]) for j in range(5)])
    loo_p = np.array(loo_p); loo_pm = np.array(loo_pm)
    return dict(joint_p_min=round(float(loo_p.min()), 4), joint_p_max=round(float(loo_p.max()), 4),
                n_joint_p_below_05=int((loo_p < .05).sum()),
                marginal_p_min=loo_pm.min(0).round(4).tolist(), marginal_p_max=loo_pm.max(0).round(4).tolist(),
                n_marg_below_05=(loo_pm < .05).sum(0).tolist(),
                most_influential_woman_index=int(np.argmax(np.abs(np.log(loo_p) - np.median(np.log(loo_p)))))), loo_p, loo_pm


def pairwise_rate_spearman(pw, n_i, Bp=50_000, seed=17):
    rates = pw / n_i[:, None]
    rng = np.random.default_rng(seed)
    pairs = []
    for (ja, la), (jb, lb) in combinations(enumerate(PATTERN_LABELS), 2):
        x, y = rates[:, ja], rates[:, jb]
        if x.std() == 0 or y.std() == 0:
            pairs.append(dict(pair=f"{la} x {lb}", r=None, p_perm=1.0)); continue
        r_obs = stats.spearmanr(x, y)[0]
        rx = stats.rankdata(x); ry = stats.rankdata(y)
        ryc = ry - ry.mean(); rxc = rx - rx.mean()
        denom = np.sqrt((rxc**2).sum() * (ryc**2).sum())
        perms = np.array([np.dot(rxc, rng.permutation(ryc)) / denom for _ in range(Bp)])
        p = (np.count_nonzero(np.abs(perms) >= abs(r_obs) - 1e-12) + 1) / (Bp + 1)
        pairs.append(dict(pair=f"{la} x {lb}", r=round(float(r_obs), 2), p_perm=round(float(p), 4)))
    ph = holm_adjust([q["p_perm"] for q in pairs])
    for q, h in zip(pairs, ph): q["p_holm"] = round(float(h), 3)
    pairs.sort(key=lambda q: q["p_perm"])
    return pairs


def pairwise_followup_robustness(pw, n_i, Bp=50_000, seed=17,
                                 pairs=((0, 4, "Haflaga x Dilug-in-Dilug"), (0, 3, "Haflaga x Week-Dilug"))):
    rates = pw / n_i[:, None]
    rng = np.random.default_rng(seed)
    strata = np.digitize(n_i, np.quantile(n_i, [.2, .4, .6, .8]))

    def strat_perm(y):
        out = y.copy()
        for k in np.unique(strata):
            idx = np.flatnonzero(strata == k); out[idx] = rng.permutation(y[idx])
        return out

    def partial_spearman(x, y, z):
        rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)
        ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz); ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
        return float(np.corrcoef(ex, ey)[0, 1])

    fu = {"spearman_n_vs_rate": [round(float(stats.spearmanr(n_i, rates[:, j])[0]), 2) if rates[:, j].std() > 0 else None
                                 for j in range(5)],
          "strata_sizes": np.bincount(strata).tolist(), "pairs": []}
    for ja, jb, lab in pairs:
        x, y = rates[:, ja], rates[:, jb]
        if x.std() == 0 or y.std() == 0:
            fu["pairs"].append(dict(pair=lab, note="degenerate")); continue
        r = float(stats.spearmanr(x, y)[0])
        perms = np.array([stats.spearmanr(x, strat_perm(y))[0] for _ in range(Bp)])
        p_s = (np.count_nonzero(np.abs(perms) >= abs(r) - 1e-12) + 1) / (Bp + 1)
        fu["pairs"].append(dict(pair=lab, r=round(r, 3), partial_r_given_n=round(partial_spearman(x, y, n_i), 3),
                                p_perm_stratified_by_followup_quintile=round(float(p_s), 4)))
    return fu


def pairwise_binary(pw, Bp=50_000, seed=17):
    flags = (pw > 0).astype(int)
    rng = np.random.default_rng(seed)
    binp = []
    for (ja, la), (jb, lb) in combinations(enumerate(PATTERN_LABELS), 2):
        x, y = flags[:, ja], flags[:, jb]
        ct = np.array([[((x==0)&(y==0)).sum(), ((x==0)&(y==1)).sum()], [((x==1)&(y==0)).sum(), ((x==1)&(y==1)).sum()]])
        if (ct.sum(0) == 0).any() or (ct.sum(1) == 0).any():
            binp.append(dict(pair=f"{la} x {lb}", chi2=None, p_perm=1.0, p_fisher=1.0)); continue
        chi2 = stats.chi2_contingency(ct, correction=False)[0]
        perms = np.empty(Bp)
        for r in range(Bp):
            yp = rng.permutation(y)
            ctp = np.array([[((x==0)&(yp==0)).sum(), ((x==0)&(yp==1)).sum()], [((x==1)&(yp==0)).sum(), ((x==1)&(yp==1)).sum()]])
            perms[r] = stats.chi2_contingency(ctp, correction=False)[0] if ctp.min() >= 0 else 0
        p = (np.count_nonzero(perms >= chi2 - 1e-12) + 1) / (Bp + 1)
        binp.append(dict(pair=f"{la} x {lb}", chi2=round(float(chi2), 2), p_perm=round(float(p), 4),
                         p_fisher=round(float(stats.fisher_exact(ct)[1]), 4)))
    ph = holm_adjust([q["p_perm"] for q in binp])
    for q, h in zip(binp, ph): q["p_holm"] = round(float(h), 3)
    binp.sort(key=lambda q: q["p_perm"])
    return binp


def summarise(ds: Dataset, R, Bp=50_000, pair_seed=17, log=print, pairwise=True):
    """Every reported statistic from a replicate archive ``R`` (dict or npz)."""
    obs = np.asarray(R["observed"], float)
    nulls = {"G": np.asarray(R["perm"], float), "iid": np.asarray(R["multi"], float), "W": np.asarray(R["within"], float)}
    B = int(R["B"])
    out = {"dataset": ds.name, "B": B, "seed": int(R["seed"]), "n_cycles": ds.n_cycles, "n_women": ds.n_women,
           "observed": obs.astype(int).tolist(), "description": ds.describe()}
    if log:
        log(f"B={B}, seed={out['seed']}, cycles={ds.n_cycles}, women={ds.n_women}")
        log("Observed: " + str(dict(zip(PATTERN_LABELS, obs.astype(int)))))
    for m, arr in nulls.items():
        out[m] = summarise_null(arr, obs, B, label=f"H_{m}", log=log)
    out["max_p_diff_G_iid"] = round(float(np.max(np.abs(np.array(out["G"]["p_two"]) - np.array(out["iid"]["p_two"])))), 3)
    out.update(run_details(ds))
    if log: log(f"\nHaflaga runs with H=30: {out['haflaga_runs_at_30']}; DiD runs with |d|=1: {out['did_runs_abs_d_1']}")
    pw = ds.C.count_by_woman(ds.H, ds.n_women)
    n_i = ds.n_i
    out["sd_comparison"] = sd_comparison(ds, pw)
    if log: log("SD comparison: " + str(out["sd_comparison"]))
    out["women_with_pattern"] = dict(zip(PATTERN_LABELS, (pw > 0).sum(0).tolist()))
    if log: log("Women establishing each pattern: " + str(out["women_with_pattern"]))
    out["influence"] = influence(ds, pw)
    if log: log("Influence: " + str(out["influence"]))
    if "within_by_woman" in (R.files if hasattr(R, "files") else R):
        wbw = np.asarray(R["within_by_woman"], np.int32)
        out["loo"], loo_p, loo_pm = leave_one_out(obs, nulls["W"], wbw, pw)
        out["_loo_joint_p"] = loo_p.round(5).tolist()
        if log: log("LOO: " + str(out["loo"]))
    if pairwise:
        out["pairwise_rate_spearman"] = pairwise_rate_spearman(pw, n_i, Bp, pair_seed)
        if log:
            log("\nPairwise rate-based Spearman:")
            for q in out["pairwise_rate_spearman"]:
                log(f"  {q['pair']:<32} r={q['r']}  p_perm={q['p_perm']:.4f}  Holm={q['p_holm']:.3f}")
        out["pairwise_followup_robustness"] = pairwise_followup_robustness(pw, n_i, Bp, pair_seed)
        if log: log("\nFollow-up robustness: " + str(out["pairwise_followup_robustness"]))
        out["pairwise_binary"] = pairwise_binary(pw, Bp, pair_seed)
        if log:
            log("\nPairwise binary (chi2 perm):")
            for q in out["pairwise_binary"][:4]:
                log(f"  {q['pair']:<32} chi2={q['chi2']} p_perm={q['p_perm']:.4f} Holm={q['p_holm']:.3f} Fisher={q['p_fisher']:.4f}")
    return out


# ── Simulation validation (size, power, truncation) ────────────────────────────
ALPHA = 0.05


def holm_reject(p, alpha=ALPHA):
    p = np.asarray(p); m = len(p); order = np.argsort(p)
    thresh = alpha / (m - np.arange(m))
    rej = np.zeros(m, bool)
    for k, idx in enumerate(order):
        if p[idx] <= thresh[k]:
            rej[idx] = True
        else:
            break
    return rej


class Validator:
    """Holds one dataset and reproduces sim_validation.py's Parts A, B and C."""

    def __init__(self, ds: Dataset, B_REF=20_000, N_FRESH=10_000, R_ALT=400, B_IN=2_000):
        self.ds = ds
        self.B_REF, self.N_FRESH, self.R_ALT, self.B_IN = B_REF, N_FRESH, R_ALT, B_IN
        vals, cnt = np.unique(ds.H, return_counts=True)
        self.vals, self.probs = vals, cnt / cnt.sum()
        self.LO, self.HI = int(ds.H.min()), int(ds.H.max())
        self.w_mean = np.array([ds.H[s:f].mean() for s, f in ds.bounds])
        self.w_sd = np.array([ds.H[s:f].std(ddof=1) for s, f in ds.bounds])

    def draw(self, null, rng):
        if null == "G":
            return perm_global(self.ds.H, rng)
        if null == "iid":
            return multinomial_iid(self.ds.H, rng, self.vals, self.probs)
        return perm_within(self.ds.H, self.ds.bounds, rng)

    def part_a(self, null, seed):
        C, B_REF, N_FRESH = self.ds.C, self.B_REF, self.N_FRESH
        rng = np.random.default_rng(seed)
        R1 = np.array([C.count(self.draw(null, rng)) for _ in range(B_REF)], float)
        R2 = np.array([C.count(self.draw(null, rng)) for _ in range(B_REF)], float)
        F = np.array([C.count(self.draw(null, rng)) for _ in range(N_FRESH)], float)
        up = (np.array([(R1[:, j][None, :] >= F[:, j][:, None]).sum(1) for j in range(5)]).T + 1) / (B_REF + 1)
        lo = (np.array([(R1[:, j][None, :] <= F[:, j][:, None]).sum(1) for j in range(5)]).T + 1) / (B_REF + 1)
        pm = np.minimum(1, 2 * np.minimum(up, lo))
        marg_rate = (pm <= ALPHA).mean(0)
        holm_fwer = np.mean([holm_reject(p).any() for p in pm])
        mu, Sig = R1.mean(0), np.cov(R1.T)
        Sinv = np.linalg.inv(Sig)
        dm = lambda X: np.sqrt(np.einsum("ni,ij,nj->n", X - mu, Sinv, X - mu))
        ref1, ref2, dF = np.sort(dm(R1)), np.sort(dm(R2)), dm(F)
        p_plug = (B_REF - np.searchsorted(ref1, dF, side="left") + 1) / (B_REF + 1)
        p_split = (B_REF - np.searchsorted(ref2, dF, side="left") + 1) / (B_REF + 1)
        return dict(null=null, marg_rate=marg_rate, holm_fwer=holm_fwer,
                    joint_plug_05=(p_plug <= .05).mean(), joint_plug_01=(p_plug <= .01).mean(),
                    joint_split_05=(p_split <= .05).mean(), joint_split_01=(p_split <= .01).mean(),
                    cond=np.linalg.cond(Sig), n_fresh=N_FRESH,
                    mc_se_05=float(np.sqrt(.05 * .95 / N_FRESH)))

    def gen_ar1(self, rho, rng):
        H, bounds = self.ds.H, self.ds.bounds
        out = np.empty(len(H), np.int64)
        for i, (s, f) in enumerate(bounds):
            n = f - s
            sd = max(self.w_sd[i], 0.5)
            e = rng.standard_normal(n) * sd
            x = np.empty(n)
            x[0] = e[0]
            c = np.sqrt(1 - rho ** 2)
            for k in range(1, n):
                x[k] = rho * x[k - 1] + c * e[k]
            out[s:f] = np.clip(np.rint(self.w_mean[i] + x), self.LO, self.HI)
        return out

    def gen_persist(self, q, rng):
        H, bounds = self.ds.H, self.ds.bounds
        out = np.empty(len(H), np.int64)
        for s, f in bounds:
            n = f - s
            pool = H[s:f]
            x = rng.choice(pool, size=n, replace=True)
            rep = rng.random(n) < q
            for k in range(1, n):
                if rep[k]:
                    x[k] = x[k - 1]
            out[s:f] = x
        return out

    def analyse_hw(self, X, rng):
        C, bounds = self.ds.C, self.ds.bounds
        obs = C.count(X).astype(float)
        null = np.array([C.count(perm_within(X, bounds, rng)) for _ in range(self.B_IN)], float)
        pm = np.array([two_sided_p(null[:, j], obs[j]) for j in range(5)])
        DM, pj, _, _ = mahalanobis_test(null, obs)
        return obs, null.mean(0), pm, pj

    def part_b(self, family, level, seed):
        rng = np.random.default_rng(seed)
        gen = (lambda: self.gen_ar1(level, rng)) if family == "ar1" else (lambda: self.gen_persist(level, rng))
        obs_all, mu_all, pm_all, pj_all = [], [], [], []
        for _ in range(self.R_ALT):
            obs, mu, pm, pj = self.analyse_hw(gen(), rng)
            obs_all.append(obs); mu_all.append(mu); pm_all.append(pm); pj_all.append(pj)
        pm_all = np.array(pm_all)
        return dict(family=family, level=level,
                    obs_mean=np.mean(obs_all, 0), null_mean=np.mean(mu_all, 0),
                    marg_power=(pm_all <= ALPHA).mean(0),
                    holm_power=np.mean([holm_reject(p).any() for p in pm_all]),
                    joint_power=(np.array(pj_all) <= ALPHA).mean(),
                    pj=np.array(pj_all), pm=pm_all)

    def part_c(self, seed, max_cycles=12, B=20_000):
        dt = self.ds.truncated(max_cycles)
        rng = np.random.default_rng(seed)
        obs = dt.C.count(dt.H).astype(float)
        null = np.array([dt.C.count(perm_within(dt.H, dt.bounds, rng)) for _ in range(B)], float)
        pm = np.array([two_sided_p(null[:, j], obs[j]) for j in range(5)])
        DM, pj, _, _ = mahalanobis_test(null, obs)
        return dict(n_cycles=dt.n_cycles, n_women=dt.n_women, obs=obs, mu=null.mean(0), sd=null.std(0, ddof=1),
                    z=(obs - null.mean(0)) / null.std(0, ddof=1), p=pm, holm=holm_adjust(pm), DM=DM, pj=pj)


# multiprocessing helpers: the Validator is installed once per worker
_VAL = None


def _init_worker(ds_args, kw):
    global _VAL
    _VAL = Validator(Dataset(*ds_args), **kw)


def _part_b_job(args):
    return _VAL.part_b(*args)


def run_validation(ds: Dataset, n_proc=3, log=print, B_REF=20_000, N_FRESH=10_000, R_ALT=400, B_IN=2_000,
                   ar1_levels=(0.0, 0.25, 0.5, 0.75), persist_levels=(0.0, 0.1, 0.2, 0.3),
                   seeds_a=(100, 101, 102), seed_c=300):
    """Parts A (size), B (H_W power) and C (first-12 truncation) with the paper's seeds."""
    import time
    from multiprocessing import Pool
    t0 = time.time()
    kw = dict(B_REF=B_REF, N_FRESH=N_FRESH, R_ALT=R_ALT, B_IN=B_IN)
    V = Validator(ds, **kw)
    jobs_b = [("ar1", r, 1000 + i) for i, r in enumerate(ar1_levels)] + \
             [("persist", q, 2000 + i) for i, q in enumerate(persist_levels)]
    with Pool(n_proc, initializer=_init_worker, initargs=((ds.name, ds.H, ds.wid, ds.bounds), kw)) as pool:
        res_b_async = pool.map_async(_part_b_job, jobs_b)
        res_a = [V.part_a(n, s) for n, s in zip(["G", "iid", "W"], seeds_a)]
        if log: log(f"Part A done {time.time()-t0:.0f}s")
        res_c = V.part_c(seed_c)
        if log: log(f"Part C done {time.time()-t0:.0f}s")
        res_b = res_b_async.get()
    if log: log(f"Part B done {time.time()-t0:.0f}s")
    return dict(part_a=res_a, part_b=res_b, part_c=res_c, B_REF=B_REF, N_FRESH=N_FRESH, R_ALT=R_ALT, B_IN=B_IN)


def print_validation(res, log=print):
    log("\n=== Part A: empirical size at alpha=.05 ===")
    for r in res["part_a"]:
        log(f"H_{r['null']:<4} marginal {np.round(r['marg_rate'],3)}  Holm-FWER {r['holm_fwer']:.3f}  "
            f"joint plug-in {r['joint_plug_05']:.3f} (.01: {r['joint_plug_01']:.3f})  "
            f"joint split {r['joint_split_05']:.3f} (.01: {r['joint_split_01']:.3f})  cond {r['cond']:.1f}")
    log("\n=== Part B: H_W power (alpha=.05) ===")
    for r in res["part_b"]:
        log(f"{r['family']} {r['level']:<5} obs mean {np.round(r['obs_mean'],1)} null mean {np.round(r['null_mean'],1)}"
            f"\n      marginal power {np.round(r['marg_power'],3)}  Holm-any {r['holm_power']:.3f}  joint {r['joint_power']:.3f}")
    c = res["part_c"]
    log("\n=== Part C: truncation to first 12 cycles ===")
    log(f"n cycles {c['n_cycles']}; obs {c['obs']}; mu {np.round(c['mu'],2)}; z {np.round(c['z'],2)}; "
        f"p {np.round(c['p'],3)}; D_M {c['DM']:.2f} p {c['pj']:.3f}")
