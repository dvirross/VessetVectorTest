"""Derive every reported statistic from results/null_replicates.npz.

Prints a human-readable summary and writes results/summary.json.
"""
import sys, os, json
import numpy as np
from scipy import stats
from itertools import combinations

sys.path.insert(0, os.path.dirname(__file__))
from patterns import (load_data, weekly_anchor_set, Counter, two_sided_p,
                      mahalanobis_test, PATTERN_LABELS)

ROOT = os.path.join(os.path.dirname(__file__), "..")
H, wid, bounds = load_data(os.path.join(ROOT, "data", "FilteredData.csv"))
n_w = len(bounds)
C = Counter(wid, weekly_anchor_set(H))
R = np.load(os.path.join(ROOT, "results", "null_replicates.npz"))
obs = R["observed"].astype(float)
nulls = {"G": R["perm"].astype(float), "iid": R["multi"].astype(float), "W": R["within"].astype(float)}
wbw = R["within_by_woman"].astype(np.int32)
B = int(R["B"])
out = {"B": B, "seed": int(R["seed"]), "n_cycles": int(R["n_cycles"]), "n_women": n_w,
       "observed": obs.astype(int).tolist()}


def holm(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for k, i in enumerate(o):
        running = max(running, (m - k) * p[i])
        adj[i] = min(1.0, running)
    return adj


def fmt_p(p):
    return "<.001" if p < .001 else f"{p:.3f}"


print(f"B={B}, seed={out['seed']}, cycles={out['n_cycles']}, women={n_w}")
print("Observed:", dict(zip(PATTERN_LABELS, obs.astype(int))))

for m, arr in nulls.items():
    mu, sd = arr.mean(0), arr.std(0, ddof=1)
    z = (obs - mu) / sd
    p = np.array([two_sided_p(arr[:, j], obs[j]) for j in range(5)])
    up = np.array([(np.count_nonzero(arr[:, j] >= obs[j]) + 1) / (B + 1) for j in range(5)])
    n_ge = [int(np.count_nonzero(arr[:, j] >= obs[j])) for j in range(5)]
    adj = holm(p)
    DM, pj, dm_null, Sig = mahalanobis_test(arr, obs)
    n_exc = int(np.count_nonzero(dm_null >= DM))
    half = B // 2
    DM_s, pj_s, _, _ = mahalanobis_test(arr[half:], obs, ref=arr[:half])
    DM3, pj3, dm3, _ = mahalanobis_test(arr[:, [0, 3, 4]], obs[[0, 3, 4]])
    corr = np.corrcoef(arr.T); np.fill_diagonal(corr, 0)
    out[m] = dict(mean=mu.round(2).tolist(), sd=sd.round(2).tolist(), z=z.round(2).tolist(),
                  p_two=p.round(4).tolist(), p_upper=up.round(5).tolist(), n_ge_obs=n_ge,
                  holm=adj.round(4).tolist(), DM=round(DM, 3), p_DM=round(pj, 4), n_exceed_DM=n_exc,
                  DM_split=round(DM_s, 3), p_DM_split=round(pj_s, 4),
                  DM3=round(DM3, 3), p_DM3=round(pj3, 4), n_exceed_DM3=int(np.count_nonzero(dm3 >= DM3)),
                  cond=round(float(np.linalg.cond(Sig)), 2), max_abs_corr=round(float(np.abs(corr).max()), 3),
                  pct95=np.percentile(arr, 95, axis=0).tolist(), pct975=np.percentile(arr, 97.5, axis=0).tolist())
    print(f"\n── H_{m} ──")
    for j, l in enumerate(PATTERN_LABELS):
        print(f"{l:<15} obs {int(obs[j]):>3}  mean {mu[j]:6.2f}  sd {sd[j]:5.2f}  z {z[j]:+5.2f}  "
              f"p2 {p[j]:.4f} ({fmt_p(p[j])})  Holm {adj[j]:.4f}  #>=obs {n_ge[j]}")
    print(f"D_M = {DM:.3f}  p = {pj:.4f} ({n_exc} of {B} exceed)  | split-batch D_M = {DM_s:.3f} p = {pj_s:.4f}")
    print(f"3-pattern D_M = {DM3:.3f} p = {pj3:.4f} ({out[m]['n_exceed_DM3']} exceed) | cond {out[m]['cond']} | max|r| {out[m]['max_abs_corr']}")

# max p-value difference between H_G and H_iid
out["max_p_diff_G_iid"] = round(float(np.max(np.abs(np.array(out["G"]["p_two"]) - np.array(out["iid"]["p_two"])))), 3)

# ── run details: Haflaga value 30, DiD |d| = 1 ────────────────────────────────
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
out["haflaga_runs_at_30"] = haf30; out["did_runs_abs_d_1"] = did1
print(f"\nHaflaga runs with H=30: {haf30}; DiD runs with |d|=1: {did1}")

# ── per-woman counts, SD comparison ───────────────────────────────────────────
pw = C.count_by_woman(H, n_w)
n_i = np.array([f - s for s, f in bounds])
sd_i = np.array([H[s:f].std(ddof=1) for s, f in bounds])
haf_w = pw[:, 0] > 0
none_w = pw.sum(1) == 0
a, b_ = sd_i[haf_w], sd_i[none_w]
U, pU = stats.mannwhitneyu(a, b_, alternative="two-sided")
t, pt = stats.ttest_ind(a, b_, equal_var=False)
se = np.sqrt(a.var(ddof=1)/len(a) + b_.var(ddof=1)/len(b_))
dfw = se**4 / ((a.var(ddof=1)/len(a))**2/(len(a)-1) + (b_.var(ddof=1)/len(b_))**2/(len(b_)-1))
ci = (a.mean()-b_.mean()) + np.array([-1, 1]) * stats.t.ppf(.975, dfw) * se
sp = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b_)-1)*b_.var(ddof=1)) / (len(a)+len(b_)-2))
out["sd_comparison"] = dict(n_haf=int(len(a)), mean_sd_haf=round(a.mean(), 2), sd_between_haf=round(a.std(ddof=1), 2),
                            n_none=int(len(b_)), mean_sd_none=round(b_.mean(), 2), sd_between_none=round(b_.std(ddof=1), 2),
                            U=float(U), p_U=float(pU), t=round(float(t), 2), df=round(float(dfw), 1), p_t=float(pt),
                            ci=ci.round(2).tolist(), cohen_d=round(float((a.mean()-b_.mean())/sp), 2))
print("SD comparison:", out["sd_comparison"])
out["women_with_pattern"] = dict(zip(PATTERN_LABELS, (pw > 0).sum(0).tolist()))
print("Women establishing each pattern:", out["women_with_pattern"])

# ── influence of follow-up length ─────────────────────────────────────────────
long = n_i >= 20
out["influence"] = dict(n_women_ge20=int(long.sum()), share_cycles_ge20=round(float(n_i[long].sum()/n_i.sum()), 3),
                        share_events_ge20=(pw[long].sum(0)/pw.sum(0)).round(3).tolist(),
                        spearman_n_vs_events=[round(float(stats.spearmanr(n_i, pw[:, j])[0]), 2) for j in range(5)])
print("Influence:", out["influence"])

# ── leave-one-woman-out for H_W ───────────────────────────────────────────────
loo_p, loo_pm = [], []
W = nulls["W"]
for i in range(n_w):
    o_i = obs - pw[i]
    null_i = W - wbw[:, i, :]
    _, pj_i, _, _ = mahalanobis_test(null_i, o_i)
    loo_p.append(pj_i)
    loo_pm.append([two_sided_p(null_i[:, j], o_i[j]) for j in range(5)])
loo_p = np.array(loo_p); loo_pm = np.array(loo_pm)
out["loo"] = dict(joint_p_min=round(float(loo_p.min()), 4), joint_p_max=round(float(loo_p.max()), 4),
                  n_joint_p_below_05=int((loo_p < .05).sum()),
                  marginal_p_min=loo_pm.min(0).round(4).tolist(), marginal_p_max=loo_pm.max(0).round(4).tolist(),
                  n_marg_below_05=(loo_pm < .05).sum(0).tolist())
print("LOO:", out["loo"])

# ── pairwise rate-based Spearman with permutation p (B=50,000, seed 17) ───────
rates = pw / n_i[:, None]
rng = np.random.default_rng(17)
Bp = 50_000
pairs = []
for (ja, la), (jb, lb) in combinations(enumerate(PATTERN_LABELS), 2):
    x, y = rates[:, ja], rates[:, jb]
    r_obs = stats.spearmanr(x, y)[0]
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    # Spearman = Pearson on ranks; permute one rank vector
    ryc = ry - ry.mean(); rxc = rx - rx.mean()
    denom = np.sqrt((rxc**2).sum() * (ryc**2).sum())
    perms = np.array([np.dot(rxc, rng.permutation(ryc)) / denom for _ in range(Bp)])
    p = (np.count_nonzero(np.abs(perms) >= abs(r_obs) - 1e-12) + 1) / (Bp + 1)
    pairs.append(dict(pair=f"{la} x {lb}", r=round(float(r_obs), 2), p_perm=round(float(p), 4)))
ph = holm([q["p_perm"] for q in pairs])
for q, h in zip(pairs, ph): q["p_holm"] = round(float(h), 3)
pairs.sort(key=lambda q: q["p_perm"])
out["pairwise_rate_spearman"] = pairs
print("\nPairwise rate-based Spearman:")
for q in pairs: print(f"  {q['pair']:<32} r={q['r']:+.2f}  p_perm={q['p_perm']:.4f}  Holm={q['p_holm']:.3f}")

# binary robustness: Fisher exact + permutation of binary flags
flags = (pw > 0).astype(int)
binp = []
rng = np.random.default_rng(17)
for (ja, la), (jb, lb) in combinations(enumerate(PATTERN_LABELS), 2):
    x, y = flags[:, ja], flags[:, jb]
    ct = np.array([[((x==0)&(y==0)).sum(), ((x==0)&(y==1)).sum()], [((x==1)&(y==0)).sum(), ((x==1)&(y==1)).sum()]])
    chi2 = stats.chi2_contingency(ct, correction=False)[0]
    perms = np.empty(Bp)
    for r in range(Bp):
        yp = rng.permutation(y)
        ctp = np.array([[((x==0)&(yp==0)).sum(), ((x==0)&(yp==1)).sum()], [((x==1)&(yp==0)).sum(), ((x==1)&(yp==1)).sum()]])
        perms[r] = stats.chi2_contingency(ctp, correction=False)[0] if ctp.min() >= 0 else 0
    p = (np.count_nonzero(perms >= chi2 - 1e-12) + 1) / (Bp + 1)
    binp.append(dict(pair=f"{la} x {lb}", chi2=round(float(chi2), 2), p_perm=round(float(p), 4), p_fisher=round(float(stats.fisher_exact(ct)[1]), 4)))
ph = holm([q["p_perm"] for q in binp])
for q, h in zip(binp, ph): q["p_holm"] = round(float(h), 3)
binp.sort(key=lambda q: q["p_perm"])
out["pairwise_binary"] = binp
print("\nPairwise binary (chi2 perm):")
for q in binp[:4]: print(f"  {q['pair']:<32} chi2={q['chi2']:.2f} p_perm={q['p_perm']:.4f} Holm={q['p_holm']:.3f} Fisher={q['p_fisher']:.4f}")

with open(os.path.join(ROOT, "results", "summary.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print("\nWrote results/summary.json")
