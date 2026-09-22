"""Two checks on the within-woman null H_W using the cached per-woman replicate counts.

A. Stratified by within-woman regularity: women are split into terciles (and a median
   split) of the SD of their intervals; marginal and joint H_W tests are recomputed per
   stratum by summing per-woman observed and replicate counts. The SD is a function of the
   multiset, which H_W preserves, so strata are fixed under the null and the test is exact.

B. Concentration of events among women: are the women who produce the most events in the
   data also the ones expected to under H_W?  (i) Spearman correlation across women between
   observed event totals and H_W expected totals; (ii) overlap of the top-k sets;
   (iii) per-woman two-sided Monte Carlo p for her total, with the count of women below .05
   compared with the null expectation and Holm/BH corrections; (iv) a concentration statistic,
   the share of all events held by the top-k women *of each replicate*, compared with the
   observed share (two-sided Monte Carlo p).

Usage: python scripts/per_woman_hw.py [DATA_CSV NULL_NPZ OUT_JSON]
  default: data/FilteredData.csv results/null_replicates.npz results/per_woman_hw.json
  Utah:    external_replication/data/sequences.csv external_replication/results/null_replicates.npz results/per_woman_hw_utah.json
"""
import json, os, sys
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from patterns import load_data, weekly_anchor_set, Counter, two_sided_p, holm_adjust, mahalanobis_test, PATTERN_LABELS

ROOT = os.path.join(os.path.dirname(__file__), "..")
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "FilteredData.csv")
NPZ = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "results", "null_replicates.npz")
OUT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "results", "per_woman_hw.json")
H, wid, bounds = load_data(DATA)
n_w = len(bounds)
C = Counter(wid, weekly_anchor_set(H))
obs_w = C.count_by_woman(H, n_w).astype(int)                 # (n_women, 5)
rep = np.load(NPZ)
W = rep["within_by_woman"].astype(np.int16)                   # (B, 118, 5)
B = W.shape[0]
assert np.array_equal(obs_w.sum(0), rep["observed"]), "per-woman observed counts do not sum to stored observed"
assert W.shape[1] == n_w, "replicate array has a different number of women"
# ordering check: per-woman null means must be zero for every woman whose multiset cannot produce a pattern
# (constant sequences aside); more usefully, they must track the women's own cycle counts in this ordering
_mu = W.sum(2).mean(0); print(f"ordering check: Spearman(per-woman H_W mean total, n_i) = {stats.spearmanr(_mu, np.array([f - s for s, f in bounds]))[0]:.3f} (expect strongly positive)")
sd_w = np.array([H[s:f].std(ddof=1) for s, f in bounds]); n_i = np.array([f - s for s, f in bounds])
out = {"B": int(B), "n_women": n_w}

def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); r = np.empty(m); r[o] = np.arange(1, m + 1)
    q = np.minimum.accumulate((p * m / r)[o][::-1])[::-1]; res = np.empty(m); res[o] = np.minimum(q, 1); return res

def stratum_test(mask, name):
    o = obs_w[mask].sum(0); r = W[:, mask, :].sum(1)
    mu, sd = r.mean(0), r.std(0, ddof=1)
    z = (o - mu) / sd; p = np.array([two_sided_p(r[:, j], o[j]) for j in range(5)])
    res = dict(name=name, n_women=int(mask.sum()), n_cycles=int(n_i[mask].sum()), sd_range=[float(sd_w[mask].min()), float(sd_w[mask].max())],
               observed=o.tolist(), null_mean=mu.round(2).tolist(), null_sd=sd.round(2).tolist(),
               z=z.round(2).tolist(), p_two=p.round(4).tolist(), holm=holm_adjust(p).round(3).tolist())
    try:
        DM, pj, _, S = mahalanobis_test(r, o); res.update(DM=round(DM, 3), p_DM=round(pj, 4), cond=round(float(np.linalg.cond(S)), 1))
    except np.linalg.LinAlgError:
        res.update(DM=None, p_DM=None, cond=None)
    return res

print("=== A. H_W by within-woman regularity (SD of intervals)")
order = np.argsort(sd_w); terc = np.array_split(order, 3)
strata = []
for lab, idx in zip(["most regular tercile", "middle tercile", "least regular tercile"], terc):
    m = np.zeros(n_w, bool); m[idx] = True; strata.append(stratum_test(m, lab))
med = np.median(sd_w)
strata.append(stratum_test(sd_w <= med, "below-median SD (more regular)")); strata.append(stratum_test(sd_w > med, "above-median SD (less regular)"))
for s in strata:
    print(f"{s['name']:32s} n={s['n_women']:3d} women, {s['n_cycles']:4d} cycles, SD {s['sd_range'][0]:.2f}-{s['sd_range'][1]:.2f}")
    print(f"   obs {s['observed']}  mu {s['null_mean']}  z {s['z']}  p {s['p_two']}  Holm {s['holm']}  D_M {s['DM']} p {s['p_DM']} cond {s['cond']}")
out["strata"] = strata

print("\n=== B. Per-woman concentration of events")
tot_o = obs_w.sum(1); tot_r = W.sum(2)                       # (118,), (B,118)
mu_w = tot_r.mean(0); sd_wn = tot_r.std(0, ddof=1)
rho, prho = stats.spearmanr(tot_o, mu_w); print(f"Spearman(observed total per woman, H_W expected total): rho = {rho:.3f} (p = {prho:.2e})")
rho_sd, _ = stats.spearmanr(tot_o, -sd_w); rho_n, _ = stats.spearmanr(tot_o, n_i)
print(f"  observed total vs regularity (-SD): rho = {rho_sd:.3f}; vs n_i: rho = {rho_n:.3f}")
out["spearman_obs_vs_expected"] = dict(rho=round(float(rho), 3), p=float(prho)); out["spearman_obs_vs_negSD"] = round(float(rho_sd), 3); out["spearman_obs_vs_n"] = round(float(rho_n), 3)
for k in (5, 10, 20):
    top_o = set(np.argsort(-tot_o)[:k]); top_e = set(np.argsort(-mu_w)[:k])
    print(f"  top-{k} by observed vs top-{k} by expected: overlap {len(top_o & top_e)}/{k}")
    out[f"top{k}_overlap"] = len(top_o & top_e)
# per-woman two-sided p for her own total
p_w = np.array([two_sided_p(tot_r[:, i], tot_o[i]) for i in range(n_w)])
z_w = np.where(sd_wn > 0, (tot_o - mu_w) / np.where(sd_wn > 0, sd_wn, 1), 0.0)
n05 = int((p_w < .05).sum()); print(f"  women with per-woman two-sided p < .05: {n05} of {n_w} (null expectation ~{0.05*n_w:.1f}); min Holm {holm_adjust(p_w).min():.3f}, min BH q {bh(p_w).min():.3f}")
print(f"  per-woman z: mean {z_w.mean():.3f}, SD {z_w.std(ddof=1):.3f}, max {z_w.max():.2f}, min {z_w.min():.2f}; women with z > 2: {(z_w>2).sum()}, z < -2: {(z_w<-2).sum()}")
out["per_woman"] = dict(n_p_below_05=n05, expected_below_05=round(0.05 * n_w, 1), min_holm=round(float(holm_adjust(p_w).min()), 3), min_bh=round(float(bh(p_w).min()), 3),
                        z_mean=round(float(z_w.mean()), 3), z_sd=round(float(z_w.std(ddof=1)), 3), z_max=round(float(z_w.max()), 2), z_min=round(float(z_w.min()), 2),
                        n_z_gt2=int((z_w > 2).sum()), n_z_lt_m2=int((z_w < -2).sum()))
# top women table
top = np.argsort(-tot_o)[:10]
print("  top-10 women by observed events: obs total | H_W mean (sd) | z | p | n_i | SD | rank of expected")
rank_e = np.empty(n_w, int); rank_e[np.argsort(-mu_w)] = np.arange(1, n_w + 1)
rows = []
for i in top:
    rows.append(dict(woman=int(i), obs=int(tot_o[i]), per_pattern=obs_w[i].tolist(), mu=round(float(mu_w[i]), 2),
                     sd=round(float(sd_wn[i]), 2), z=round(float(z_w[i]), 2), p=round(float(p_w[i]), 3), n=int(n_i[i]),
                     sd_H=round(float(sd_w[i]), 2), expected_rank=int(rank_e[i])))
    r = rows[-1]; print(f"   #{r['woman']:3d}: {r['obs']:2d} {r['per_pattern']} | {r['mu']:5.2f} ({r['sd']:.2f}) | z {r['z']:+.2f} | p {r['p']:.3f} | n {r['n']:2d} | SD {r['sd_H']:.2f} | exp-rank {r['expected_rank']}")
out["top10"] = rows
# concentration statistic: share of all events held by the top-k women of each replicate
for k in (6, 12, 24):
    def share(v): s = np.sort(v, axis=-1)[..., ::-1]; return s[..., :k].sum(-1) / np.maximum(v.sum(-1), 1)
    so = share(tot_o); sr = share(tot_r); p = two_sided_p(sr, so)
    print(f"  share of events held by top-{k} women: observed {so:.3f}; H_W null mean {sr.mean():.3f} (SD {sr.std(ddof=1):.3f}); two-sided p = {p:.3f}")
    out[f"share_top{k}"] = dict(observed=round(float(so), 3), null_mean=round(float(sr.mean()), 3), null_sd=round(float(sr.std(ddof=1)), 3), p_two=round(float(p), 3))
# Gini of per-woman totals
def gini(v):
    v = np.sort(v, axis=-1).astype(float); n = v.shape[-1]; idx = np.arange(1, n + 1)
    return (2 * (idx * v).sum(-1) / (n * np.maximum(v.sum(-1), 1e-12))) - (n + 1) / n
go, gr = gini(tot_o), gini(tot_r); pg = two_sided_p(gr, go)
print(f"  Gini of per-woman event totals: observed {go:.3f}; null mean {gr.mean():.3f} (SD {gr.std(ddof=1):.3f}); two-sided p = {pg:.3f}")
out["gini"] = dict(observed=round(float(go), 3), null_mean=round(float(gr.mean()), 3), null_sd=round(float(gr.std(ddof=1)), 3), p_two=round(float(pg), 3))
# women with zero events: observed vs null
zo = int((tot_o == 0).sum()); zr = (tot_r == 0).sum(1); print(f"  women with no event of any type: observed {zo}; null mean {zr.mean():.1f} (SD {zr.std(ddof=1):.1f}); two-sided p = {two_sided_p(zr, zo):.3f}")
out["zero_event_women"] = dict(observed=zo, null_mean=round(float(zr.mean()), 1), null_sd=round(float(zr.std(ddof=1)), 1), p_two=round(float(two_sided_p(zr, zo)), 3))
json.dump(out, open(OUT, "w"), indent=1)
print("saved", OUT)
