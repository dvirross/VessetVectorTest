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
          persistence each cycle repeats the previous interval exactly with
                      probability q, otherwise is drawn from the woman's own
                      observed multiset (q = 0 is exchangeable: size check).
Part C  Sensitivity of the H_W analysis to unequal follow-up: each woman is
        truncated to her first 12 recorded cycles.

Writes results/sim_validation.npz and prints a summary.
"""
import sys, os, time
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(__file__))
from patterns import (load_data, weekly_anchor_set, Counter, perm_global,
                      perm_within, multinomial_iid, two_sided_p, mahalanobis_test,
                      PATTERN_LABELS)

ROOT = os.path.join(os.path.dirname(__file__), "..")
B_REF = 20_000       # reference batch size for Part A
N_FRESH = 10_000     # fresh pseudo-observed datasets for Part A
R_ALT = 400          # datasets per alternative scenario (Part B)
B_IN = 2_000         # within-woman permutations per simulated dataset (Part B)
ALPHA = 0.05

H, wid, bounds = load_data(os.path.join(ROOT, "data", "FilteredData.csv"))
n_w = len(bounds)
weekly = weekly_anchor_set(H)
C = Counter(wid, weekly)
vals, cnt = np.unique(H, return_counts=True)
probs = cnt / cnt.sum()
LO, HI = int(H.min()), int(H.max())
w_mean = np.array([H[s:f].mean() for s, f in bounds])
w_sd = np.array([H[s:f].std(ddof=1) for s, f in bounds])
w_n = np.array([f - s for s, f in bounds])


def holm_reject(p, alpha=ALPHA):
    p = np.asarray(p)
    m = len(p)
    order = np.argsort(p)
    thresh = alpha / (m - np.arange(m))
    rej = np.zeros(m, bool)
    for k, idx in enumerate(order):
        if p[idx] <= thresh[k]:
            rej[idx] = True
        else:
            break
    return rej


# ── Part A ─────────────────────────────────────────────────────────────────────
def draw(null, rng):
    if null == "G":
        return perm_global(H, rng)
    if null == "iid":
        return multinomial_iid(H, rng, vals, probs)
    return perm_within(H, bounds, rng)


def part_a(null, seed):
    rng = np.random.default_rng(seed)
    R1 = np.array([C.count(draw(null, rng)) for _ in range(B_REF)], float)
    R2 = np.array([C.count(draw(null, rng)) for _ in range(B_REF)], float)
    F = np.array([C.count(draw(null, rng)) for _ in range(N_FRESH)], float)
    # marginal two-sided p-values against R1
    up = (np.array([(R1[:, j][None, :] >= F[:, j][:, None]).sum(1) for j in range(5)]).T + 1) / (B_REF + 1)
    lo = (np.array([(R1[:, j][None, :] <= F[:, j][:, None]).sum(1) for j in range(5)]).T + 1) / (B_REF + 1)
    pm = np.minimum(1, 2 * np.minimum(up, lo))
    marg_rate = (pm <= ALPHA).mean(0)
    holm_fwer = np.mean([holm_reject(p).any() for p in pm])
    # joint: plug-in (R1 for mu/Sigma and reference), split (R1 mu/Sigma, R2 reference)
    mu, Sig = R1.mean(0), np.cov(R1.T)
    Sinv = np.linalg.inv(Sig)
    dm = lambda X: np.sqrt(np.einsum("ni,ij,nj->n", X - mu, Sinv, X - mu))
    ref1, ref2, dF = np.sort(dm(R1)), np.sort(dm(R2)), dm(F)
    p_plug = (B_REF - np.searchsorted(ref1, dF, side="left") + 1) / (B_REF + 1)
    p_split = (B_REF - np.searchsorted(ref2, dF, side="left") + 1) / (B_REF + 1)
    return dict(null=null, marg_rate=marg_rate, holm_fwer=holm_fwer,
                joint_plug_05=(p_plug <= .05).mean(), joint_plug_01=(p_plug <= .01).mean(),
                joint_split_05=(p_split <= .05).mean(), joint_split_01=(p_split <= .01).mean(),
                cond=np.linalg.cond(Sig))


# ── Part B ─────────────────────────────────────────────────────────────────────
def gen_ar1(rho, rng):
    out = np.empty(len(H), np.int64)
    for i, (s, f) in enumerate(bounds):
        n = f - s
        sd = max(w_sd[i], 0.5)
        e = rng.standard_normal(n) * sd
        x = np.empty(n)
        x[0] = e[0]
        c = np.sqrt(1 - rho ** 2)
        for k in range(1, n):
            x[k] = rho * x[k - 1] + c * e[k]
        out[s:f] = np.clip(np.rint(w_mean[i] + x), LO, HI)
    return out


def gen_persist(q, rng):
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


def analyse_hw(X, rng):
    obs = C.count(X).astype(float)
    null = np.array([C.count(perm_within(X, bounds, rng)) for _ in range(B_IN)], float)
    pm = np.array([two_sided_p(null[:, j], obs[j]) for j in range(5)])
    DM, pj, _, _ = mahalanobis_test(null, obs)
    return obs, null.mean(0), pm, pj


def part_b_worker(args):
    family, level, seed = args
    rng = np.random.default_rng(seed)
    gen = (lambda: gen_ar1(level, rng)) if family == "ar1" else (lambda: gen_persist(level, rng))
    obs_all, mu_all, pm_all, pj_all = [], [], [], []
    for _ in range(R_ALT):
        obs, mu, pm, pj = analyse_hw(gen(), rng)
        obs_all.append(obs); mu_all.append(mu); pm_all.append(pm); pj_all.append(pj)
    pm_all = np.array(pm_all)
    return dict(family=family, level=level,
                obs_mean=np.mean(obs_all, 0), null_mean=np.mean(mu_all, 0),
                marg_power=(pm_all <= ALPHA).mean(0),
                holm_power=np.mean([holm_reject(p).any() for p in pm_all]),
                joint_power=(np.array(pj_all) <= ALPHA).mean(),
                pj=np.array(pj_all), pm=pm_all)


# ── Part C ─────────────────────────────────────────────────────────────────────
def part_c(seed, max_cycles=12, B=20_000):
    keep = np.concatenate([np.arange(s, min(f, s + max_cycles)) for s, f in bounds])
    Ht = H[keep]
    widt = wid[keep]
    starts = np.flatnonzero(np.r_[True, widt[1:] != widt[:-1]])
    bt = list(zip(starts.tolist(), np.r_[starts[1:], len(Ht)].tolist()))
    Ct = Counter(widt, weekly)
    rng = np.random.default_rng(seed)
    obs = Ct.count(Ht).astype(float)
    null = np.array([Ct.count(perm_within(Ht, bt, rng)) for _ in range(B)], float)
    pm = np.array([two_sided_p(null[:, j], obs[j]) for j in range(5)])
    DM, pj, _, _ = mahalanobis_test(null, obs)
    return dict(n_cycles=len(Ht), obs=obs, mu=null.mean(0), sd=null.std(0, ddof=1),
                z=(obs - null.mean(0)) / null.std(0, ddof=1), p=pm, DM=DM, pj=pj)


if __name__ == "__main__":
    t0 = time.time()
    jobs_b = [("ar1", r, 1000 + i) for i, r in enumerate([0.0, 0.25, 0.5, 0.75])] + \
             [("persist", q, 2000 + i) for i, q in enumerate([0.0, 0.1, 0.2, 0.3])]
    with Pool(3) as pool:
        res_b_async = pool.map_async(part_b_worker, jobs_b)
        res_a = [part_a(n, 100 + k) for k, n in enumerate(["G", "iid", "W"])]
        print(f"Part A done {time.time()-t0:.0f}s", flush=True)
        res_c = part_c(300)
        print(f"Part C done {time.time()-t0:.0f}s", flush=True)
        res_b = res_b_async.get()
    print(f"Part B done {time.time()-t0:.0f}s", flush=True)

    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    np.savez_compressed(os.path.join(ROOT, "results", "sim_validation.npz"),
                        part_a=np.array(res_a, dtype=object), part_b=np.array(res_b, dtype=object),
                        part_c=np.array([res_c], dtype=object), labels=np.array(PATTERN_LABELS),
                        B_REF=B_REF, N_FRESH=N_FRESH, R_ALT=R_ALT, B_IN=B_IN)

    print("\n=== Part A: empirical size at alpha=.05 ===")
    for r in res_a:
        print(f"H_{r['null']:<4} marginal {np.round(r['marg_rate'],3)}  Holm-FWER {r['holm_fwer']:.3f}  "
              f"joint plug-in {r['joint_plug_05']:.3f} (.01: {r['joint_plug_01']:.3f})  "
              f"joint split {r['joint_split_05']:.3f} (.01: {r['joint_split_01']:.3f})  cond {r['cond']:.1f}")
    print("\n=== Part B: H_W power (alpha=.05) ===")
    for r in res_b:
        print(f"{r['family']} {r['level']:<5} obs mean {np.round(r['obs_mean'],1)} null mean {np.round(r['null_mean'],1)}"
              f"\n      marginal power {np.round(r['marg_power'],3)}  Holm-any {r['holm_power']:.3f}  joint {r['joint_power']:.3f}")
    print("\n=== Part C: truncation to first 12 cycles ===")
    print(f"n cycles {res_c['n_cycles']}; obs {res_c['obs']}; mu {np.round(res_c['mu'],2)}; z {np.round(res_c['z'],2)}; "
          f"p {np.round(res_c['p'],3)}; D_M {res_c['DM']:.2f} p {res_c['pj']:.3f}")
