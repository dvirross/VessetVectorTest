"""Build external_replication/notebooks/replication.ipynb programmatically.

Like scripts/build_notebook.py for the main analysis: edit this builder, not the notebook, then
    python external_replication/scripts/build_replication_notebook.py
    (cd external_replication/notebooks && jupyter nbconvert --to notebook --execute --inplace replication.ipynb)
The notebook loads the cached replicates in ../results (produced by run_external.py) and re-derives
every table and figure of the replication from them; the heavy Monte Carlo is not re-run.
"""
import nbformat as nbf
import os

nb = nbf.v4.new_notebook()
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))
code = lambda s: cells.append(nbf.v4.new_code_cell(s.strip("\n")))

md(r"""
# External replication — Utah Creighton Model cohort

Companion notebook for the external replication of *Null-Model Choice for Rule-Defined Pattern Counts in
Clustered Longitudinal Sequences: A Case Study of Jewish Halachic Menstrual Anticipation Rules* (D. Ross).

**Data**: Stanford, J. B., & Najmabadi, S. (2023). *Menstrual Cycles Length of Women in the USA and Canada,
1990–2013*. The Hive, University of Utah. https://doi.org/10.7278/S50d-4gxs-s4hj (CC BY-NC). Three Creighton
Model FertilityCare cohorts (CMFS 1990–96, TTP 2003–06, CEIBA 2009–13); women aged 18–40 with regular
bleeding and no known subfertility.

Counting rules are imported from the repository's `scripts/patterns.py` (the same code as the main analysis);
the analysis engine is `../scripts/analysis_engine.py`, verified to reproduce every cached Fehring result bit for
bit. Heavy Monte Carlo output is cached in `../results/` by `../scripts/run_external.py`; this notebook
re-derives all tables and figures from those replicates.

| Step | Script | Output |
|---|---|---|
| raw deposit → analysed sequences | `../scripts/prepare_utah.py` | `../data/sequences.csv`, `women.csv`, `preprocessing.json` |
| nulls, summary, sensitivity, validation, figures | `../scripts/run_external.py` | `../results/` |
""")

code(r"""
import sys, os, json
sys.path.insert(0, os.path.abspath('../../scripts')); sys.path.insert(0, os.path.abspath('../scripts'))
import numpy as np, pandas as pd, scipy.stats as stats
import matplotlib.pyplot as plt
from IPython.display import display, Markdown
from patterns import (count_all_loop, perm_within, perm_global, multinomial_iid, two_sided_p, holm_adjust,
                      mahalanobis_test, PATTERN_LABELS)
import analysis_engine as E
import external_figures as F

RES, DATA = '../results', '../data'
FIG = os.path.join(RES, 'figures')
pd.set_option('display.width', 160)
""")

md(r"""
## 1. Source file and preprocessing

Rules fixed before any pattern was counted (see `../reports/utah_preprocessing_and_deviations.md`):
cycle length $L$ = recorded `cycle_length` (asserted equal to end − start + 1); $H = L + 1$; the 180 conception
cycles (no recorded length, all terminal) excluded; a woman's record split at every gap between consecutive
recorded cycles (never bridged); each woman contributes her longest run of ≥ 5 consecutive cycles; no
cycle-length cutoff.
""")
code(r"""
raw = pd.read_csv(os.path.join(DATA, 'raw', 'CrMcyclelength_share.csv')); raw.columns = [c.lstrip('﻿') for c in raw.columns]
prep = json.load(open(os.path.join(DATA, 'preprocessing.json')))
print(f"raw deposit: {len(raw):,} rows, {raw.new_id.nunique()} women; recorded lengths {raw.cycle_length.min():.0f}-{raw.cycle_length.max():.0f}; "
      f"conception flags {raw.conception_cycle.value_counts().to_dict()}")
keys = ['conception_cycles_excluded', 'cycles_with_missing_conception_flag_kept', 'gap_breaks', 'gap_breaks_with_cycle_number_jump',
        'women_with_eligible_run', 'women_without_eligible_run', 'women_with_several_eligible_runs', 'n_women_analysed', 'n_cycles_analysed',
        'cycles_per_woman', 'L_min', 'L_max', 'n_cycles_L_below_18', 'n_cycles_L_above_54', 'age', 'cohort_era_counts']
display(pd.Series({k: prep[k] for k in keys}).to_frame('value'))
print('run-length histogram (all runs of consecutive recorded cycles):', prep['run_length_histogram'])
""")
code(r"""
ds = E.Dataset.from_csv(os.path.join(DATA, 'sequences.csv'), 'utah')
women = pd.read_csv(os.path.join(DATA, 'women.csv'))
desc = ds.describe()
display(pd.Series(desc).to_frame('utah'))
print('Week anchors W (H with L a multiple of 7, within the observed range):', desc['weekly_anchors'])
from IPython.display import Image
F.fig_descriptives(ds, FIG); display(Image(os.path.join(FIG, 'figE1_descriptives.png')))
""")

md(r"""
## 2. Counting rules and equivalence check

Same five rules as the main analysis (Haflaga ≥ 3 equal; Dilug constant non-zero first difference; Week two
consecutive **equal** values at a same-weekday anchor $H \in \{22, 29, 36, \ldots\}$; Week-Dilug two consecutive
equal values at a one-weekday-progression anchor $H \in \{23, 30, 37, \ldots\}$ — in these data 40 of the 42 events are
at 30; Dilug-in-Dilug constant non-zero second difference); distinct maximal runs, never across
women. The vectorised counter is checked against the reference loop implementation.
""")
code(r"""
observed = ds.observed()
assert np.array_equal(observed, count_all_loop(ds.H, ds.bounds, ds.weekly))
rng = np.random.default_rng(0)
for _ in range(300):
    for x in (perm_global(ds.H, rng), perm_within(ds.H, ds.bounds, rng), multinomial_iid(ds.H, rng)):
        assert np.array_equal(ds.C.count(x), count_all_loop(x, ds.bounds, ds.weekly))
print('Vectorised and loop counters agree on the observed data and 900 null sequences (run_external.py checks 9,000).')
print('Observed counts:', dict(zip(PATTERN_LABELS, observed.tolist())))
pw = ds.C.count_by_woman(ds.H, ds.n_women); n_i = ds.n_i
print('Women establishing each pattern at least once:', dict(zip(PATTERN_LABELS, (pw > 0).sum(0).tolist())))
""")

md(r"""
## 3. Null replicates

$H_G$ (global permutation), $H_{iid}$ (i.i.d. draws from the pooled empirical distribution) and $H_W$
(within-woman permutation without replacement), $B = 50{,}000$, seed 17, one shared RNG stream — identical
code and constants to the main analysis. Uncomment the last line to regenerate (≈ 1 min).
""")
code(r"""
R = np.load(os.path.join(RES, 'null_replicates.npz'))
assert np.array_equal(R['observed'], observed) and int(R['n_cycles']) == ds.n_cycles
B = int(R['B']); obs = observed.astype(float)
arrs = {'H_G': R['perm'].astype(float), 'H_iid': R['multi'].astype(float), 'H_W': R['within'].astype(float)}
wbw = R['within_by_woman'].astype(np.int32)
print(f'Loaded B={B:,} replicates per null (seed {int(R["seed"])}).')
# R = E.run_nulls(ds, 50_000, 17)   # regenerate
""")

md(r"""
## 4. Marginal tests
Two-sided Monte Carlo $p = \min\{1, 2\min((b^+ + 1)/(B+1), (b^- + 1)/(B+1))\}$; Holm across the five patterns
within each null.
""")
code(r"""
rows = []
for name, arr in arrs.items():
    mu, sd = arr.mean(0), arr.std(0, ddof=1)
    p = np.array([two_sided_p(arr[:, j], obs[j]) for j in range(5)]); adj = holm_adjust(p)
    for j, l in enumerate(PATTERN_LABELS):
        rows.append(dict(Null=name, Pattern=l, Obs=int(obs[j]), Mean=round(mu[j], 2), SD=round(sd[j], 2), z=round((obs[j]-mu[j])/sd[j], 2),
                         p_two_sided=round(p[j], 4), p_Holm=round(adj[j], 4), n_ge_obs=int((arr[:, j] >= obs[j]).sum()), n_le_obs=int((arr[:, j] <= obs[j]).sum())))
marg = pd.DataFrame(rows)
display(marg.style.format({'p_two_sided': '{:.4f}', 'p_Holm': '{:.4f}'}).set_caption('Utah: marginal results'))
summ = json.load(open(os.path.join(RES, 'summary.json')))
F.fig_null_pmfs(summ, R, FIG); display(Image(os.path.join(FIG, 'figE2_null_distributions.png')))
F.fig_zscores(summ, FIG); display(Image(os.path.join(FIG, 'figE3_zscores.png')))
""")

md(r"""
## 5. Joint five-pattern Mahalanobis test (primary), split-batch check, exploratory three-pattern subspace
""")
code(r"""
joint = []
for name, arr in arrs.items():
    DM, p, dm_null, Sig = mahalanobis_test(arr, obs)
    DMs, ps, _, _ = mahalanobis_test(arr[B//2:], obs, ref=arr[:B//2])
    DM3, p3, _, _ = mahalanobis_test(arr[:, [0, 3, 4]], obs[[0, 3, 4]])
    corr = np.corrcoef(arr.T); np.fill_diagonal(corr, 0)
    joint.append(dict(Null=name, D_M=round(DM, 3), p=round(p, 4), n_exceed=int((dm_null >= DM).sum()), D_M_split=round(DMs, 3), p_split=round(ps, 4),
                      D_M_3pattern=round(DM3, 3), p_3pattern=round(p3, 4), cond=round(np.linalg.cond(Sig), 2), max_abs_corr=round(np.abs(corr).max(), 3)))
display(pd.DataFrame(joint).set_index('Null'))
print('H_W null covariance:'); display(pd.DataFrame(np.cov(arrs['H_W'].T), index=PATTERN_LABELS, columns=PATTERN_LABELS).round(2))
F.fig_dm(summ, R, FIG); display(Image(os.path.join(FIG, 'figE4_dm_nulldist.png')))
""")

md(r"""
## 6. Leave-one-woman-out ($H_W$) and follow-up influence
Exact: each woman's contribution is subtracted from the observed vector and from every replicate (additive
per-woman counts), so no re-simulation is needed.
""")
code(r"""
loo, loo_p, loo_pm = E.leave_one_out(obs, arrs['H_W'], wbw, pw)
print({k: v for k, v in loo.items()})
q80 = np.quantile(n_i, .8); long = n_i >= q80
print(f'top follow-up quintile (n_i >= {q80:.0f}): {long.sum()} women hold {n_i[long].sum()/n_i.sum():.1%} of cycles and '
      f'{dict(zip(PATTERN_LABELS, np.round(pw[long].sum(0)/pw.sum(0), 2)))} of events')
print('Spearman(n_i, events per woman):', dict(zip(PATTERN_LABELS, [round(stats.spearmanr(n_i, pw[:, j])[0], 2) for j in range(5)])))
""")

md(r"""
## 7. Sensitivity analyses
* **Harmonised support** 18 ≤ L ≤ 54 (the observed range of the Marquette data; a support restriction, not the
  Marquette inclusion criterion): out-of-range cycles removed, sequences split at those points, blocks of ≥ 5 kept.
* **Truncation** to each woman's first 12 cycles (from the validation run, Part C).
""")
code(r"""
sens = json.load(open(os.path.join(RES, 'sensitivity_L18-54.json')))
print('18-54:', sens['info']); print('observed counts', sens['observed'])
display(pd.DataFrame({m: dict(D_M=sens[m]['DM'], p=sens[m]['p_DM'], D_M_split=sens[m]['DM_split'], p_split=sens[m]['p_DM_split'],
                               p_two_sided=sens[m]['p_two'], Holm=sens[m]['holm']) for m in ['G', 'iid', 'W']}).T)
V = np.load(os.path.join(RES, 'sim_validation.npz'), allow_pickle=True)
c = V['part_c'][0]
print(f"first-12 truncation: {c['n_cycles']} cycles, {c['n_women']} women; obs {c['obs'].astype(int).tolist()}; H_W mean {np.round(c['mu'],2).tolist()}; "
      f"z {np.round(c['z'],2).tolist()}; p {np.round(c['p'],3).tolist()}; D_M {c['DM']:.2f} p {c['pj']:.3f}")
""")

md(r"""
## 8. Simulation validation (same sizes as the paper)
Part A: size of every procedure (10,000 fresh datasets per null vs an independent batch of 20,000 replicates).
Part B: $H_W$ power against AR(1) and persistence alternatives (400 datasets × 2,000 permutations per scenario).
""")
code(r"""
A = list(V['part_a']); Bp = list(V['part_b'])
display(pd.DataFrame([dict(null='H_'+r['null'], **{f'marg_{l}': round(r['marg_rate'][j], 3) for j, l in enumerate(PATTERN_LABELS)},
                           Holm_FWER=round(r['holm_fwer'], 3), joint_plugin=round(r['joint_plug_05'], 3), joint_split=round(r['joint_split_05'], 3),
                           joint_plugin_01=round(r['joint_plug_01'], 4), joint_split_01=round(r['joint_split_01'], 4)) for r in A]).set_index('null'))
display(pd.DataFrame([dict(alternative=f"{r['family']} {r['level']}", **{l: round(r['marg_power'][j], 3) for j, l in enumerate(PATTERN_LABELS)},
                           Holm_any=round(r['holm_power'], 3), joint=round(r['joint_power'], 3)) for r in Bp]).set_index('alternative'))
val = dict(part_a=A, part_b=Bp, part_c=c)
F.fig_power(val, FIG); display(Image(os.path.join(FIG, 'figE6_power.png')))
""")

md(r"""
## 9. Exploratory pairwise analyses (Supporting-Information type)
Per-woman rates, Spearman, two-sided permutation $p$ (50,000, seed 17), Holm over the ten pairs; binary
ever/never and follow-up-stratified robustness. Between-woman associations only.
""")
code(r"""
display(pd.DataFrame(summ['pairwise_rate_spearman']).set_index('pair'))
display(pd.DataFrame(summ['pairwise_binary']).set_index('pair').head(4))
print(summ['pairwise_followup_robustness'])
print('SD comparison (descriptive; expected by construction):', summ['sd_comparison'])
""")

md(r"""
## 10. Side-by-side with the Marquette (Fehring) analysis and decomposition of the global $z$
$z_G = (\mu_W - \mu_G)/\sigma_G + (O - \mu_W)/\sigma_G$: a between-woman heterogeneity gap plus the departure of
the observed count from its own within-woman expectation.
""")
code(r"""
ref = json.load(open('../../results/summary.json')); ref.setdefault('dataset', 'fehring')
ref.setdefault('description', E.Dataset.from_csv('../../data/FilteredData.csv', 'fehring').describe())
tab = []
for m in ['G', 'iid', 'W']:
    for lab, S in [('Fehring', ref), ('Utah', summ)]:
        tab.append(dict(Null='H_'+m, dataset=lab, **{f'z_{l}': S[m]['z'][j] for j, l in enumerate(PATTERN_LABELS)}, D_M=S[m]['DM'], p=S[m]['p_DM']))
display(pd.DataFrame(tab).set_index(['Null', 'dataset']))
for lab, S in [('Fehring', ref), ('Utah', summ)]:
    o = np.array(S['observed'], float); mG, sG, mW = (np.array(S['G'][k]) for k in ('mean', 'sd', 'mean')); mW = np.array(S['W']['mean'])
    print(f"{lab:8s} heterogeneity gap {np.round((mW-mG)/sG, 2)}  within-woman residual {np.round((o-mW)/sG, 2)}  z_G {np.round((o-mG)/sG, 2)}")
F.fig_comparison(ref, summ, FIG); display(Image(os.path.join(FIG, 'figE5_comparison_z.png')))
display(Markdown(open(os.path.join(RES, 'comparison.md')).read()))
""")

md(r"""
## 10b. Regularity strata and per-woman concentration under $H_W$ (main article SI, Section S6)
Same checks as for the Marquette data (`scripts/per_woman_hw.py` in the repository root): $H_W$ tests within
terciles/halves of within-woman SD, and whether the women producing the most events are those the null expects to.
""")
code(r"""
from per_woman_hw import run as per_woman_checks
pw_res = per_woman_checks(os.path.join(DATA, 'sequences.csv'), os.path.join(RES, 'null_replicates.npz'))
""")

md(r"""
## 10c. Sensitivity to the run-selection rule (main article SI, Table S8)
The pre-specified rule keeps each woman's single longest run of $\geq 5$ consecutive cycles. Following the 2026-11-01
preprocessing audit (`../reports/2026-11-01_preprocessing_audit_and_open_questions.md`), `../scripts/all_runs_sensitivity.py`
re-segments the raw deposit, verifies that the baseline reproduces `sequences.csv` and the cached replicates exactly,
and analyses three alternative specifications (all runs $\geq 5$ with $H_W$ within each run; all runs $\geq 5$ with
$H_W$ pooling each woman's values over her runs; all runs $\geq 3$). About 5 min; the cached
`../results/all_runs_sensitivity.json` is displayed (set `RUN_SENSITIVITY = True` to recompute).
""")
code(r"""
RUN_SENSITIVITY = False
if RUN_SENSITIVITY:
    import all_runs_sensitivity; all_runs_sensitivity.main()
ar = json.load(open(os.path.join(RES, 'all_runs_sensitivity.json')))
print(f"B = {ar['B']:,}, seed {ar['seed']}; segmentation: {ar['segmentation']}")
print('baseline reproduces cached replicates:', ar['baseline'].get('reproduces_cached_replicates'))
rows = []
for spec in ('baseline', 'all5_run', 'all5_woman', 'all3_run'):
    r = ar[spec]
    rows.append(dict(spec=spec, runs=r['n_runs'], women=r['n_women'], cycles=r['n_cycles'], observed=r['observed'],
                     **{f'DM_{m}': r[m]['DM'] for m in ('G', 'iid', 'W')}, **{f'p_{m}': r[m]['p_DM'] for m in ('G', 'iid', 'W')}))
display(pd.DataFrame(rows).set_index('spec'))
""")

md(r"""
## 10d. Directional and generic ordering statistics under $H_W$ (main article SI, Section S7; exploratory)
Same statistics as for the Marquette data (`scripts/directional_ordering.py` in the repository root): directional
summaries ($T_{tot}$, $T_z$) from the cached replicates under each null, and three generic serial statistics (lag-1
autocorrelation, mean absolute consecutive difference, proportion of equal consecutive intervals) against fresh
within-woman permutations ($B = 50{,}000$, seed 17; Holm across the three), repeated in the women whose every cycle
lies in 18–54 days. Recomputed here and checked against the cached `../../results/directional_ordering_utah.json`.
""")
code(r"""
from directional_ordering import run as directional_ordering
do_res = directional_ordering(os.path.join(DATA, 'sequences.csv'), os.path.join(RES, 'null_replicates.npz'))
cached = json.load(open('../../results/directional_ordering_utah.json'))
for key in ('directional', 'ordering_HW', 'ordering_HW_women_all_L18_54'):
    assert do_res[key] == cached[key], key
print('matches results/directional_ordering_utah.json')
""")

md(r"""
## 10e. What the replication can and cannot show (main article SI, Section S8)
Utah rows of the heterogeneity-only benchmark (`scripts/heterogeneity_benchmark.py`; 500 synthetic cohorts with
per-woman i.i.d. intervals, so $H_W$ is true and $H_G$ false by construction) and of the independent re-implementation
(`scripts/independent_check.py`, no project imports). Both are computed by the main notebook; the cached JSON files
are displayed here.
""")
code(r"""
hb = json.load(open('../../results/heterogeneity_benchmark.json'))['Utah']
print(f"benchmark (Utah): Haflaga z_G median {hb['haflaga_zG_q025_50_975'][1]} (95% {hb['haflaga_zG_q025_50_975'][0]}–{hb['haflaga_zG_q025_50_975'][2]}, "
      f"min {hb['min_haflaga_zG']}); D_M(G) median {hb['DM_G_q025_50_975'][1]}, D_M(W) median {hb['DM_W_q025_50_975'][1]}; "
      f"reject H_G {hb['reject_G_05']:.3f}, reject H_W {hb['reject_W_05']:.3f}")
print(f"observed Utah: Haflaga z_G {summ['G']['z'][0]:.2f}, D_M(G) {summ['G']['DM']:.2f}, D_M(W) {summ['W']['DM']:.2f}")
ic = json.load(open('../../results/independent_check.json'))['Utah']
assert ic['observed'] == list(summ['observed'])
print(f"independent re-implementation (Utah): counts identical {ic['observed']}; D_M(G) {ic['G']['DM']:.2f} vs {summ['G']['DM']:.2f}; "
      f"D_M(W) {ic['W']['DM']:.2f} vs {summ['W']['DM']:.2f} (p {ic['W']['p_DM']:.3f} vs {summ['W']['p_DM']:.3f})")
""")

md(r"""
## 11. Provenance
""")
code(r"""
prov = json.load(open(os.path.join(RES, 'provenance.json')))
display(pd.Series({k: v for k, v in prov.items() if k != 'preprocessing'}).to_frame('value'))
""")

nb["cells"] = cells
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "notebooks", "replication.ipynb")
nbf.write(nb, out)
print("wrote", os.path.relpath(out))
