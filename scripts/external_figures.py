"""Figures for an external-replication run (same palette and conventions as the paper's
notebook figures). All functions take plain arrays / dicts produced by analysis_engine."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from patterns import PATTERN_LABELS

C_PERM, C_MULTI, C_STRAT, C_OBS = "#2E86AB", "#3DAA5C", "#E07B39", "#C73E1D"
INK, MUTED = "#222222", "#666666"
NULL_STYLE = {"G": (C_PERM, "Global permutation ($H_G$)"), "iid": (C_MULTI, "Multinomial ($H_{iid}$)"),
              "W": (C_STRAT, "Within-woman permutation ($H_W$)")}
plt.rcParams.update({"figure.dpi": 150, "font.family": "serif", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titlesize": 11, "axes.labelsize": 10, "font.size": 9,
                     "axes.edgecolor": "#999999", "grid.color": "#e6e6e6"})


def _save(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, f"{name}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, f"{name}.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def _fmt_p(p):
    return "$p<.001$" if p < .001 else f"$p={p:.3f}$"


def fig_descriptives(ds, out_dir, name="figE1_descriptives"):
    L = ds.H - 1
    n_i = ds.n_i
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    ax.hist(L, bins=range(int(L.min()), int(L.max()) + 2), color=C_PERM, alpha=0.75, edgecolor="white")
    ax.axvline(L.mean(), color=INK, lw=1.2, ls="--")
    ax.set_xlabel("Cycle length $L$ (days)"); ax.set_ylabel("Cycles")
    ax.set_title(f"{ds.name}: {ds.n_cycles:,} cycles, {ds.n_women} women (mean {L.mean():.1f}, SD {L.std(ddof=1):.1f})")
    ax = axes[1]
    ax.hist(n_i, bins=range(int(n_i.min()), int(n_i.max()) + 2), color=C_STRAT, alpha=0.75, edgecolor="white")
    ax.set_xlabel("Recorded cycles per woman $n_i$"); ax.set_ylabel("Women")
    ax.set_title(f"Follow-up: mean {n_i.mean():.1f}, median {np.median(n_i):.0f}, range {n_i.min()}–{n_i.max()}")
    fig.tight_layout(); _save(fig, out_dir, name)


def fig_null_pmfs(summary, R, out_dir, name="figE2_null_distributions"):
    obs = np.asarray(R["observed"])
    arrs = {"G": np.asarray(R["perm"]), "iid": np.asarray(R["multi"]), "W": np.asarray(R["within"])}
    fig, axes = plt.subplots(1, 5, figsize=(16, 3.6), sharey=False)
    for j, (ax, lab) in enumerate(zip(axes, PATTERN_LABELS)):
        lo = min(a[:, j].min() for a in arrs.values()); hi = max(max(a[:, j].max() for a in arrs.values()), obs[j])
        xs = np.arange(lo, hi + 1)
        for k, (m, (col, _)) in enumerate(NULL_STYLE.items()):
            pmf = np.bincount(arrs[m][:, j].astype(int) - lo, minlength=len(xs)) / len(arrs[m])
            ax.step(xs, pmf, where="mid", color=col, lw=1.6)
            ax.text(0.98, 0.95 - 0.11 * k, f"$H_{{{m}}}$ {_fmt_p(summary[m]['p_two'][j])}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8, color=col)
        ax.axvline(obs[j], color=C_OBS, lw=1.8, ls="--")
        ax.set_title(f"{lab} (observed {int(obs[j])})"); ax.set_xlabel("Count"); ax.grid(axis="y")
        if j == 0: ax.set_ylabel("Null probability")
    handles = [Line2D([0], [0], color=c, lw=2, label=l) for c, l in NULL_STYLE.values()] + \
              [Line2D([0], [0], color=C_OBS, lw=1.8, ls="--", label="Observed count")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle(f"{summary['dataset']}: null distributions of the five pattern counts (two-sided Monte Carlo $p$)", y=1.02)
    fig.tight_layout(); _save(fig, out_dir, name)


def fig_zscores(summary, out_dir, name="figE3_zscores"):
    x = np.arange(5); w = 0.26
    fig, ax = plt.subplots(figsize=(9, 3.8))
    for k, (m, (col, lab)) in enumerate(NULL_STYLE.items()):
        z = np.asarray(summary[m]["z"])
        bars = ax.bar(x + (k - 1) * w, z, w, color=col, alpha=0.9, label=lab)
        for b, hz, p in zip(bars, z, summary[m]["holm"]):
            if p < .05:
                ax.text(b.get_x() + b.get_width() / 2, hz + (0.12 if hz >= 0 else -0.3), "*", ha="center", color=INK, fontsize=12)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.axhline(1.96, color=MUTED, lw=0.8, ls=":"); ax.axhline(-1.96, color=MUTED, lw=0.8, ls=":")
    ax.set_xticks(x); ax.set_xticklabels(PATTERN_LABELS); ax.set_ylabel("$z = (O-\\mu)/\\sigma$")
    ax.set_title(f"{summary['dataset']}: standardised departures under the three nulls (* Holm-adjusted $p<.05$)")
    ax.legend(frameon=False, ncol=3, loc="upper right"); ax.grid(axis="y")
    fig.tight_layout(); _save(fig, out_dir, name)


def fig_dm(summary, R, out_dir, name="figE4_dm_nulldist"):
    from patterns import mahalanobis_test
    obs = np.asarray(R["observed"], float)
    arrs = {"G": np.asarray(R["perm"], float), "iid": np.asarray(R["multi"], float), "W": np.asarray(R["within"], float)}
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
    for ax, (m, (col, lab)) in zip(axes, NULL_STYLE.items()):
        DM, p, dm_null, _ = mahalanobis_test(arrs[m], obs)
        ax.hist(dm_null, bins=60, color=col, alpha=0.8, edgecolor="none", density=True)
        ax.axvline(DM, color=C_OBS, lw=1.8, ls="--")
        ax.set_title(f"{lab}\n$D_M={DM:.2f}$, {_fmt_p(p)}"); ax.set_xlabel("$D_M$ of null replicates"); ax.grid(axis="y")
    axes[0].set_ylabel("Density")
    fig.suptitle(f"{summary['dataset']}: five-pattern Mahalanobis distance, observed vs null", y=1.04)
    fig.tight_layout(); _save(fig, out_dir, name)


def fig_comparison(summ_ref, summ_ext, out_dir, name="figE5_comparison_z"):
    """z-scores of the reference (Fehring) and external datasets side by side, per null."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8), sharey=True)
    x = np.arange(5); w = 0.38
    for ax, (m, (col, lab)) in zip(axes, NULL_STYLE.items()):
        ax.bar(x - w / 2, summ_ref[m]["z"], w, color=col, alpha=0.5, label=summ_ref["dataset"], edgecolor=col)
        ax.bar(x + w / 2, summ_ext[m]["z"], w, color=col, alpha=1.0, label=summ_ext["dataset"], hatch="//", edgecolor="white")
        ax.axhline(0, color=MUTED, lw=0.8); ax.axhline(1.96, color=MUTED, lw=0.8, ls=":"); ax.axhline(-1.96, color=MUTED, lw=0.8, ls=":")
        ax.set_xticks(x); ax.set_xticklabels([l.replace("-", "-\n") if len(l) > 10 else l for l in PATTERN_LABELS], fontsize=8)
        ax.set_title(f"{lab}\n$D_M$ {summ_ref[m]['DM']:.2f} ({_fmt_p(summ_ref[m]['p_DM'])}) vs {summ_ext[m]['DM']:.2f} ({_fmt_p(summ_ext[m]['p_DM'])})", fontsize=9)
        ax.legend(frameon=False, fontsize=8); ax.grid(axis="y")
    axes[0].set_ylabel("$z$")
    fig.suptitle("Standardised departures: reference vs external dataset", y=1.04)
    fig.tight_layout(); _save(fig, out_dir, name)


def fig_power(val, out_dir, name="figE6_power"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, fam, xl in zip(axes, ["ar1", "persist"], ["AR(1) coefficient $\\rho$", "Repetition probability $q$"]):
        rs = [r for r in val["part_b"] if r["family"] == fam]
        xs = [r["level"] for r in rs]
        ax.plot(xs, [r["joint_power"] for r in rs], "o-", color=C_OBS, lw=2, label="Joint $D_M$ (5 patterns)")
        ax.plot(xs, [r["holm_power"] for r in rs], "s-", color=C_STRAT, lw=2, label="Holm family (any rejection)")
        for j, (lab, ls) in enumerate(zip(PATTERN_LABELS, [":", "--", "-.", ":", "--"])):
            ax.plot(xs, [r["marg_power"][j] for r in rs], ls, color=MUTED, lw=1, label=lab)
        ax.axhline(0.05, color=MUTED, lw=0.8, ls="-."); ax.set_xlabel(xl); ax.set_ylim(0, 1.02); ax.grid(axis="y")
    axes[0].set_ylabel("Rejection rate at $\\alpha=.05$"); axes[1].legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle("$H_W$ power against within-woman dependence alternatives", y=1.03)
    fig.tight_layout(); _save(fig, out_dir, name)


def make_all(ds, summary, R, out_dir, summ_ref=None, val=None):
    fig_descriptives(ds, out_dir)
    fig_null_pmfs(summary, R, out_dir)
    fig_zscores(summary, out_dir)
    fig_dm(summary, R, out_dir)
    if summ_ref is not None:
        fig_comparison(summ_ref, summary, out_dir)
    if val is not None:
        fig_power(val, out_dir)
