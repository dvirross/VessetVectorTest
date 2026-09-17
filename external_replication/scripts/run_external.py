"""One-command external replication of the full analysis on a new set of sequences.

    python external_replication/scripts/run_external.py [NAME] [--csv PATH] [--B 50000] [--seed 17]
                                                          [--nproc 3] [--skip-validation] [--quick]

Input : external_replication/data/sequences.csv (columns ClientID, CycleNumber, LengthofCycle;
        the same schema and contiguity requirement as data/FilteredData.csv), or --csv PATH.
Output: external_replication/results/ for NAME = utah (default); external_replication/results/NAME otherwise:
        null_replicates.npz      B replicates under H_G, H_iid, H_W (+ per-woman H_W counts)
        summary.json             every statistic reported for Fehring, for this dataset
        sensitivity_L18-54.json  harmonised-support sensitivity (18 <= L <= 54, sequences split at gaps)
        sim_validation.npz/.txt  size (Part A), H_W power (Part B), first-12 truncation (Part C)
        comparison.json/.md      side-by-side with the cached Fehring results
        figures/                 figE1..figE6
        provenance.json, run.log
Nothing outside external_replication/ is written. --quick runs a tiny smoke test.
"""
import argparse, hashlib, json, os, platform, subprocess, sys, time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
from analysis_engine import (Dataset, run_nulls, summarise, summarise_null, run_validation,  # noqa: E402
                             print_validation, holm_adjust)
from patterns import PATTERN_LABELS  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
EXT = os.path.join(ROOT, "external_replication")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return None


def to_jsonable(o):
    if isinstance(o, dict):
        return {k: to_jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


def comparison_tables(ref, ext, ref_val, ext_val, sens, out_md):
    """Markdown + dict comparing reference (Fehring) and external summaries."""
    comp = {"reference": ref["dataset"] if "dataset" in ref else "fehring", "external": ext["dataset"]}
    lines = [f"# {ext['dataset']} vs Fehring — side-by-side\n"]
    # descriptives
    rd, ed = ref.get("description", {}), ext["description"]
    lines += ["## Dataset characteristics\n", "| | Fehring | " + ext["dataset"] + " |", "|---|---|---|",
              f"| women | {ref['n_women']} | {ext['n_women']} |", f"| cycles | {ref['n_cycles']} | {ext['n_cycles']} |"]
    for k, lab in [("cycles_per_woman_mean", "cycles/woman mean"), ("cycles_per_woman_sd", "cycles/woman SD"),
                   ("cycles_per_woman_min", "cycles/woman min"), ("cycles_per_woman_max", "cycles/woman max"),
                   ("L_mean", "L mean"), ("L_sd", "L SD"), ("L_median", "L median"), ("L_min", "L min"), ("L_max", "L max"),
                   ("share_L_ge35", "share L >= 35"), ("weekly_anchors", "Week anchors (H)")]:
        lines.append(f"| {lab} | {rd.get(k, '')} | {ed.get(k, '')} |")
    lines += ["", "## Observed five-pattern counts\n", "| Pattern | Fehring | " + ext["dataset"] + " |", "|---|---|---|"]
    for j, l in enumerate(PATTERN_LABELS):
        lines.append(f"| {l} | {ref['observed'][j]} | {ext['observed'][j]} |")
    comp["observed"] = {"fehring": ref["observed"], "external": ext["observed"]}
    for m, title in [("G", "H_G global permutation"), ("iid", "H_iid multinomial"), ("W", "H_W within-woman permutation")]:
        lines += ["", f"## Marginal tests under {title}\n",
                  "| Pattern | Fehring μ (σ) | z | p | Holm | " + ext["dataset"] + " μ (σ) | z | p | Holm |", "|---|---|---|---|---|---|---|---|---|"]
        for j, l in enumerate(PATTERN_LABELS):
            r, e = ref[m], ext[m]
            lines.append(f"| {l} | {r['mean'][j]} ({r['sd'][j]}) | {r['z'][j]:+.2f} | {r['p_two'][j]:.4f} | {r['holm'][j]:.3f} | "
                         f"{e['mean'][j]} ({e['sd'][j]}) | {e['z'][j]:+.2f} | {e['p_two'][j]:.4f} | {e['holm'][j]:.3f} |")
        comp[m] = {"fehring": {k: ref[m][k] for k in ["mean", "sd", "z", "p_two", "holm", "DM", "p_DM", "DM_split", "p_DM_split", "DM3", "p_DM3", "cond"]},
                   "external": {k: ext[m][k] for k in ["mean", "sd", "z", "p_two", "holm", "DM", "p_DM", "DM_split", "p_DM_split", "DM3", "p_DM3", "cond"]}}
    lines += ["", "## Joint tests\n", "| Null | Statistic | Fehring | " + ext["dataset"] + " |", "|---|---|---|---|"]
    for m in ["G", "iid", "W"]:
        for k, lab in [("DM", "5-pattern D_M"), ("p_DM", "5-pattern p"), ("DM_split", "split-batch D_M"), ("p_DM_split", "split-batch p"),
                       ("DM3", "3-pattern D_M (exploratory)"), ("p_DM3", "3-pattern p"), ("cond", "covariance condition number")]:
            lines.append(f"| H_{m} | {lab} | {ref[m][k]} | {ext[m][k]} |")
    if "loo" in ref and "loo" in ext:
        lines += ["", "## Leave-one-woman-out (H_W)\n", "| | Fehring | " + ext["dataset"] + " |", "|---|---|---|",
                  f"| joint p range | {ref['loo']['joint_p_min']}–{ref['loo']['joint_p_max']} | {ext['loo']['joint_p_min']}–{ext['loo']['joint_p_max']} |",
                  f"| # LOO joint p < .05 | {ref['loo']['n_joint_p_below_05']} | {ext['loo']['n_joint_p_below_05']} |",
                  f"| marginal p min | {ref['loo']['marginal_p_min']} | {ext['loo']['marginal_p_min']} |",
                  f"| marginal p max | {ref['loo']['marginal_p_max']} | {ext['loo']['marginal_p_max']} |"]
        comp["loo"] = {"fehring": ref["loo"], "external": ext["loo"]}
    if sens:
        lines += ["", "## Harmonised-support sensitivity (18 <= L <= 54), external dataset\n",
                  f"blocks {sens['info']['n_blocks']} from {sens['info']['women_kept']} women, {sens['info']['n_cycles']} cycles; "
                  f"{sens['info']['cycles_removed']} cycles removed; {sens['info']['women_split_into_several_blocks']} women split at a gap.\n",
                  "| Null | 5-pattern D_M | p | split p | marginal two-sided p (Haf, Dil, Wk, W-D, DiD) | Holm |", "|---|---|---|---|---|---|"]
        for m in ["G", "iid", "W"]:
            s = sens[m]
            lines.append(f"| H_{m} | {s['DM']} | {s['p_DM']} | {s['p_DM_split']} | {s['p_two']} | {s['holm']} |")
        comp["sensitivity_L18_54"] = sens
    if ref_val is not None and ext_val is not None:
        lines += ["", "## Simulation validation\n", "| | Fehring | " + ext["dataset"] + " |", "|---|---|---|"]
        for ra, ea in zip(ref_val["part_a"], ext_val["part_a"]):
            lines.append(f"| size H_{ra['null']}: marginal | {np.round(ra['marg_rate'], 3).tolist()} | {np.round(ea['marg_rate'], 3).tolist()} |")
            lines.append(f"| size H_{ra['null']}: Holm / joint plug-in / split | {ra['holm_fwer']:.3f} / {ra['joint_plug_05']:.3f} / {ra['joint_split_05']:.3f} | "
                         f"{ea['holm_fwer']:.3f} / {ea['joint_plug_05']:.3f} / {ea['joint_split_05']:.3f} |")
        for rb, eb in zip(ref_val["part_b"], ext_val["part_b"]):
            lines.append(f"| power {rb['family']} {rb['level']}: Holm / joint | {rb['holm_power']:.3f} / {rb['joint_power']:.3f} | {eb['holm_power']:.3f} / {eb['joint_power']:.3f} |")
        rc, ec = ref_val["part_c"], ext_val["part_c"]
        lines.append(f"| first-12 truncation: cycles, D_M, p | {rc['n_cycles']}, {rc['DM']:.2f}, {rc['pj']:.3f} | {ec['n_cycles']}, {ec['DM']:.2f}, {ec['pj']:.3f}"
                     + (" (inapplicable: no woman has > 12 cycles)" if ec.get("inapplicable") else "") + " |")
        comp["validation"] = {"fehring": to_jsonable({k: ref_val[k] for k in ["part_a", "part_c"]} | {"part_b": [{k: v for k, v in r.items() if k not in ("pj", "pm")} for r in ref_val["part_b"]]}),
                              "external": to_jsonable({k: ext_val[k] for k in ["part_a", "part_c"]} | {"part_b": [{k: v for k, v in r.items() if k not in ("pj", "pm")} for r in ext_val["part_b"]]})}
    if "pairwise_rate_spearman" in ext:
        lines += ["", "## Exploratory pairwise (per-woman rates, Spearman, permutation p, Holm over 10 pairs)\n",
                  "| Pair | Fehring r / p / Holm | " + ext["dataset"] + " r / p / Holm |", "|---|---|---|"]
        rmap = {q["pair"]: q for q in ref.get("pairwise_rate_spearman", [])}
        for q in ext["pairwise_rate_spearman"]:
            r = rmap.get(q["pair"], {})
            lines.append(f"| {q['pair']} | {r.get('r')} / {r.get('p_perm')} / {r.get('p_holm')} | {q['r']} / {q['p_perm']} / {q['p_holm']} |")
    with open(out_md, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return comp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name", nargs="?", default="utah")
    ap.add_argument("--csv")
    ap.add_argument("--B", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--nproc", type=int, default=3)
    ap.add_argument("--skip-validation", action="store_true")
    ap.add_argument("--quick", action="store_true", help="smoke test: tiny B and validation sizes")
    a = ap.parse_args()

    csv = a.csv or os.path.join(EXT, "data", "sequences.csv")
    out_dir = os.path.join(EXT, "results") if a.name == "utah" else os.path.join(EXT, "results", a.name)
    os.makedirs(out_dir, exist_ok=True)
    logf = open(os.path.join(out_dir, "run.log"), "w")

    def log(msg=""):
        print(msg, flush=True); logf.write(str(msg) + "\n"); logf.flush()

    B = 500 if a.quick else a.B
    Bp = 500 if a.quick else 50_000
    val_kw = dict(B_REF=400, N_FRESH=200, R_ALT=4, B_IN=100) if a.quick else {}
    t0 = time.time()
    ds = Dataset.from_csv(csv, a.name)
    log(f"== {a.name}: {ds.n_women} women, {ds.n_cycles} cycles from {os.path.relpath(csv, ROOT)}")
    log(json.dumps(ds.describe(), indent=1))
    n_checked = ds.verify_counting(n_perm=30 if a.quick else 3000)
    log(f"Vectorised counter agrees with the reference loop on the observed data and {n_checked:,} null sequences.")

    # 1. null replicates
    R = run_nulls(ds, B, a.seed, log=log)
    np.savez_compressed(os.path.join(out_dir, "null_replicates.npz"), **R)
    log(f"null replicates written ({time.time()-t0:.0f}s)")

    # 2. full summary
    summ = summarise(ds, R, Bp=Bp, pair_seed=a.seed, log=log)
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(to_jsonable(summ), fh, indent=1)

    # 3. harmonised support sensitivity
    ds_s, info = ds.restricted_support(18, 54, min_cycles=5)
    sens = {"info": info}
    if ds_s.n_cycles > 0 and ds_s.n_women >= 2:
        log(f"\n== harmonised support 18<=L<=54: {info}")
        Rs = run_nulls(ds_s, B, a.seed, store_by_woman=False, log=None)
        for m, key in [("G", "perm"), ("iid", "multi"), ("W", "within")]:
            sens[m] = summarise_null(Rs[key], Rs["observed"], B, label=f"H_{m} (18-54)", log=log)
        sens["observed"] = Rs["observed"].tolist()
        sens["description"] = ds_s.describe()
    with open(os.path.join(out_dir, "sensitivity_L18-54.json"), "w") as fh:
        json.dump(to_jsonable(sens), fh, indent=1)

    # 4. validation
    val = None
    if not a.skip_validation:
        val = run_validation(ds, n_proc=a.nproc, log=log, **val_kw)
        if ds.n_i.max() <= 12:
            val["part_c"]["inapplicable"] = True
            log("Part C: first-12 truncation inapplicable (no woman has more than 12 cycles); reported unchanged.")
        np.savez_compressed(os.path.join(out_dir, "sim_validation.npz"),
                            part_a=np.array(val["part_a"], dtype=object), part_b=np.array(val["part_b"], dtype=object),
                            part_c=np.array([val["part_c"]], dtype=object), labels=np.array(PATTERN_LABELS),
                            B_REF=val["B_REF"], N_FRESH=val["N_FRESH"], R_ALT=val["R_ALT"], B_IN=val["B_IN"])
        print_validation(val, log=log)

    # 5. comparison with Fehring
    ref = json.load(open(os.path.join(ROOT, "results", "summary.json")))
    ref.setdefault("dataset", "fehring")
    ref.setdefault("description", Dataset.from_csv(os.path.join(ROOT, "data", "FilteredData.csv"), "fehring").describe())
    ref_val = None
    rv_path = os.path.join(ROOT, "results", "sim_validation.npz")
    if val is not None and os.path.exists(rv_path):
        rv = np.load(rv_path, allow_pickle=True)
        ref_val = dict(part_a=list(rv["part_a"]), part_b=list(rv["part_b"]), part_c=rv["part_c"][0])
    comp = comparison_tables(ref, summ, ref_val, val, sens if "G" in sens else None, os.path.join(out_dir, "comparison.md"))
    with open(os.path.join(out_dir, "comparison.json"), "w") as fh:
        json.dump(to_jsonable(comp), fh, indent=1)

    # 6. figures
    import external_figures as F
    F.make_all(ds, summ, R, os.path.join(out_dir, "figures"), summ_ref=ref, val=val)

    # 7. provenance
    import scipy, pandas, matplotlib
    prov = dict(dataset=a.name, sequences_csv=os.path.relpath(csv, ROOT), sequences_sha256=sha256(csv),
                B=B, seed=a.seed, pairwise_B=Bp, validation=val_kw or dict(B_REF=20_000, N_FRESH=10_000, R_ALT=400, B_IN=2_000),
                validation_seeds=dict(part_a=[100, 101, 102], part_b_ar1=[1000, 1001, 1002, 1003], part_b_persist=[2000, 2001, 2002, 2003], part_c=300),
                python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, pandas=pandas.__version__,
                matplotlib=matplotlib.__version__, code_commit=git_commit(), quick=a.quick,
                runtime_s=round(time.time() - t0), timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    pp = os.path.join(os.path.dirname(csv), "preprocessing.json")
    if os.path.exists(pp):
        prov["preprocessing"] = json.load(open(pp))
    with open(os.path.join(out_dir, "provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=1)
    log(f"\nDone in {time.time()-t0:.0f}s -> {os.path.relpath(out_dir, ROOT)}")


if __name__ == "__main__":
    main()
