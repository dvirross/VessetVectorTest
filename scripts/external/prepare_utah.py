"""Build the analysed external sequence file from the University of Utah deposit
"Menstrual Cycles Length of Women in the USA and Canada, 1990-2013"
(Stanford & Najmabadi; Hive, DOI 10.7278/S50d-4gxs-s4hj; CC BY-NC).

Pre-specified rules (fixed on 2026-09-16, before the file was seen):
  1. Cycle length L = start of the next recorded cycle - start of this cycle (days), per
     woman, in chronological order. Haflaga value H = L + 1 as in scripts/patterns.py.
  2. Two successive records are *consecutive* when the next start equals this record's
     end + 1 (if end dates exist); otherwise when next start - start equals the recorded
     cycle length (if a length column exists). A woman's record is split at every gap;
     gaps are never bridged.
  3. Each run of at least MIN_CYCLES consecutive cycles is eligible; a woman contributes
     her single longest eligible run (earliest on ties), as in the Fehring analysis where
     each woman is one contiguous series.
  4. No cycle-length cutoff. Records with unparseable dates or non-positive lengths are
     dropped and counted. Exact duplicate rows (same woman, same start date) are removed
     and counted.
  5. Output schema equals data/FilteredData.csv's analysed columns (ClientID, CycleNumber,
     LengthofCycle) so that patterns.load_data() and its contiguity assertion apply.

Usage: python scripts/external/prepare_utah.py RAW_FILE [--id COL --start COL --end COL --length COL]
Writes data/external/utah/sequences.csv and data/external/utah/preprocessing.json.
"""
import argparse, hashlib, json, os, sys
from collections import Counter
import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
OUT_DIR = os.path.join(ROOT, "data", "external", "utah")
MIN_CYCLES = 5


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def guess(cols, keys):
    for c in cols:
        lc = c.lower()
        if any(k in lc for k in keys):
            return c
    return None


def read_any(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path)
    if ext == ".sav":
        return pd.read_spss(path)
    if ext in (".dta",):
        return pd.read_stata(path)
    return pd.read_csv(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw")
    ap.add_argument("--id"); ap.add_argument("--start"); ap.add_argument("--end"); ap.add_argument("--length")
    ap.add_argument("--cohort", help="optional cohort column to report")
    a = ap.parse_args()
    raw = read_any(a.raw)
    cols = list(raw.columns)
    id_c = a.id or guess(cols, ["subject", "woman", "participant", "user", "id"])
    st_c = a.start or guess(cols, ["start", "begin", "onset", "first_day", "firstday"])
    en_c = a.end or guess(cols, ["end", "last_day", "lastday", "stop"])
    ln_c = a.length or guess(cols, ["length", "cycle_len", "cyclelength", "days"])
    print("columns:", cols)
    print(f"using id={id_c!r} start={st_c!r} end={en_c!r} length={ln_c!r}")
    if id_c is None or st_c is None:
        sys.exit("Could not identify id/start columns; pass --id/--start explicitly.")

    log = dict(source_file=os.path.basename(a.raw), source_sha256=sha256(a.raw), n_rows_raw=int(len(raw)),
               n_women_raw=int(raw[id_c].nunique()), columns=cols, id_col=id_c, start_col=st_c, end_col=en_c,
               length_col=ln_c, min_cycles=MIN_CYCLES, rules="see module docstring")
    df = raw.copy()
    df["_start"] = pd.to_datetime(df[st_c], errors="coerce")
    bad_dates = int(df["_start"].isna().sum())
    df = df[df["_start"].notna()]
    dups = int(df.duplicated([id_c, "_start"]).sum())
    df = df.drop_duplicates([id_c, "_start"]).sort_values([id_c, "_start"]).reset_index(drop=True)
    if en_c:
        df["_end"] = pd.to_datetime(df[en_c], errors="coerce")
    log.update(rows_unparseable_start=bad_dates, duplicate_rows_removed=dups)

    seqs, per_woman = [], []
    n_gap_breaks = n_nonpos = 0
    run_len_hist = Counter()
    for wid_, g in df.groupby(id_c, sort=False):
        starts = g["_start"].to_numpy()
        L = (np.diff(starts) / np.timedelta64(1, "D")).astype(float)          # length of cycles 1..n-1
        # consecutive flag between record k and k+1
        if en_c and g["_end"].notna().all():
            ends = g["_end"].to_numpy()
            consec = (starts[1:] - ends[:-1]) / np.timedelta64(1, "D") == 1
            rule = "next_start == end + 1"
        elif ln_c:
            rec = pd.to_numeric(g[ln_c], errors="coerce").to_numpy()[:-1]
            consec = np.isclose(L, rec)
            rule = "next_start - start == recorded length"
        else:
            consec = np.ones(len(L), bool)
            rule = "no gap information: successive records treated as consecutive"
        runs, cur = [], []
        for k in range(len(L)):
            if L[k] <= 0 or not np.isfinite(L[k]):
                n_nonpos += 1
                if cur: runs.append(cur); cur = []
                continue
            if cur and not consec[k - 1]:
                n_gap_breaks += 1
                runs.append(cur); cur = []
            cur.append(int(round(L[k])))
        if cur: runs.append(cur)
        for r in runs: run_len_hist[len(r)] += 1
        elig = [r for r in runs if len(r) >= MIN_CYCLES]
        per_woman.append(dict(id=str(wid_), n_records=int(len(g)), n_runs=len(runs), n_eligible=len(elig)))
        if elig:
            best = max(elig, key=len)   # max keeps the first maximum -> earliest on ties
            seqs.append((str(wid_), best))
    log.update(gap_rule=rule, gap_breaks=n_gap_breaks, nonpositive_or_missing_lengths=n_nonpos,
               women_with_eligible_run=len(seqs), women_without_eligible_run=len(per_woman) - len(seqs),
               women_with_several_eligible_runs=int(sum(w["n_eligible"] > 1 for w in per_woman)),
               run_length_histogram={int(k): v for k, v in sorted(run_len_hist.items())})
    rows = []
    for wid_, s in seqs:
        for k, L_ in enumerate(s, 1):
            rows.append((wid_, k, L_))
    out = pd.DataFrame(rows, columns=["ClientID", "CycleNumber", "LengthofCycle"])
    if a.cohort and a.cohort in raw.columns:
        coh = raw.drop_duplicates(id_c).set_index(id_c)[a.cohort]
        log["cohort_counts_analysed"] = out.drop_duplicates("ClientID")["ClientID"].map(lambda x: coh.get(x, coh.get(int(x)) if str(x).isdigit() else None)).value_counts(dropna=False).to_dict()
    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_csv(os.path.join(OUT_DIR, "sequences.csv"), index=False)
    log.update(n_cycles_analysed=int(len(out)), n_women_analysed=int(out["ClientID"].nunique()),
               L_min=int(out["LengthofCycle"].min()), L_max=int(out["LengthofCycle"].max()),
               output_sha256=sha256(os.path.join(OUT_DIR, "sequences.csv")))
    with open(os.path.join(OUT_DIR, "preprocessing.json"), "w") as fh:
        json.dump(log, fh, indent=1)
    print(json.dumps({k: v for k, v in log.items() if k != "columns"}, indent=1))


if __name__ == "__main__":
    main()
