"""Derive data/FilteredData.csv from the archived Marquette source file.

Re-implements, step for step, the preprocessing originally carried out in the
author's notebook ``Filtering.ipynb`` (https://github.com/dvirross/PhD), so that the
raw -> filtered step of the analysis is reproducible.

Steps (row and woman counts for the archived file used by the author):
  0. archived file                                   1,665 cycles, 159 women
     - every LengthofCycle is present and lies in 18..54; no cycle-length filter is applied
  1. drop women whose (ClientID, CycleNumber) pairs are duplicated, except nfp8107
     (nfp8106, nfp8109, nfp8114: 16 rows)            1,649 cycles, 156 women
  2. keep women with more than 4 recorded cycles     1,562 cycles, 118 women
  3. nfp8107 appears twice (two near-identical copies of cycles 1-8, differing in the
     recorded length of cycle 2: 27 vs 26 days); the first copy is dropped, keeping the
     second, as in the original notebook                 1,554 cycles, 118 women

Usage:
    python scripts/filter_raw.py data/RawData.csv [OUT.csv]
If OUT.csv is omitted the result is compared with data/FilteredData.csv and the script
exits non-zero on any difference.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MIN_CYCLES = 5          # women with fewer than 5 recorded cycles are dropped
KEEP_DUP_WOMAN = "nfp8107"


def filter_raw(raw: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
    def log(msg):
        if verbose:
            print(msg)

    log(f"archived file: {len(raw):,} cycles, {raw['ClientID'].nunique()} women; "
        f"LengthofCycle {raw['LengthofCycle'].min()}-{raw['LengthofCycle'].max()}, "
        f"missing = {raw['LengthofCycle'].isna().sum()}")

    dup = raw[raw.duplicated(["ClientID", "CycleNumber"], keep=False)]
    dup_women = sorted(set(dup["ClientID"]) - {KEEP_DUP_WOMAN})
    d = raw[~raw["ClientID"].isin(dup_women)].reset_index(drop=True)
    log(f"step 1: removed {len(dup_women)} women with duplicated cycle records "
        f"{dup_women} ({len(raw) - len(d)} rows) -> {len(d):,} cycles, "
        f"{d['ClientID'].nunique()} women")

    f = d.groupby("ClientID").filter(lambda x: len(x) >= MIN_CYCLES).reset_index(drop=True)
    log(f"step 2: kept women with >= {MIN_CYCLES} cycles: removed "
        f"{d['ClientID'].nunique() - f['ClientID'].nunique()} women ({len(d) - len(f)} rows) "
        f"-> {len(f):,} cycles, {f['ClientID'].nunique()} women")

    # step 3: nfp8107 has two copies of cycles 1..8 (rows in file order); drop the first
    idx = f.index[f["ClientID"] == KEEP_DUP_WOMAN]
    n_dup = f.loc[idx].duplicated(["ClientID", "CycleNumber"], keep="last").sum()
    first_copy = idx[:n_dup]
    f = f.drop(first_copy).reset_index(drop=True)
    log(f"step 3: dropped the first of two copies of {KEEP_DUP_WOMAN}'s cycles "
        f"({len(first_copy)} rows) -> {len(f):,} cycles, {f['ClientID'].nunique()} women")

    assert not f.duplicated(["ClientID", "CycleNumber"]).any()
    for _, g in f.groupby("ClientID", sort=False):
        assert (g["CycleNumber"].to_numpy() == range(1, len(g) + 1)).all(), "non-contiguous"
    return f


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    out = filter_raw(pd.read_csv(sys.argv[1]))
    if len(sys.argv) > 2:
        out.to_csv(sys.argv[2], index=False)
        print(f"written {sys.argv[2]}")
    else:
        ref = pd.read_csv(ROOT / "data" / "FilteredData.csv")
        same = out.equals(ref)
        print("identical to data/FilteredData.csv:", same)
        sys.exit(0 if same else 1)
