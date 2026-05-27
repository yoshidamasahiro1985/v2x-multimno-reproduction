"""Dataset characterization statistics for the manuscript Datasets section.

Reproducibility-packaging Rule 3: every descriptive number stated in the
"Datasets and Methodology" prose (row counts, durations, session/date
counts, co-presence, sampling cadence, service-class sizes, MNO splits)
must be *printed* by this script, never hand-transcribed.

Canonical logic:
  - Vienna sessions  = breaks in the union timeline with gap > 60 s.
  - Vienna dates     = distinct calendar dates of `time`.
  - Berlin V2X class = measured_qos == 'delay' (CAM/400 kbps), dual-operator
    simultaneous runs {8,9,10}.

Run:  python dataset_stats.py
Out:  stdout + numbers/dataset_stats.txt
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C


def _haversine_m(lat1, lon1, lat2, lon2):
    """Great-circle distance in metres between matched coordinate arrays."""
    r = 6_371_000.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlmb = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(a))

SESSION_GAP_SEC = 60.0          # union-timeline gap separating drives (src 05)
DUAL_OP_RUNS = (8, 9, 10)       # Berlin V2X delay-class dual-op runs (src 16)

_lines: list[str] = []


def p(s: str = "") -> None:
    print(s)
    _lines.append(s)


def h(s: str) -> None:
    p("")
    p("=" * 70)
    p(s)
    p("=" * 70)


# ---------------------------------------------------------------------
def vienna() -> None:
    h("VIENNA 4G/5G drive-test (phone LTE) — wiedner2026vienna")
    df = pd.read_parquet(
        C.VIENNA_LTE,
        columns=["time", "operator", "rsrp_dbm", "sinr_db"],
    ).sort_values("time").reset_index(drop=True)

    p(f"Total LTE phone rows           : {len(df):,}")
    vc = df["operator"].value_counts()
    for op in ("A", "B", "C"):
        p(f"  operator {op} rows              : {vc.get(op, 0):,} "
          f"({100*vc.get(op,0)/len(df):.1f}%)")

    # Sessions = union-timeline gaps > 60 s
    union_dt = df["time"].diff().dt.total_seconds().fillna(0.0)
    df["session"] = (union_dt > SESSION_GAP_SEC).cumsum().astype(int)
    n_sessions = df["session"].nunique()
    n_dates = df["time"].dt.normalize().nunique()
    p(f"Sessions (union gap > {SESSION_GAP_SEC:.0f} s)   : {n_sessions}")
    p(f"Distinct calendar dates        : {n_dates}")

    # Active measurement time = sum of per-session spans
    spans = df.groupby("session")["time"].agg(lambda s: s.max() - s.min())
    active_h = spans.sum().total_seconds() / 3600.0
    p(f"Active measurement time (h)    : {active_h:.1f}")

    # Co-presence (per calendar date)
    df["date"] = df["time"].dt.normalize()
    by_date = df.groupby("date")["operator"].agg(lambda s: set(s.unique()))
    all3 = sum(1 for s in by_date if {"A", "B", "C"} <= s)
    ac = sum(1 for s in by_date if {"A", "C"} <= s)
    p(f"Dates with all of A,B,C present: {all3}")
    p(f"Dates with both A&C present    : {ac}  (dependence-pair co-presence)")

    # Native cadence: per-operator median inter-sample dt within sessions
    p("Native sampling cadence (median dt within sessions):")
    for op in ("A", "B", "C"):
        o = df[df["operator"] == op].sort_values("time")
        dt = o.groupby("session")["time"].diff().dt.total_seconds()
        dt = dt[(dt > 0) & (dt <= SESSION_GAP_SEC)]
        med = dt.median()
        p(f"  operator {op}: median dt = {med:.3f} s  (~{1/med:.2f} Hz)")

    # DEF_A inputs present (RSRP & SINR)
    p(f"RSRP present (non-NaN)         : {100*df['rsrp_dbm'].notna().mean():.1f}%")
    p(f"SINR present (non-NaN)         : {100*df['sinr_db'].notna().mean():.1f}%")


# ---------------------------------------------------------------------
def berlin() -> None:
    h("BERLIN V2X (cellular) — hernangomez2023berlin")
    df = pd.read_parquet(
        C.BERLIN_CELLULAR,
        columns=["operator", "measured_qos", "measurement", "ping_ms",
                 "Latitude", "Longitude"],
    )
    p(f"Total cellular rows            : {len(df):,}")
    vc = df["operator"].value_counts()
    for op in sorted(vc.index):
        p(f"  operator {op} rows              : {vc[op]:,} "
          f"({100*vc[op]/len(df):.1f}%)")

    qc = df["measured_qos"].value_counts()
    p("Service classes (measured_qos):")
    for q in qc.index:
        p(f"  {q:<10}: {qc[q]:,}")
    p(f"ping_ms null rate (all classes): {100*df['ping_ms'].isna().mean():.1f}%")

    # V2X delay class, dual-op simultaneous runs {8,9,10}
    delay = df[df["measured_qos"] == "delay"]
    core = delay[delay["measurement"].isin(DUAL_OP_RUNS)]
    p(f"V2X delay class total rows     : {len(delay):,}")
    p(f"V2X delay, dual-op runs {DUAL_OP_RUNS}: {len(core):,}")
    cvc = core["operator"].value_counts()
    for op in sorted(cvc.index):
        p(f"  run-8/9/10 operator {op} rows    : {cvc[op]:,}")
    runs_present = sorted(int(r) for r in delay["measurement"].unique())
    p(f"All delay-class runs present   : {runs_present}")

    # Convoy co-location: per-second op1<->op2 great-circle separation (runs 8-10)
    c = core.dropna(subset=["Latitude", "Longitude"]).copy()
    c["sec"] = c.index.floor("s")
    pos = (c.groupby(["sec", "operator"])[["Latitude", "Longitude"]]
             .mean().unstack("operator"))
    pos = pos.dropna()
    sep = _haversine_m(pos[("Latitude", 1)], pos[("Longitude", 1)],
                       pos[("Latitude", 2)], pos[("Longitude", 2)])
    p(f"Convoy op1<->op2 separation (m): median {sep.median():.1f}, "
      f"p90 {sep.quantile(0.9):.1f}  (n={len(sep)} matched seconds)")


# ---------------------------------------------------------------------
def donext() -> None:
    h("DoNext (mobile) — schippers2025donext")
    df = pd.read_csv(
        C.DONEXT_CELL, sep=";",
        usecols=["network", "MNO", "rsrp", "sinr", "ss_rsrp"],
    )
    p(f"Total mobile cell rows         : {len(df):,}")
    vc = df["MNO"].value_counts()
    for mno in sorted(vc.index, key=str):
        p(f"  MNO {mno} rows                 : {vc[mno]:,} "
          f"({100*vc[mno]/len(df):.1f}%)")
    nc = df["network"].value_counts()
    p("RAT (network) split:")
    for net in nc.index:
        p(f"  {net:<8}: {nc[net]:,} ({100*nc[net]/len(df):.1f}%)")
    p(f"5G-NSA samples (ss_rsrp non-NaN): {df['ss_rsrp'].notna().sum():,} "
      f"({100*df['ss_rsrp'].notna().mean():.1f}%)")


# ---------------------------------------------------------------------
def main() -> None:
    C.ensure_dirs()
    C.check_data()
    p("DATASET CHARACTERIZATION — manuscript Datasets section (Rule 3)")
    p(f"(seed {C.SEED}; DEF_A floors RSRP<{C.RSRP_FLOOR_DBM} dBm, "
      f"SINR<{C.SINR_FLOOR_DB} dB)")
    vienna()
    berlin()
    donext()
    out = C.NUMBERS_DIR / "dataset_stats.txt"
    out.write_text("\n".join(_lines) + "\n")
    print(f"\n[written] {out}")


if __name__ == "__main__":
    main()
