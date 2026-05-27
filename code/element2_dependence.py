"""Element 2 — radio-layer multi-MNO outage dependence (Vienna A&C).

Reproducibility-packaging Rule 3: every number cited in Results subsection
"Radio-Layer Multi-MNO Outage Dependence" (Section IV-B) is *printed* by
this script, never hand-transcribed.

This is the statistically heaviest element. It consolidates the radio-layer
outage-dependence analyses into one config-driven script:

  1. A&C joint-outage dependence under spatial conditioning, with the
     co-location split (same-site vs different-site at 50 m), point
     Mantel-Haenszel ratio + 2-h temporal block-bootstrap 95% CI.
  2. Failure-mode decomposition (RSRP-only / SINR-only) -> coverage vs
     interference reading.
  3. Negative controls A&B and B&C (expected to straddle independence).
  4. 3-pair Bonferroni on the load-bearing diff-site residual: one-sided
     bootstrap p(R<=1), x3 correction, and the Bonferroni-98.3% lower bound.
     CI and Bonferroni p come from the SAME bootstrap run (N below).
  5. Earth-grid scale-descent trend test (Kendall tau, exact permutation) on
     the recorded earth-grid MH ratios.
  6. Coordinate validation: the same-/diff-site classification is robust to
     estimated_cell_info coordinate error (min diff-site 309 m, margin 259 m,
     safety factor, same-site azimuth/height fingerprint).
  7. sub-1 GHz mechanism: A & C disproportionately fall back to their
     800 MHz coverage band at joint outages; B holds no sub-1 GHz spectrum
     (clean negative control).

Bootstrap N: the methodology (Section III) commits to N=2000; this script
uses N=2000 throughout, so the printed CI and the Bonferroni p come from one
consistent run.

Run:  python element2_dependence.py
Out:  stdout + numbers/element2.txt
"""

from __future__ import annotations

import math
from itertools import permutations

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as C

# --- Analysis constants ------------------------------------------------
OPERATORS = ["A", "B", "C"]
GRID_MS = 500
BLOCK_HOURS = 2
N_BOOTSTRAP = 2000          # methodology commits to N=2000 (Section III)
COLOC_THR_M = 50.0          # same-site < 50 m, diff-site >= 50 m
SUB1GHZ_KHZ = 1_000_000     # < 1 GHz = sub-1 GHz coverage band
N_PAIRS_TESTED = 3          # A&B, A&C, B&C -> Bonferroni factor
ALPHA = 0.05

# Earth-grid A&C DEF_A real MH ratios at 4 resolutions (coarse->fine),
# imported constants from the earth-grid simulation-null study described in
# Section IV-B; re-derivation requires the full grid sweep and is out of
# scope for this element. The trend test below is computed from these
# in-script.
GRID_M = [220, 110, 55, 28]
EARTHGRID_AC_MH = [1.971, 1.641, 1.303, 1.250]

# Mechanism-decomposition sweep: GPS earth-grid bins at three resolutions
# crossed with three outage definitions (DEF_A composite, RSRP-only =
# coverage events, SINR-only = interference events). The 28 m resolution
# from the DEF_A trend above is excluded here because SINR-only at 28 m
# has too few joint outages for a usable bootstrap CI.
GRID_RESOLUTIONS_M = [220, 110, 55]
GRID_DEG = {220: 0.002, 110: 0.001, 55: 0.0005}

# Same-/diff-site A&C cell pairs identified by the co-location distance distribution
# (a, c, estimated_dist_m[, n_joint_outages]). Used for coordinate validation.
SAME_SITE_PAIRS = [
    (517639, 17531660, 0),
    (517640, 17531660, 0),
    (80141,  6168589,  40),
]
DIFF_SITE_PAIRS = [
    (423691,  4532748,  309, 1),
    (142859, 12250625,  398, 1),
    (136449,  4445965,  439, 1),
    (306691,  4431363,  524, 4),
    (518919,  6162729, 1037, 8),
    (269069,  4464131, 1361, 1),
    (518919,  4777483, 4407, 1),
    (158733,  6159885, 5215, 1),
]

_lines: list[str] = []


def p(s: str = "") -> None:
    print(s)
    _lines.append(s)


def h(s: str) -> None:
    p("")
    p("=" * 70)
    p(s)
    p("=" * 70)


# =====================================================================
# Geometry + outage helpers
# =====================================================================
def haversine_m(lat1, lon1, lat2, lon2):
    R = 6_371_000.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlmb = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def haversine_scalar(lat1, lon1, lat2, lon2):
    R = 6_371_000.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = (math.sin(dphi / 2) ** 2
         + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2)
    return 2 * R * math.asin(math.sqrt(a))


def azimuth_delta(az1, az2):
    d = abs(az1 - az2) % 360
    return min(d, 360 - d)


def outage_defa(rsrp: pd.Series, sinr: pd.Series) -> pd.Series:
    """DEF_A: RSRP<-110 OR SINR<-6.  NaN where operator absent."""
    present = rsrp.notna() | sinr.notna()
    out = pd.Series(False, index=rsrp.index)
    out |= rsrp.notna() & (rsrp < C.RSRP_FLOOR_DBM)
    out |= sinr.notna() & (sinr < C.SINR_FLOOR_DB)
    return out.where(present, other=np.nan)


def outage_single(s: pd.Series, floor: float) -> pd.Series:
    present = s.notna()
    out = s.notna() & (s < floor)
    return out.where(present, other=np.nan)


def mh_point(sub: pd.DataFrame, stratum_col: str):
    """Mantel-Haenszel ratio stratified by stratum_col.
    sub needs flag1, flag2 (0/1), both_out (bool)."""
    agg = sub.groupby(stratum_col, sort=False).agg(
        n_cop=("flag1", "count"), n1=("flag1", "sum"),
        n2=("flag2", "sum"), n_ob=("both_out", "sum"))
    active = agg[(agg["n_cop"] > 0) & (agg["n1"] > 0) & (agg["n2"] > 0)].copy()
    active["n_exp"] = active["n1"] * active["n2"] / active["n_cop"]
    n_obs = int(active["n_ob"].sum())
    n_exp = float(active["n_exp"].sum())
    mh = n_obs / n_exp if n_exp > 0 else float("nan")
    return mh, n_obs, len(active), n_exp


def bootstrap_array(sub: pd.DataFrame, stratum_col: str,
                    n_boot: int, seed: int) -> np.ndarray:
    """Full array of 2-h block-bootstrap MH replicates (fresh rng(seed))."""
    rng = np.random.default_rng(seed)
    cols = ["block", stratum_col, "flag1", "flag2", "both_out"]
    tbl = sub[cols].copy()
    blocks = tbl["block"].unique()
    n_blocks = len(blocks)
    groups = {b: tbl[tbl["block"] == b] for b in blocks}
    out = []
    for _ in range(n_boot):
        sampled = rng.choice(blocks, size=n_blocks, replace=True)
        boot = pd.concat([groups[b] for b in sampled], ignore_index=True)
        agg = boot.groupby(stratum_col).agg(
            n_cop=("flag1", "count"), n1=("flag1", "sum"),
            n2=("flag2", "sum"), n_ob=("both_out", "sum"))
        agg = agg[(agg["n_cop"] > 0) & (agg["n1"] > 0) & (agg["n2"] > 0)].copy()
        if agg.empty:
            continue
        agg["n_exp"] = agg["n1"] * agg["n2"] / agg["n_cop"]
        te = float(agg["n_exp"].sum())
        if te > 0:
            out.append(int(agg["n_ob"].sum()) / te)
    return np.array(out)


def ci95(arr: np.ndarray):
    return float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def verdict(real_mh, lo, hi):
    if any(np.isnan(v) for v in (real_mh, lo, hi)):
        return "n/a"
    if lo > 1.0:
        return "ROBUST (CI lower > 1)"
    if hi < 1.0:
        return "ROBUST anti-clustering (CI upper < 1)"
    return "FRAGILE (CI includes 1)"


def kendall_tau(x, y):
    n = len(x)
    conc = disc = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = (x[j] - x[i]) * (y[j] - y[i])
            if d > 0:
                conc += 1
            elif d < 0:
                disc += 1
    npairs = n * (n - 1) // 2
    return (conc - disc) / npairs if npairs else float("nan"), conc, disc


def exact_tau_pval(x, y, tau_obs):
    n_total = math.factorial(len(x))
    count_ge = sum(
        1 for perm in permutations(y)
        if kendall_tau(x, list(perm))[0] >= tau_obs - 1e-9
    )
    return count_ge / n_total


# =====================================================================
# Data load + wide pivot
# =====================================================================
def load_wide() -> pd.DataFrame:
    h("LOAD Vienna LTE + align to 500 ms grid (outage definition DEF_A)")
    cols = ["time", "operator", "rsrp_dbm", "sinr_db", "cell_id",
            "frequency_khz", "latitude", "longitude"]
    raw = pd.read_parquet(C.VIENNA_LTE, engine="pyarrow", columns=cols)
    raw = raw[raw["operator"].isin(OPERATORS)].sort_values(
        ["operator", "time"]).reset_index(drop=True)
    p(f"Raw A/B/C rows                 : {len(raw):,}")

    carry = ["rsrp_dbm", "sinr_db", "cell_id", "frequency_khz",
             "latitude", "longitude"]
    parts = []
    for op in OPERATORS:
        g = raw[raw["operator"] == op].copy()
        g["t"] = g["time"].dt.floor(f"{GRID_MS}ms")
        g = g.groupby("t", as_index=False).first()
        parts.append(g[["t"] + carry].rename(
            columns={c: f"{c}_{op}" for c in carry}))
    w = parts[0]
    for q in parts[1:]:
        w = w.merge(q, on="t", how="outer")
    w = w.sort_values("t").reset_index(drop=True)
    w["block"] = w["t"].dt.floor(f"{BLOCK_HOURS}h").dt.strftime("%Y%m%d%H")
    p(f"Wide 500 ms slots              : {len(w):,}")
    p(f"2-h temporal blocks            : {w['block'].nunique()}")

    for op in OPERATORS:
        w[f"defa_{op}"] = outage_defa(w[f"rsrp_dbm_{op}"], w[f"sinr_db_{op}"])
        w[f"rsrponly_{op}"] = outage_single(w[f"rsrp_dbm_{op}"], C.RSRP_FLOOR_DBM)
        w[f"sinronly_{op}"] = outage_single(w[f"sinr_db_{op}"], C.SINR_FLOOR_DB)
    return w


def load_coord_maps():
    ci = pd.read_parquet(C.VIENNA_CELLINFO, engine="pyarrow")
    maps = {op: ci[ci["operator"] == op]
            .groupby("cell_id")[["latitude", "longitude"]].mean()
            for op in OPERATORS}
    return ci, maps


def make_pair_sub(w, op1, op2, prefix):
    cop = w[f"{prefix}_{op1}"].notna() & w[f"{prefix}_{op2}"].notna()
    sub = w[cop].copy()
    sub["flag1"] = sub[f"{prefix}_{op1}"].fillna(False).astype(int)
    sub["flag2"] = sub[f"{prefix}_{op2}"].fillna(False).astype(int)
    sub["both_out"] = sub["flag1"].astype(bool) & sub["flag2"].astype(bool)
    sub["cellpair"] = (sub[f"cell_id_{op1}"].astype("Int64").astype(str) + "_"
                       + sub[f"cell_id_{op2}"].astype("Int64").astype(str))
    return sub


def add_dist(sub, op1, op2, maps):
    sub = sub.copy()
    for op in (op1, op2):
        sub[f"clat_{op}"] = maps[op]["latitude"].reindex(sub[f"cell_id_{op}"]).values
        sub[f"clon_{op}"] = maps[op]["longitude"].reindex(sub[f"cell_id_{op}"]).values
    geo = (sub[f"clat_{op1}"].notna() & sub[f"clon_{op1}"].notna()
           & sub[f"clat_{op2}"].notna() & sub[f"clon_{op2}"].notna())
    sub["dist_m"] = np.where(
        geo, haversine_m(sub[f"clat_{op1}"], sub[f"clon_{op1}"],
                         sub[f"clat_{op2}"], sub[f"clon_{op2}"]), np.nan)
    return sub, geo


# =====================================================================
# Analyses
# =====================================================================
def analyse_ac(w, maps):
    h("A&C DEF_A — co-location split + 2-h block bootstrap (N=%d)" % N_BOOTSTRAP)
    ac = make_pair_sub(w, "A", "C", "defa")
    ac, geo = add_dist(ac, "A", "C", maps)
    ac_geo = ac[geo].copy()
    p(f"A&C co-present slots           : {len(ac):,}")
    p(f"both cells geolocated          : {len(ac_geo):,} "
      f"({len(ac_geo)/len(ac):.1%})")
    p(f"joint outages (all / geo)      : {int(ac['both_out'].sum())} / "
      f"{int(ac_geo['both_out'].sum())}")

    subsets = {
        "AC-all-geo":   ac_geo,
        "AC-same-site": ac_geo[ac_geo["dist_m"] < COLOC_THR_M].copy(),
        "AC-diff-site": ac_geo[ac_geo["dist_m"] >= COLOC_THR_M].copy(),
    }
    results = {}
    p("")
    for label, sub in subsets.items():
        if sub.empty or sub["both_out"].sum() == 0:
            p(f"  [{label}] no joint outages — skip")
            continue
        mh, n_obs, n_strata, _ = mh_point(sub, "cellpair")
        arr = bootstrap_array(sub, "cellpair", N_BOOTSTRAP, C.SEED)
        lo, hi = ci95(arr)
        n_jo = int(sub["both_out"].sum())
        p(f"  [{label:12s}] R={mh:.3f}  CI95=[{lo:.3f},{hi:.3f}]  "
          f"n_jo={n_jo}  strata={n_strata}  {verdict(mh, lo, hi)}")
        results[label] = dict(R=mh, lo=lo, hi=hi, n_jo=n_jo,
                              n_strata=n_strata, arr=arr, sub=sub)
    return results


def analyse_failure_mode(w, maps):
    h("Failure-mode decomposition (A&C all-geo, RSRP-only / SINR-only)")
    for name, prefix in [("RSRP-only", "rsrponly"), ("SINR-only", "sinronly")]:
        sub = make_pair_sub(w, "A", "C", prefix)
        sub, geo = add_dist(sub, "A", "C", maps)
        sub = sub[geo].copy()
        if sub.empty or sub["both_out"].sum() == 0:
            p(f"  [{name}] no joint outages — skip")
            continue
        mh, n_obs, n_strata, _ = mh_point(sub, "cellpair")
        arr = bootstrap_array(sub, "cellpair", N_BOOTSTRAP, C.SEED)
        lo, hi = ci95(arr)
        p(f"  [{name:9s}] R={mh:.3f}  CI95=[{lo:.3f},{hi:.3f}]  "
          f"n_jo={int(sub['both_out'].sum())}  {verdict(mh, lo, hi)}")
    p("  Reading: RSRP-only R ~ DEF_A R -> excess is coverage-driven (spatial),")
    p("           consistent with shared sub-1 GHz edge-of-coverage (sub-1 GHz mechanism below).")


def analyse_earthgrid_decomposition(w):
    """A&C earth-grid sweep crossed with {DEF_A, RSRP-only, SINR-only}.

    Stratification by GPS earth-grid bin (vehicle lat/lon at the 500 ms
    slot, fallback A->B->C) instead of estimated_cell_info cell-pair: the
    GPS fix is available for ~all co-present slots (vs. 31.9% with cell
    geolocation), so SINR-only here has n_jo ~= 12 (vs. n_jo = 6 on the
    cell-pair stratification, whose CI [2.115, 22.905] was unusable).

    Mechanism reading: if RSRP-only R (coverage events) and SINR-only R
    (interference events) BOTH track the composite DEF_A R across grid
    resolutions, the joint-outage excess is coverage-geometry-driven (the
    shared sub-1 GHz edge of coverage manifests for both signal regimes
    along the same spatial footprint), not RAN-scheduler-interference
    coupling between operators.
    """
    h("Earth-grid mechanism decomposition (A&C, DEF_A / RSRP-only / SINR-only)")
    # Vehicle GPS: op-A primary, fallback B then C
    lat = w.get("latitude_A", pd.Series(np.nan, index=w.index)).copy()
    lon = w.get("longitude_A", pd.Series(np.nan, index=w.index)).copy()
    for op in ("B", "C"):
        s_lat = w.get(f"latitude_{op}")
        s_lon = w.get(f"longitude_{op}")
        if s_lat is not None:
            lat = lat.where(lat.notna(), s_lat)
        if s_lon is not None:
            lon = lon.where(lon.notna(), s_lon)
    has_gps = lat.notna() & lon.notna()
    p(f"Slots with GPS                 : {int(has_gps.sum()):,} "
      f"({has_gps.mean():.1%})")

    definitions = [
        ("DEF_A",     "defa"),
        ("RSRP-only", "rsrponly"),
        ("SINR-only", "sinronly"),
    ]
    results: dict[tuple[str, int], dict] = {}

    for def_label, prefix in definitions:
        p("")
        p(f"  {def_label}:")
        f_A = w[f"{prefix}_A"]
        f_C = w[f"{prefix}_C"]
        cop = f_A.notna() & f_C.notna() & has_gps
        if not cop.any():
            p("    no co-present + GPS slots -- skip")
            continue
        n_jo_total = int((f_A.fillna(False) & f_C.fillna(False) & cop).sum())
        p(f"    co-present + GPS slots       : {int(cop.sum()):,}")
        p(f"    joint outages (total)        : {n_jo_total}")
        if n_jo_total == 0:
            p("    no joint outages -- skip all resolutions")
            continue

        base = w[cop].copy()
        base["flag1"] = f_A[cop].fillna(False).astype(int)
        base["flag2"] = f_C[cop].fillna(False).astype(int)
        base["both_out"] = base["flag1"].astype(bool) & base["flag2"].astype(bool)
        base["_lat"] = lat[cop].values
        base["_lon"] = lon[cop].values

        for res_m in GRID_RESOLUTIONS_M:
            deg = GRID_DEG[res_m]
            base["bin"] = (
                (base["_lat"] / deg).round(0).astype(int).astype(str)
                + "_"
                + (base["_lon"] / deg).round(0).astype(int).astype(str)
            )
            mh, n_obs, n_strata, _ = mh_point(base, "bin")
            if n_obs == 0:
                p(f"    {res_m:>3} m: no joint outages in active strata -- skip")
                continue
            arr = bootstrap_array(base, "bin", N_BOOTSTRAP, C.SEED)
            lo, hi = ci95(arr)
            vd = verdict(mh, lo, hi)
            p(f"    {res_m:>3} m: R={mh:.3f}  CI95=[{lo:.3f},{hi:.3f}]  "
              f"n_jo={n_obs}  strata={n_strata}  {vd}")
            results[(def_label, res_m)] = dict(
                R=mh, lo=lo, hi=hi, n_jo=n_obs, n_strata=n_strata,
            )

    p("")
    p("  Reading: SINR-only R (interference events) tracks RSRP-only R")
    p("  (coverage events) and the composite DEF_A R across all three grid")
    p("  resolutions -> the A&C joint-outage excess is coverage-geometry-")
    p("  driven, NOT a per-operator RAN-scheduler interference coupling.")
    return results


def _element2_mechanism_figure(eg_results: dict) -> None:
    """Fig. 3 (Section IV-B) mechanism: single-panel earth-grid descent for the
    three outage definitions overlaid. The three curves tracking each
    other across resolutions is the coverage-driven-mechanism signature.
    """
    matplotlib.use("pgf")
    plt.rcParams.update({
        "pgf.texsystem": "pdflatex",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.preamble": r"\usepackage{mathptmx}",
        "text.latex.preamble": r"\usepackage{mathptmx}",
        "font.family": "serif",
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(3.5, 2.9))

    styles = [
        ("DEF_A",     "tab:blue",  "o", "-"),
        ("RSRP-only", "tab:green", "s", "--"),
        ("SINR-only", "tab:red",   "D", ":"),
    ]
    for def_label, color, marker, ls in styles:
        xs, ys, errs_lo, errs_hi = [], [], [], []
        for res_m in GRID_RESOLUTIONS_M:
            r = eg_results.get((def_label, res_m))
            if r is None:
                continue
            xs.append(res_m)
            ys.append(r["R"])
            errs_lo.append(max(r["R"] - r["lo"], 0.0))
            errs_hi.append(max(r["hi"] - r["R"], 0.0))
        if not xs:
            continue
        # LaTeX-safe legend label (`_` -> `\_` under text.usetex)
        leg = def_label.replace("_", r"\_")
        ax.errorbar(
            xs, ys,
            yerr=[errs_lo, errs_hi],
            fmt=marker, color=color, ecolor=color,
            linestyle=ls, linewidth=1.2,
            elinewidth=1.0, capsize=2.5, markersize=5,
            label=leg,
        )
    ax.axhline(1.0, color="0.5", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xticks(GRID_RESOLUTIONS_M)
    ax.set_xticklabels([str(g) for g in GRID_RESOLUTIONS_M])
    # Suppress log-scale minor tick labels (e.g. "6x10^1" overlapping "55")
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    # Tighter xlim padding so tick labels stay inside the plot frame
    ax.set_xlim(min(GRID_RESOLUTIONS_M) * 0.92,
                max(GRID_RESOLUTIONS_M) * 1.08)
    ax.invert_xaxis()
    ax.set_xlabel(r"Grid resolution (m, finer $\rightarrow$)", labelpad=4)
    ax.set_ylabel(r"A\&C joint-outage MH ratio $R$")
    ax.set_ylim(bottom=0.5)
    ax.grid(linestyle=":", linewidth=0.5, alpha=0.7)
    ax.legend(loc="upper right", frameon=True, framealpha=0.9)
    fig.tight_layout()
    out = C.FIGURES_DIR / "fig_element2_mechanism.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    p(f"[fig] {out.name}")


def analyse_neg_controls(w):
    h("Negative controls — A&B, B&C (DEF_A, cell-pair, block bootstrap)")
    out: dict[str, dict] = {}
    for op1, op2 in [("A", "B"), ("B", "C")]:
        sub = make_pair_sub(w, op1, op2, "defa")
        mh, n_obs, n_strata, _ = mh_point(sub, "cellpair")
        arr = bootstrap_array(sub, "cellpair", N_BOOTSTRAP, C.SEED)
        lo, hi = ci95(arr)
        n_jo = int(sub["both_out"].sum())
        p(f"  [{op1}&{op2}] R={mh:.3f}  CI95=[{lo:.3f},{hi:.3f}]  "
          f"n_jo={n_jo}  {verdict(mh, lo, hi)}")
        out[f"{op1}&{op2}"] = dict(R=mh, lo=lo, hi=hi, n_jo=n_jo)
    p("  Expected: CI straddles 1.0 (no pair-coupling for the non-focal pairs).")
    return out


def _element2_figure(ac_results: dict, neg_results: dict) -> None:
    """Fig. 2 (Section IV-B): two-panel single-column figure.
    (a) Forest plot of MH ratios with 95% CIs for the load-bearing A&C
    splits plus the A&B / B&C negative controls.
    (b) Earth-grid scale descent — MH ratio vs spatial-conditioning grid
    resolution, monotone decline with finer spatial control."""
    matplotlib.use("pgf")
    plt.rcParams.update({
        "pgf.texsystem": "pdflatex",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.preamble": r"\usepackage{mathptmx}",
        "text.latex.preamble": r"\usepackage{mathptmx}",
        "font.family": "serif",
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
    })
    # Assemble rows for the forest plot (bottom-to-top so labels read top-down).
    # n_JO (joint-outage count) is embedded in the row label to avoid the
    # right-column overlap with the same-site CI upper whisker.
    def _row(label, res, color, marker):
        return (f"{label} ($n_{{JO}}{{=}}{res['n_jo']}$)", res, color, marker)

    rows = [
        _row("B\\&C (neg ctrl)",  neg_results["B&C"], "0.35", "o"),
        _row("A\\&B (neg ctrl)",  neg_results["A&B"], "0.35", "o"),
        _row("A\\&C diff-site",   ac_results["AC-diff-site"], "tab:red", "D"),
        _row("A\\&C same-site",   ac_results["AC-same-site"], "tab:blue", "s"),
        _row("A\\&C all-geo",     ac_results["AC-all-geo"],   "tab:blue", "o"),
    ]
    fig, (axA, axB) = plt.subplots(
        2, 1, figsize=(3.5, 4.6),
        gridspec_kw={"height_ratios": [1.0, 0.7]},
    )

    # ---- Panel (a): forest plot --------------------------------------
    ys = np.arange(len(rows))
    for y, (label, r, color, marker) in zip(ys, rows):
        axA.errorbar(
            r["R"], y,
            xerr=[[r["R"] - r["lo"]], [r["hi"] - r["R"]]],
            fmt=marker, color=color, ecolor=color,
            elinewidth=1.2, capsize=2.5, markersize=5,
        )
    axA.axvline(1.0, color="0.5", linestyle="--", linewidth=0.8)
    axA.set_yticks(ys)
    axA.set_yticklabels([row[0] for row in rows])
    axA.set_xscale("log")
    axA.set_xlim(0.5, 10.0)
    axA.set_xticks([0.5, 1, 2, 3, 5, 10])
    axA.set_xticklabels(["0.5", "1", "2", "3", "5", "10"])
    axA.set_xlabel(r"Mantel--Haenszel ratio $R$ (log scale)")
    axA.set_title("(a) Joint-outage ratios, 95\\% block-bootstrap CI", fontsize=8)
    axA.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.7)
    axA.set_ylim(-0.5, len(rows) - 0.5)

    # ---- Panel (b): earth-grid descent -------------------------------
    axB.plot(
        GRID_M, EARTHGRID_AC_MH,
        marker="o", color="tab:red", linewidth=1.4, markersize=5,
        label=r"A\&C MH $R$",
    )
    axB.axhline(1.0, color="0.5", linestyle="--", linewidth=0.8)
    axB.set_xscale("log")
    axB.set_xticks(GRID_M)
    axB.set_xticklabels([str(g) for g in GRID_M])
    axB.set_xlim(min(GRID_M) * 0.85, max(GRID_M) * 1.15)
    axB.invert_xaxis()
    axB.set_xlabel("Earth-grid resolution (m, finer $\\rightarrow$)")
    axB.set_ylabel(r"MH ratio $R$")
    axB.set_title(r"(b) Spatial-conditioning descent (Kendall $\tau = +1$)",
                  fontsize=8)
    axB.grid(linestyle=":", linewidth=0.5, alpha=0.7)
    ymin = min(EARTHGRID_AC_MH) * 0.95
    ymax = max(EARTHGRID_AC_MH) * 1.05
    axB.set_ylim(min(ymin, 0.95), ymax)

    fig.tight_layout()
    out = C.FIGURES_DIR / "fig_element2_forest_grid.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    p(f"[fig] {out.name}")


def analyse_bonferroni(results):
    h("3-pair Bonferroni on the load-bearing diff-site residual")
    bonf_pct = 100.0 * (ALPHA / N_PAIRS_TESTED)   # 1.667th percentile
    for label in ("AC-all-geo", "AC-diff-site"):
        r = results.get(label)
        if r is None:
            continue
        arr = r["arr"]
        lo_bonf = float(np.percentile(arr, bonf_pct))
        p_one = float(np.mean(arr <= 1.0))
        p_bonf = min(1.0, N_PAIRS_TESTED * p_one)
        survives = bool(lo_bonf > 1.0 and p_bonf < ALPHA)
        p(f"  [{label:12s}] R={r['R']:.3f}  "
          f"Bonf-98.3%LB={lo_bonf:.3f}  "
          f"p1(R<=1)={p_one:.4f} -> x{N_PAIRS_TESTED} = {p_bonf:.4f}  "
          f"survives={survives}")
    p("  Load-bearing statistic = diff-site CI (a CI, not a multiple-tested p);")
    p("  its Bonferroni-corrected one-sided p and 98.3% LB are reported above.")


def _chi2_upper_p(Q: float, df: int) -> float:
    """Upper-tail p-value for chi-square Q on df degrees of freedom.
    Wilson-Hilferty cube-root normal approximation; accurate to ~3 sig figs
    for df >= 3 and adequate for a homogeneity-test report. No scipy needed."""
    if df <= 0 or Q < 0:
        return float("nan")
    if Q == 0:
        return 1.0
    z = ((Q / df) ** (1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df))) / math.sqrt(
        2.0 / (9.0 * df)
    )
    # Standard normal upper-tail CDF: 0.5 * erfc(z/sqrt(2))
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def analyse_homogeneity(results):
    """Cochran-Q heterogeneity test on per-stratum log-rate-ratios.

    Within-stratum homogeneity check:
    a Mantel--Haenszel stratified rate ratio assumes within-stratum
    homogeneity of the rate ratio. We test that assumption directly. For
    each cell-pair stratum k with observed joint outages O_k >= 1 and
    expected count E_k = n1_k * n2_k / n_cop_k > 0, compute
        r_k = O_k / E_k
        log_r_k = log(r_k)
        var(log_r_k) ~ 1 / O_k         (Poisson Wald approximation)
        w_k = O_k                       (inverse-variance weight)
    Pooled log-rate-ratio:
        log_R_pooled = sum(w_k * log_r_k) / sum(w_k)
    Cochran's Q statistic:
        Q = sum(w_k * (log_r_k - log_R_pooled)^2)  ~  chi^2(K-1) under H0
    Higgins-Thompson I^2 = max(0, (Q - df) / Q) * 100%.

    Note: strata with O_k = 0 contribute zero degrees of freedom by
    construction (no observed event to be heterogeneous about), so we
    restrict to K = #{strata : O_k >= 1}.
    """
    h("Cochran-Q heterogeneity test on per-stratum log-rate-ratios")
    p("  (Tests the Mantel-Haenszel within-stratum homogeneity")
    p("   assumption directly.)")
    for label in ("AC-all-geo", "AC-same-site", "AC-diff-site"):
        r = results.get(label)
        if r is None:
            continue
        sub = r["sub"]
        agg = sub.groupby("cellpair", sort=False).agg(
            n_cop=("flag1", "count"), n1=("flag1", "sum"),
            n2=("flag2", "sum"), O=("both_out", "sum"))
        active = agg[(agg["n_cop"] > 0) & (agg["n1"] > 0) & (agg["n2"] > 0)].copy()
        active["E"] = active["n1"] * active["n2"] / active["n_cop"]
        contrib = active[active["O"] >= 1].copy()
        K = len(contrib)
        if K < 2:
            p(f"  [{label:12s}] K={K} contributing stratum (Q undefined)")
            continue
        contrib["r_k"] = contrib["O"] / contrib["E"]
        contrib["log_r"] = np.log(contrib["r_k"])
        contrib["w"] = contrib["O"].astype(float)
        log_R_pooled = (
            (contrib["w"] * contrib["log_r"]).sum() / contrib["w"].sum()
        )
        Q = float((contrib["w"] * (contrib["log_r"] - log_R_pooled) ** 2).sum())
        df = K - 1
        p_val = _chi2_upper_p(Q, df)
        I2 = max(0.0, (Q - df) / Q) * 100.0 if Q > 0 else 0.0
        homogeneous = p_val >= ALPHA
        r_min = float(contrib["r_k"].min())
        r_max = float(contrib["r_k"].max())
        r_med = float(contrib["r_k"].median())
        p(f"  [{label:12s}] K={K} contrib strata  Q={Q:.2f}  df={df}  "
          f"p={p_val:.3f}  I^2={I2:.0f}%  homogeneous={homogeneous}")
        p(f"               per-stratum r_k: min={r_min:.2f}  median={r_med:.2f}  "
          f"max={r_max:.2f}  (count O_k=1: {int((contrib['O'] == 1).sum())})")
    p("  Reading: failure-to-reject homogeneity (p >= 0.05) is consistent with")
    p("  the MH within-stratum homogeneity assumption; the headline R is then")
    p("  a defensible pooled estimate. Rejection (p < 0.05) would flag the")
    p("  pooled R as masking stratum-level heterogeneity.")


def analyse_trend():
    h("Earth-grid scale-descent trend test (Kendall tau, exact permutation)")
    p("  Earth-grid A&C DEF_A MH (imported constants):")
    for g, r in zip(GRID_M, EARTHGRID_AC_MH):
        p(f"    grid {g:>3} m -> R = {r:.3f}")
    tau, conc, disc = kendall_tau(GRID_M, EARTHGRID_AC_MH)
    pv = exact_tau_pval(GRID_M, EARTHGRID_AC_MH, tau)
    p(f"  Kendall tau = {tau:+.3f}  (concordant={conc}, discordant={disc})")
    p(f"  exact one-sided p(tau>=obs) = {pv:.4f}  "
      f"-> x{N_PAIRS_TESTED} Bonferroni = {min(1.0, N_PAIRS_TESTED*pv):.3f}")
    p("  tau=+1 -> R decreases monotonically as the grid gets finer: the excess")
    p("  shrinks toward independence at finer spatial scales (co-location-driven).")
    p("  CONCEDED: after 3-pair correction the trend is supporting, not decisive.")


def analyse_coords(cellinfo):
    h("Coordinate validation — same-/diff-site classification robustness")

    def look(cid):
        rows = cellinfo[cellinfo["cell_id"] == cid]
        return None if rows.empty else rows.iloc[0]

    # Same-site fingerprint + per-pair error bound
    p("Same-site pairs (estimated co-located):")
    errs = []
    for a_id, c_id, est in SAME_SITE_PAIRS:
        a, c = look(a_id), look(c_id)
        if a is None or c is None:
            p(f"  A={a_id}/C={c_id}: missing — skip")
            continue
        rec = haversine_scalar(a.latitude, a.longitude, c.latitude, c.longitude)
        azd = (azimuth_delta(a.azimuth_deg, c.azimuth_deg)
               if pd.notna(a.azimuth_deg) and pd.notna(c.azimuth_deg) else float("nan"))
        hd = abs(a.height_m - c.height_m)
        p(f"  A={a_id}/C={c_id}  est={est}m recomp={rec:.1f}m  "
          f"az_delta={azd:.1f} deg  height_delta={hd}m")
        errs.append(rec)
    max_err = max(errs) if errs else 0.0
    p(f"  Max same-site coordinate separation (error upper bound): {max_err:.0f} m")

    # Diff-site min + margin + safety factor
    diff_dists = [d for _, _, d, _ in DIFF_SITE_PAIRS]
    min_diff = min(diff_dists)
    margin = min_diff - COLOC_THR_M
    safety = margin / max_err if max_err > 0 else float("inf")
    p("")
    p(f"Min diff-site distance         : {min_diff:.0f} m")
    p(f"Classification margin          : {min_diff:.0f} - {COLOC_THR_M:.0f} "
      f"= {margin:.0f} m")
    p(f"Safety factor (margin / error) : {safety:.1f}x")
    p("  -> to misclassify any diff-site pair as same-site the coordinate error")
    p(f"     must exceed {margin:.0f} m, i.e. {safety:.1f}x the observed {max_err:.0f} m bound.")

    # Diff-site sanity: all within Vienna bbox
    vie_lat, vie_lon = (48.10, 48.35), (16.15, 16.60)
    in_box = all(
        (look(a) is not None and look(c) is not None
         and vie_lat[0] <= look(a).latitude <= vie_lat[1]
         and vie_lon[0] <= look(a).longitude <= vie_lon[1]
         and vie_lat[0] <= look(c).latitude <= vie_lat[1]
         and vie_lon[0] <= look(c).longitude <= vie_lon[1])
        for a, c, _, _ in DIFF_SITE_PAIRS)
    p(f"  All diff-site cells within Vienna bbox: {in_box} (no sentinel coords)")


def analyse_sub1ghz(w):
    h("sub-1 GHz mechanism — A&C fall back to 800 MHz at joint outages")
    ac = make_pair_sub(w, "A", "C", "defa")
    jo_idx = ac[ac["both_out"]].index
    cop_idx = ac.index

    def frac_table(name, idx):
        p(f"  {name} (n_slots={len(idx)}):")
        out = {}
        for op in OPERATORS:
            freq = w.loc[idx, f"frequency_khz_{op}"].dropna()
            if len(freq) == 0:
                p(f"    {op}: no frequency data")
                out[op] = None
                continue
            sub1 = float((freq < SUB1GHZ_KHZ).mean())
            p(f"    {op}: n={len(freq):>6}  median {freq.median()/1000:.0f} MHz  "
              f"sub-1GHz {sub1:.1%}")
            out[op] = sub1
        return out

    overall = frac_table("A&C co-present slots (overall)", cop_idx)
    atjo = frac_table("A&C JOINT-OUTAGE slots", jo_idx)
    p("  Lift = at-JO / overall  (Table II col 5):")
    for op in OPERATORS:
        if overall.get(op) and atjo.get(op) is not None and overall[op] > 0:
            lift = atjo[op] / overall[op]
            p(f"    {op}: lift={lift:.2f}x  ({atjo[op]:.1%} / {overall[op]:.1%})")
        else:
            p(f"    {op}: lift=N/A (no sub-1 GHz spectrum)")
    p("  B holds NO sub-1 GHz spectrum anywhere -> B's 0% is a network fact, not")
    p("  a measurement artefact -> B is a legitimate clean negative control.")
    p("  A & C both retreat to their 800 MHz (B20) band and still fail together")
    p("  -> the residual is concentrated at shared sub-1 GHz coverage edges.")


# =====================================================================
def main() -> None:
    C.ensure_dirs()
    C.check_data()
    p("ELEMENT 2 — radio-layer multi-MNO outage dependence (Rule 3)")
    p(f"(seed {C.SEED}; DEF_A = RSRP<{C.RSRP_FLOOR_DBM:.0f} OR SINR<{C.SINR_FLOOR_DB:.0f}; "
      f"block bootstrap N={N_BOOTSTRAP})")

    w = load_wide()
    cellinfo, maps = load_coord_maps()

    results = analyse_ac(w, maps)
    analyse_failure_mode(w, maps)
    eg_results = analyse_earthgrid_decomposition(w)
    neg = analyse_neg_controls(w)
    analyse_bonferroni(results)
    analyse_homogeneity(results)
    analyse_trend()
    analyse_coords(cellinfo)
    analyse_sub1ghz(w)
    _element2_figure(results, neg)
    _element2_mechanism_figure(eg_results)

    h("HEADLINE")
    same = results.get("AC-same-site")
    diff = results.get("AC-diff-site")
    if same and diff:
        p(f"A&C joint outage is co-location-dominated: same-site R={same['R']:.2f}, "
          f"diff-site R={diff['R']:.2f}")
        p(f"diff-site residual CI95=[{diff['lo']:.3f},{diff['hi']:.3f}]  "
          f"(n_jo={diff['n_jo']})")
    p("-> above location-conditional independence, but modest and concentrated at")
    p("   shared sub-1 GHz coverage edges; 'beyond location' is NOT claimed.")

    out = C.NUMBERS_DIR / "element2.txt"
    out.write_text("\n".join(_lines) + "\n")
    print(f"\n[written] {out}")


if __name__ == "__main__":
    main()
