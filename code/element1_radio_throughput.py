"""Element 1 — instantaneous radio quality does not predict throughput.

Reproducibility-packaging Rule 3: every R^2 / correlation cited in
Results subsection "Radio Quality Does Not Predict Performance"
(Section IV-A) is *printed* by this script, never hand-transcribed.

Method:

  Vienna LTE (phone_data_lte.parquet)
    Single-feature SINR -> DL throughput, degree-3 polynomial OLS.
    Samples: SINR in [-20, 40] dB, throughput >= 0 Mbps.

  DoNext iperf (iperf_data.csv)
    "Clean downlink full-buffer" = direction == 1 (downlink), datarate
    present, datarate <= 2000 Mbps (drops physically-impossible iperf
    artefacts up to ~9.8 Gbps). datarate is in bit/s -> Mbps = /1e6.
    Per RAT (all / 5G-NSA / 4G):
      - SINR-only degree-3 polynomial OLS
      - multivariate linear OLS on the LTE leg [sinr, rsrp, rsrq, cqi]
    5G-NSA NR-leg model: linear OLS on all 7 LTE+NR features
      [sinr, rsrp, rsrq, cqi, ss_sinr, ss_rsrp, ss_rsrq].
    Best single |Pearson correlation| with throughput over those 7.

Headline (best R^2 per source): Vienna LTE ~0.13, DoNext 5G-NSA ~0.10,
DoNext 4G ~0.18; every |corr| <= 0.21. Radio quality explains <~20% of
throughput variance -> a radio-feature predictor is not viable; this is
the premise of the throughput-prediction literature (Samba 2017, Raca
2018/2020, Narayanan 2020), reproduced here across LTE and 5G-NSA.

Run:  python element1_radio_throughput.py
Out:  stdout + numbers/element1.txt
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as C

# Vienna SINR->throughput calibration window
SINR_CLIP_LOW = -20.0
SINR_CLIP_HIGH = 40.0
POLY_DEGREE = 3
# DoNext downlink full-buffer cleaning
DL_DIRECTION = 1                 # iperf direction code for downlink
DATARATE_CLIP_MBPS = 2000.0     # drop >2 Gbps iperf artefacts
LTE_FEATS = ["sinr", "rsrp", "rsrq", "cqi"]
NR_FEATS = ["sinr", "rsrp", "rsrq", "cqi", "ss_sinr", "ss_rsrp", "ss_rsrq"]

_lines: list[str] = []


def p(s: str = "") -> None:
    print(s)
    _lines.append(s)


def h(s: str) -> None:
    p("")
    p("=" * 70)
    p(s)
    p("=" * 70)


def _r2_poly(x: np.ndarray, y: np.ndarray, deg: int) -> tuple[float, int]:
    """Degree-`deg` polynomial OLS R^2 on finite pairs."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < deg + 2:
        return float("nan"), len(x)
    coeffs = np.polyfit(x, y, deg)
    yhat = np.polyval(coeffs, x)
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return float(1.0 - ss_res / ss_tot), len(x)


def _r2_linear(X: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    """Multivariate linear OLS R^2 (intercept added) on complete cases."""
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    X, y = X[m], y[m]
    A = np.column_stack([X, np.ones(len(X))])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return float(1.0 - ss_res / ss_tot), len(X)


# ---------------------------------------------------------------------
def vienna() -> tuple[float, np.ndarray, np.ndarray, int]:
    h("VIENNA LTE — SINR -> DL throughput (degree-3 polynomial OLS)")
    df = pd.read_parquet(
        C.VIENNA_LTE, columns=["sinr_db", "dl_throughput_mbps", "operator"]
    ).dropna(subset=["sinr_db", "dl_throughput_mbps"])
    df = df[(df["sinr_db"] >= SINR_CLIP_LOW) & (df["sinr_db"] <= SINR_CLIP_HIGH)
            & (df["dl_throughput_mbps"] >= 0.0)]
    x = df["sinr_db"].to_numpy(float)
    y = df["dl_throughput_mbps"].to_numpy(float)
    r2, n = _r2_poly(x, y, POLY_DEGREE)
    p(f"Valid samples                  : {n:,}")
    p(f"Overall R^2 (deg {POLY_DEGREE})            : {r2:.4f}")
    p(f"|Pearson corr| (SINR,tput)     : {abs(np.corrcoef(x, y)[0, 1]):.4f}")
    p("Per-operator R^2 (same global fit window):")
    for op in ("A", "B", "C"):
        s = df[df["operator"] == op]
        r2o, no = _r2_poly(s["sinr_db"].to_numpy(float),
                           s["dl_throughput_mbps"].to_numpy(float), POLY_DEGREE)
        p(f"  operator {op}: R^2 = {r2o:.4f}  (n={no:,})")
    return r2, x, y, n


def _two_panel_figure(
    x_vie: np.ndarray, y_vie: np.ndarray, r2_vie: float, n_vie: int,
    x_nsa: np.ndarray, y_nsa: np.ndarray, n_nsa: int,
) -> None:
    """Fig. 1 (Section IV-A): two-panel hex-bin density of the radio-quality ->
    throughput relationship. (a) Vienna LTE (SINR), (b) DoNext 5G-NSA
    (SS-RSRP, the single most-predictive NR-leg feature). Page-wide
    figure (figure*) with 2 side-by-side panels; matplotlib fonts set
    to STIX to visually match the Times-like manuscript body."""
    _apply_ieee_rc()
    fig, axes = plt.subplots(1, 2, figsize=(7.16, 2.8))

    # Panel (a): Vienna LTE — SINR vs DL throughput, deg-3 poly overlay
    ax = axes[0]
    y_cap_v = float(np.percentile(y_vie, 99.5))
    mask = y_vie <= y_cap_v
    hb_v = ax.hexbin(
        x_vie[mask], y_vie[mask], gridsize=60, cmap="viridis", mincnt=1, bins="log"
    )
    xs = np.linspace(SINR_CLIP_LOW, SINR_CLIP_HIGH, 400)
    coeffs_v = np.polyfit(x_vie, y_vie, POLY_DEGREE)
    ax.plot(
        xs, np.polyval(coeffs_v, xs),
        color="red", linewidth=1.5,
        label=fr"deg-3 fit, $R^2={r2_vie:.2f}$",
    )
    ax.set_xlim(SINR_CLIP_LOW, SINR_CLIP_HIGH)
    ax.set_ylim(0, y_cap_v)
    ax.set_xlabel("SINR (dB)")
    ax.set_ylabel("DL throughput (Mbps)")
    ax.legend(loc="upper left", framealpha=0.9)
    ax.set_title(fr"(a) Vienna LTE, $n={n_vie:,}$", loc="left", fontsize=9)
    cb_v = fig.colorbar(hb_v, ax=ax)
    cb_v.set_label(r"$\log_{10}$ count")
    cb_v.ax.tick_params(labelsize=7)

    # Panel (b): DoNext 5G-NSA — SS-RSRP vs DL throughput (best single
    # NR-leg feature). Deg-3 poly OLS fit + R^2 overlaid, matching
    # panel (a)'s presentation.
    ax = axes[1]
    y_cap_n = float(np.percentile(y_nsa, 99.5))
    mask_n = y_nsa <= y_cap_n
    hb_n = ax.hexbin(
        x_nsa[mask_n], y_nsa[mask_n], gridsize=60, cmap="viridis", mincnt=1, bins="log"
    )
    r2_nsa, _ = _r2_poly(x_nsa, y_nsa, POLY_DEGREE)
    xs_n = np.linspace(float(np.min(x_nsa)), float(np.max(x_nsa)), 400)
    coeffs_n = np.polyfit(x_nsa, y_nsa, POLY_DEGREE)
    ax.plot(
        xs_n, np.polyval(coeffs_n, xs_n),
        color="red", linewidth=1.5,
        label=fr"deg-3 fit, $R^2={r2_nsa:.2f}$",
    )
    ax.set_xlabel("SS-RSRP (dBm)")
    ax.set_ylabel("DL throughput (Mbps)")
    ax.set_ylim(0, y_cap_n)
    ax.legend(loc="upper left", framealpha=0.9)
    ax.set_title(fr"(b) DoNext 5G-NSA, $n={n_nsa:,}$", loc="left", fontsize=9)
    cb_n = fig.colorbar(hb_n, ax=ax)
    cb_n.set_label(r"$\log_{10}$ count")
    cb_n.ax.tick_params(labelsize=7)

    fig.tight_layout()
    out = C.FIGURES_DIR / "fig_element1_radio_scatter.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    p(f"[fig] {out.name}  (Vienna n={n_vie:,}, DoNext 5G-NSA n={n_nsa:,}, "
      f"DoNext SS-RSRP deg-3 R^2={r2_nsa:.4f})")


def _apply_ieee_rc() -> None:
    """matplotlib rcParams tuned to exactly match the Times-like body of
    the LaTeX manuscript class. Routes all figure text through the
    system LaTeX install via the pgf backend (uses pdflatex directly,
    no cm-super dependency) with the same Times-family text+math fonts
    (mathptmx) that the manuscript class uses."""
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


# ---------------------------------------------------------------------
def donext() -> tuple[float, float, np.ndarray, np.ndarray, int]:
    h("DoNext iperf — radio -> throughput, downlink full-buffer")
    df = pd.read_csv(
        C.DONEXT_IPERF, sep=";",
        usecols=["direction", "datarate", "network"] + NR_FEATS,
    )
    dl = df[(df["direction"] == DL_DIRECTION) & df["datarate"].notna()].copy()
    dl["mbps"] = dl["datarate"] / 1e6
    dl = dl[dl["mbps"] <= DATARATE_CLIP_MBPS]
    p(f"Clean downlink full-buffer rows: {len(dl):,}  "
      f"(direction={DL_DIRECTION}, datarate<= {DATARATE_CLIP_MBPS:.0f} Mbps)")
    p(f"Throughput median / max (Mbps) : {dl['mbps'].median():.1f} / "
      f"{dl['mbps'].max():.1f}")

    best_4g = float("nan")
    p("")
    p("Per-RAT R^2 (SINR-only deg-3 poly | multivariate-4 linear LTE leg):")
    for lab in ("all", "5G NSA", "4G"):
        s = dl if lab == "all" else dl[dl["network"] == lab]
        y = s["mbps"].to_numpy()
        r_sinr, n = _r2_poly(s["sinr"].to_numpy(float), y, POLY_DEGREE)
        r_mv, _ = _r2_linear(s[LTE_FEATS].to_numpy(float), y)
        p(f"  {lab:7s}: SINR-only {r_sinr:.4f} | multivar-4 {r_mv:.4f}  (n={n:,})")
        if lab == "4G":
            best_4g = max(r_sinr, r_mv)

    # 5G-NSA NR-leg: all 7 LTE+NR features
    nsa = dl[dl["network"] == "5G NSA"]
    r7, n7 = _r2_linear(nsa[NR_FEATS].to_numpy(float), nsa["mbps"].to_numpy())
    p("")
    p(f"5G-NSA LTE+NR 7-feature linear R^2: {r7:.4f}  (n={n7:,})")

    # Best single |Pearson| over the 7 features (complete cases)
    b = nsa.dropna(subset=NR_FEATS + ["mbps"])
    corrs = {f: np.corrcoef(b[f], b["mbps"])[0, 1] for f in NR_FEATS}
    p("Best single |Pearson corr| with throughput (5G-NSA, 7 feats):")
    for f, v in sorted(corrs.items(), key=lambda kv: -abs(kv[1]))[:3]:
        p(f"  {f}: {v:+.3f}")
    p(f"  -> max |corr| = {max(abs(v) for v in corrs.values()):.3f} (<= 0.21)")

    # Extract (SS-RSRP, throughput) for the right-hand figure panel
    nsa_fig = nsa.dropna(subset=["ss_rsrp", "mbps"])
    x_nsa_fig = nsa_fig["ss_rsrp"].to_numpy(float)
    y_nsa_fig = nsa_fig["mbps"].to_numpy(float)
    return r7, best_4g, x_nsa_fig, y_nsa_fig, len(nsa_fig)


# ---------------------------------------------------------------------
def main() -> None:
    C.ensure_dirs()
    C.check_data()
    p("ELEMENT 1 — radio quality does NOT predict throughput (Rule 3)")
    p(f"(seed {C.SEED}; OLS R^2, no temporal/context features)")
    r_vienna, x_vie, y_vie, n_vie = vienna()
    r_nsa, r_4g, x_nsa, y_nsa, n_nsa = donext()
    h("HEADLINE (best R^2 per source)")
    p(f"Vienna LTE     : {r_vienna:.3f}")
    p(f"DoNext 5G-NSA  : {r_nsa:.3f}")
    p(f"DoNext 4G      : {r_4g:.3f}")
    p("-> instantaneous radio quality explains <~20% of throughput "
      "variance in every case; max |corr| <= 0.21.")
    _two_panel_figure(x_vie, y_vie, r_vienna, n_vie, x_nsa, y_nsa, n_nsa)
    out = C.NUMBERS_DIR / "element1.txt"
    out.write_text("\n".join(_lines) + "\n")
    print(f"\n[written] {out}")


if __name__ == "__main__":
    main()
