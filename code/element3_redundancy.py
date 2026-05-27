"""Element 3 — redundancy utility under measured dependence.

Reproducibility-packaging Rule 3: every number in "Redundancy Utility Under
Measured Dependence" (Section IV-C) is *printed* by this script, never
hand-transcribed.

Propagates the radio-layer rate-ratio family R in {1.00, 1.30, 1.47, 1.97}
from element 2 (Section IV-B) into the dual-operator availability gain via the
Mantel-Haenszel rate-ratio relation P(joint) = R * p_A * p_C, the measurement
quantity element 2 actually estimates. Shows that the measured dependence
erodes the independence-predicted gain by at most ~1.7%, so simple redundancy
retains the vast majority of its theoretical benefit.

R family (from element2_dependence.py; N=2000 block bootstrap):
  1.00 — independence baseline (P(joint) = p_A * p_C)
  1.30 — diff-site characterised point estimate
  1.47 — block-bootstrap diff-site R point estimate (18 JO)
  1.97 — block-bootstrap 97.5th-pct CI upper bound

Berlin numbers (element4_berlin.py is the authoritative source):
  budget=200ms, runs 8/9/10 pooled; JOINT=0 -> observed gain is the full
  single-operator outage rate; we report this alongside the independence-
  predicted gain to avoid the cross-framework labeling confusion that R=0%
  in the Vienna OR-rows would otherwise suggest.

Run:  python element3_redundancy.py
Out:  stdout + numbers/element3.txt
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as C

# --- Analysis constants ---------------------------------------------------
OPERATORS = ["A", "B", "C"]
GRID_MS = 500

# Rate-ratio R family from Section IV-B (Mantel-Haenszel rate ratio, not odds ratio)
R_VALUES = [1.00, 1.30, 1.47, 1.97]
R_LABELS = {
    1.00: "independence baseline",
    1.30: "diff-site char. point est.",
    1.47: "block-bootstrap diff-site R (18 JO, N=2000)",
    1.97: "block-bootstrap 97.5th-pct CI upper",
}

# Berlin: hardcoded constants (budget=200ms, runs 8/9/10 pooled).
# element4_berlin.py is the authoritative source; hardcoded here to keep
# element3 self-contained (no data-load dependency on the Berlin parquet).
_B_N_OP1_OUT   = 25
_B_N_OP1_TOTAL = 10690   # matched seconds with op1 ping present
_B_N_OP2_OUT   = 219
_B_N_TOTAL     = 10574   # matched seconds (both ops have ping)
_B_JOINT_OBS   = 0

BERLIN_P_OP1  = _B_N_OP1_OUT  / _B_N_OP1_TOTAL
BERLIN_P_OP2  = _B_N_OP2_OUT  / _B_N_TOTAL
BERLIN_JOINT  = _B_JOINT_OBS  / _B_N_TOTAL   # = 0.0


# =============================================================================
# Helpers
# =============================================================================

def joint_from_r(p_a: float, p_b: float, r_val: float) -> float:
    """P(A∩B) from marginals and Mantel-Haenszel rate ratio R.

    The rate ratio is R = P(joint) / (p_A * p_B), so P(joint) = R * p_A * p_B.
    This is the relation element 2 actually estimates (a risk/rate ratio of
    observed joint outages to the location-conditional-independence
    prediction), not the odds ratio. At small marginals (p ~ 1-2%) the two
    quantities differ by < 5% relative, but the labels must be honest.

    Capped at min(p_a, p_b) to keep P(A∩B) <= min(P(A), P(B)).
    """
    return min(r_val * p_a * p_b, min(p_a, p_b))


def gain_and_erosion(p_a: float, p_b: float, r_val: float,
                     gain_indep: float) -> tuple[float, float, float]:
    """Return (p_joint, gain_ppm, erosion_pct) for the given rate ratio R."""
    p_j  = joint_from_r(p_a, p_b, r_val)
    gain = min(p_a, p_b) - p_j          # avail_dual - avail_best
    eros = (gain_indep - gain) / gain_indep * 100.0 if gain_indep > 0 else 0.0
    return p_j, gain * 1e6, eros


# =============================================================================
# Step 1: Load Vienna LTE — compute marginal outage rates (DEF_A)
# =============================================================================
C.ensure_dirs()

print("=" * 70)
print("Element 3: Redundancy utility under measured dependence")
print("=" * 70)
print(f"\nLoading Vienna LTE: {C.VIENNA_LTE}")
C.check_data()

df = pd.read_parquet(
    C.VIENNA_LTE, engine="pyarrow",
    columns=["time", "operator", "rsrp_dbm", "sinr_db"],
)
df = df[df["operator"].isin(OPERATORS)].copy()
df["t"] = df["time"].dt.floor(f"{GRID_MS}ms")
df_al = (
    df.sort_values(["operator", "t", "time"])
      .groupby(["operator", "t"], as_index=False)
      .first()
)
df_w = df_al.set_index(["t", "operator"]).unstack("operator")
df_w.columns = [f"{m}_{op}" for m, op in df_w.columns]
df_w = df_w.reset_index()
print(f"  Wide slots: {len(df_w):,}")

marg: dict[str, float] = {}
for op in OPERATORS:
    rsrp    = df_w.get(f"rsrp_dbm_{op}")
    sinr    = df_w.get(f"sinr_db_{op}")
    present = rsrp.notna() | sinr.notna()
    out_flag = (
        (rsrp.fillna(0) < C.RSRP_FLOOR_DBM) | (sinr.fillna(0) < C.SINR_FLOOR_DB)
    ).where(present, other=pd.NA)
    n_pres = int(present.sum())
    n_out  = int(out_flag.dropna().astype(bool).sum())
    marg[op] = n_out / n_pres
    print(f"  Operator {op}: {n_out}/{n_pres} = {marg[op]*100:.4f}%")

p_a, p_c = marg["A"], marg["C"]


# =============================================================================
# Step 2: Vienna A&C availability metrics per OR
# =============================================================================
print("\n--- Vienna A&C: availability under R (rate-ratio) family ---")
print(f"  {'R':<8} {'label':<47} {'P(joint)%':>14} {'Gain(ppm)':>11} {'Erosion%':>9}")
print("  " + "-" * 93)

p_j_indep  = joint_from_r(p_a, p_c, 1.0)
gain_indep = (min(p_a, p_c) - p_j_indep) * 1e6   # ppm

rows: list[dict] = []
for r_val in R_VALUES:
    p_j, gain_ppm, eros = gain_and_erosion(p_a, p_c, r_val, gain_indep / 1e6)
    rows.append(dict(r_val=r_val, p_joint=p_j, gain_ppm=gain_ppm, erosion_pct=eros))
    print(f"  {r_val:<8.2f} {R_LABELS[r_val]:<47} {p_j*100:>13.6f}% "
          f"{gain_ppm:>11.4f} {eros:>9.2f}%")

r_indep = rows[0]
r_130   = rows[1]
r_147   = rows[2]
r_197   = rows[3]


# =============================================================================
# Step 3: Berlin (hardcoded constants)
# =============================================================================
print("\n--- Berlin app-layer (budget=200ms, runs 8/9/10 pooled) ---")
p_j_b_indep  = BERLIN_P_OP1 * BERLIN_P_OP2
avail_best_b = 1.0 - min(BERLIN_P_OP1, BERLIN_P_OP2)
gain_b_indep = ((1.0 - p_j_b_indep) - avail_best_b) * 1e6    # ppm
gain_b_obs   = ((1.0 - BERLIN_JOINT)  - avail_best_b) * 1e6  # ppm (JOINT=0)
op2_ratio    = BERLIN_P_OP2 / BERLIN_P_OP1

print(f"  op1={BERLIN_P_OP1*100:.4f}%  op2={BERLIN_P_OP2*100:.4f}%  op2/op1={op2_ratio:.0f}×")
print(f"  gain (independence-predicted): {gain_b_indep:.1f} ppm = {gain_b_indep/1e4:.4f}%")
print(f"  gain (observed, JOINT=0):       {gain_b_obs:.1f} ppm = {gain_b_obs/1e4:.4f}%")
b_obs_minus_indep = gain_b_obs - gain_b_indep
print(f"  observed minus independence:    {b_obs_minus_indep:+.1f} ppm (observed is at or above independence)")
print(f"  break-even (observed): {gain_b_obs/1e4:.4f}%")


# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 70)
print("KEY NUMBERS FOR PAPER (Section IV-C)")
print("=" * 70)
print(f"Vienna A marginal outage: {marg['A']*100:.4f}%")
print(f"Vienna C marginal outage: {marg['C']*100:.4f}%")
print(f"Vienna single-best availability: {(1-min(p_a,p_c))*100:.4f}%")
print(f"Gain under independence (R=1): {r_indep['gain_ppm']:.2f} ppm")
print(f"Gain under R=1.30: {r_130['gain_ppm']:.2f} ppm  "
      f"erosion={r_130['erosion_pct']:.2f}%")
print(f"Gain under R=1.47: {r_147['gain_ppm']:.2f} ppm  "
      f"erosion={r_147['erosion_pct']:.2f}%")
print(f"Gain under R=1.97: {r_197['gain_ppm']:.2f} ppm  "
      f"erosion={r_197['erosion_pct']:.2f}%")
retains = 100.0 - r_197["erosion_pct"]
print(f"  → retains {retains:.1f}% of independence-predicted benefit at CI upper")
print(f"Break-even Vienna (R=1.30): {r_130['gain_ppm']:.2f} ppm  "
      f"({r_130['gain_ppm']/1e4:.4f}%)")
print(f"Berlin gain (independence-predicted): {gain_b_indep:.1f} ppm  "
      f"({gain_b_indep/1e4:.4f}%)")
print(f"Berlin gain (observed, JOINT=0): {gain_b_obs:.1f} ppm  "
      f"({gain_b_obs/1e4:.4f}%)")
print(f"Berlin break-even (observed): {gain_b_obs/1e4:.4f}%")


# =============================================================================
# Write numbers/element3.txt
# =============================================================================
lines_out = [
    f"Vienna A marginal outage: {marg['A']*100:.4f}%",
    f"Vienna C marginal outage: {marg['C']*100:.4f}%",
    f"Vienna single-best availability: {(1-min(p_a,p_c))*100:.4f}%",
    f"gain independence (R=1): {r_indep['gain_ppm']:.2f} ppm",
    f"gain R=1.30: {r_130['gain_ppm']:.2f} ppm",
    f"erosion R=1.30: {r_130['erosion_pct']:.2f}%",
    f"gain R=1.47: {r_147['gain_ppm']:.2f} ppm",
    f"erosion R=1.47: {r_147['erosion_pct']:.2f}%",
    f"gain R=1.97: {r_197['gain_ppm']:.2f} ppm",
    f"erosion R=1.97: {r_197['erosion_pct']:.2f}%",
    f"retains at R=1.97 CI upper: {retains:.1f}%",
    f"break-even Vienna R=1.30: {r_130['gain_ppm']:.2f} ppm ({r_130['gain_ppm']/1e4:.4f}%)",
    f"Berlin gain independence-predicted: {gain_b_indep:.1f} ppm ({gain_b_indep/1e4:.4f}%)",
    f"Berlin gain observed JOINT=0: {gain_b_obs:.1f} ppm ({gain_b_obs/1e4:.4f}%)",
    f"Berlin observed minus independence: {b_obs_minus_indep:+.1f} ppm",
    f"Berlin break-even (observed): {gain_b_obs/1e4:.4f}%",
    f"Berlin op2/op1 ratio: {op2_ratio:.0f}x",
]

out_path = C.NUMBERS_DIR / "element3.txt"
out_path.write_text("\n".join(lines_out) + "\n")
print(f"\nNumbers written to {out_path}")


# =============================================================================
# Fig. 4 (Section IV-C): availability-gain erosion vs R
# =============================================================================
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
rs_curve = np.linspace(1.0, 2.2, 121)
gains_curve = np.array([
    gain_and_erosion(p_a, p_c, r, gain_indep / 1e6)[1] for r in rs_curve
])
fig, ax = plt.subplots(figsize=(3.5, 2.6))
ax.plot(rs_curve, gains_curve, color="tab:blue", linewidth=1.4,
        label="Vienna A\\&C, $P_{\\mathrm{joint}}=R\\,p_A p_C$")
mark_rs = [1.00, 1.30, 1.47, 1.97]
mark_lbls = ["independence", "earth-grid 55 m (1.30)",
             "diff-site $R$ (1.47)", "diff-site CI upper (1.97)"]
mark_colors = ["0.4", "tab:green", "tab:orange", "tab:red"]
for r_val, lbl, col in zip(mark_rs, mark_lbls, mark_colors):
    g = gain_and_erosion(p_a, p_c, r_val, gain_indep / 1e6)[1]
    eros = gain_and_erosion(p_a, p_c, r_val, gain_indep / 1e6)[2]
    ax.scatter([r_val], [g], color=col, s=28, zorder=5, edgecolors="black",
               linewidths=0.4, label=f"$R$={r_val:.2f}: erosion={eros:.2f}\\%")
ax.axhline(gain_indep, color="0.6", linestyle=":", linewidth=0.7)
ax.set_xlabel(r"Mantel--Haenszel joint-outage rate ratio $R$")
ax.set_ylabel("Avail.\\ gain over single-best (ppm)")
ax.set_xlim(0.97, 2.03)
ax.set_ylim(bottom=min(gains_curve) * 0.985)
ax.grid(linestyle=":", linewidth=0.5, alpha=0.7)
ax.legend(loc="lower left", framealpha=0.9)
fig.tight_layout()
out_fig = C.FIGURES_DIR / "fig_element3_gain_vs_or.pdf"
fig.savefig(out_fig, bbox_inches="tight")
plt.close(fig)
print(f"[fig] {out_fig.name}")
print("Done.")
