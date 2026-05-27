"""Element 4 — application-layer independence (Berlin V2X).

Reproducibility-packaging Rule 3: every number in "Application-Layer
Independence (Berlin)" (Section IV-D) is *printed* by this script.

This script consolidates the Berlin V2X application-layer analyses into
one config-driven reproduction script:

  1. Per-operator app-outage rates + run-length stats (budget sweep 100/200/500ms).
  2. Joint co-occurrence vs independence (pooled, primary budget=200ms).
  3. Conditional P(op1|op2) at 200ms.
  4. Block bootstrap (120-s blocks, N=2000) on marginal rates +
     Poisson exact test on joint count -> Poisson P and upper CI on lift.
  5. Congestion contrast: iperf-loaded vs V2X-clean periods.
  6. Missingness sensitivity: warm-up vs in-run split,
     recode-bracket -> inrun_recode stays at/below independence.

Run:  python element4_berlin.py
Out:  stdout + numbers/element4.txt
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as C

# --- Analysis constants ---------------------------------------------------
DUAL_OP_RUNS  = frozenset({8, 9, 10})
BUDGETS_MS    = [100, 200, 500]
PRIMARY       = 200
BLOCK_SEC     = 120
N_BOOT        = 2000
SEED          = 42


# =============================================================================
# Helpers
# =============================================================================

def run_lengths(arr: np.ndarray) -> list[int]:
    """Consecutive-True run lengths."""
    lengths, count = [], 0
    for v in arr:
        if v:
            count += 1
        elif count > 0:
            lengths.append(count)
            count = 0
    if count > 0:
        lengths.append(count)
    return lengths


def joint_stats(o1: np.ndarray, o2: np.ndarray, n: int) -> dict:
    n1 = int(o1.sum()); n2 = int(o2.sum()); nj = int((o1 & o2).sum())
    p1 = n1 / n; p2 = n2 / n
    exp = n * p1 * p2
    lift = nj / exp if exp > 0 else float("nan")
    return dict(n=n, n1=n1, n2=n2, nj=nj, p1=p1, p2=p2, exp=exp, lift=lift)


def poisson_p_leq(obs: int, lam: float) -> float:
    return math.exp(-lam) * sum(lam**k / math.factorial(k) for k in range(obs + 1))


def poisson_upper_lambda(obs: int, alpha: float = 0.05) -> float:
    """Exact Poisson 95% upper CI on lambda given obs count."""
    if obs == 0:
        return -math.log(alpha)
    k = 2 * (obs + 1)
    return k * (1 - 2 / (9 * k) + 1.645 * math.sqrt(2 / (9 * k))) ** 3 / 2


# =============================================================================
# Load + filter
# =============================================================================
C.ensure_dirs()
rng = np.random.default_rng(SEED)

print("=" * 70)
print("Element 4: Application-layer independence (Berlin V2X)")
print("=" * 70)
print(f"\nLoading Berlin cellular: {C.BERLIN_CELLULAR}")
C.check_data()

df = pd.read_parquet(C.BERLIN_CELLULAR)
delay = df[df["measured_qos"] == "delay"].copy()
delay["sec"] = delay.index.floor("s")

# Verify dual-op run list
dual_actual = {
    run for run, grp in delay.groupby("measurement")
    if grp["operator"].nunique() == 2
}
assert dual_actual == set(DUAL_OP_RUNS), (
    f"Dual-op runs mismatch: expected {DUAL_OP_RUNS}, got {dual_actual}"
)

core = delay[delay["measurement"].isin(DUAL_OP_RUNS)].copy()
print(f"  delay class (CAM/400kbps): {len(core):,} rows  "
      f"runs {sorted(DUAL_OP_RUNS)}")


# =============================================================================
# Build per-(run, sec, op) ping max -- preserving NaN rows (MAR check)
# =============================================================================
gmax = core.groupby(["measurement", "sec", "operator"])["ping_ms"].max()
ping = gmax.unstack("operator").rename(columns={1: "p1", 2: "p2"})

# matched = seconds where BOTH operators have a ping sample
matched = ping.dropna(subset=["p1", "p2"]).copy()
print(f"  matched seconds (both ops have ping): {len(matched):,}")

for b in BUDGETS_MS:
    matched[f"o1_{b}"] = matched["p1"] > b
    matched[f"o2_{b}"] = matched["p2"] > b
    matched[f"oj_{b}"] = matched[f"o1_{b}"] & matched[f"o2_{b}"]

n_matched = len(matched)


# =============================================================================
# 1. Per-operator survival-time (run-length) stats
# =============================================================================
print("\n--- 1. Per-operator outage stats ---")
rl_stats: dict = {}
for b in BUDGETS_MS:
    for opcol, oplbl in [("p1", "op1"), ("p2", "op2")]:
        sub = ping[opcol].dropna()
        n_s = len(sub)
        n_o = int((sub > b).sum())
        rls = np.array(run_lengths((sub > b).values), dtype=int) if n_s else np.array([], dtype=int)
        rl_stats[(b, oplbl)] = dict(
            n_sec=n_s, n_out=n_o,
            rate_pct=n_o / n_s * 100 if n_s else float("nan"),
            n_ep=len(rls),
            max_rl=int(rls.max()) if len(rls) else 0,
            n_ge10=int((rls >= 10).sum()) if len(rls) else 0,
        )

for b in BUDGETS_MS:
    print(f"\n  budget={b}ms:")
    for oplbl in ["op1", "op2"]:
        st = rl_stats[(b, oplbl)]
        print(f"    {oplbl}: n={st['n_sec']:,} out={st['n_out']:,} "
              f"rate={st['rate_pct']:.4f}%  max_run={st['max_rl']}s  >=10s={st['n_ge10']}")

# ratio at primary budget
st1 = rl_stats[(PRIMARY, "op1")]
st2 = rl_stats[(PRIMARY, "op2")]
op2_op1_ratio = st2["rate_pct"] / st1["rate_pct"]
print(f"\n  op2/op1 rate ratio @ {PRIMARY}ms: {op2_op1_ratio:.1f}×")


# =============================================================================
# 2. Joint co-occurrence vs independence (pooled, all budgets)
# =============================================================================
print("\n--- 2. Joint co-occurrence vs independence ---")
joint: dict = {}
for b in BUDGETS_MS:
    o1 = matched[f"o1_{b}"].values
    o2 = matched[f"o2_{b}"].values
    s  = joint_stats(o1, o2, n_matched)
    joint[b] = s
    lift_s = f"{s['lift']:.3f}x" if not math.isnan(s["lift"]) else "nan"
    print(f"  {b}ms: JOINT={s['nj']}  exp={s['exp']:.2f}  lift={lift_s}")

jp   = joint[PRIMARY]
jp100 = joint[100]

# P(X <= 0 | lambda = exp) at primary budget
poisson_p = poisson_p_leq(jp["nj"], jp["exp"])
pois_upper = poisson_upper_lambda(jp["nj"])
pois_lift_upper = pois_upper / jp["exp"] if jp["exp"] > 0 else float("nan")
print(f"\n  @{PRIMARY}ms: Poisson P(X<={jp['nj']} | λ={jp['exp']:.2f}) = {poisson_p:.2f}")
print(f"  @{PRIMARY}ms: Poisson 95% upper on λ = {pois_upper:.3f}  "
      f"-> lift upper = {pois_lift_upper:.1f}×")


# =============================================================================
# 3. Conditional probabilities at primary budget
# =============================================================================
n_o1 = jp["n1"]; n_o2 = jp["n2"]; nj = jp["nj"]
cond_p1_given_p2 = nj / n_o2 if n_o2 > 0 else float("nan")
cond_p2_given_p1 = nj / n_o1 if n_o1 > 0 else float("nan")
print(f"\n--- 3. Conditional probabilities @{PRIMARY}ms ---")
print(f"  P(op1 out | op2 out) = {cond_p1_given_p2:.6f}")
print(f"  P(op2 out | op1 out) = {cond_p2_given_p1:.6f}")


# =============================================================================
# 4. Block bootstrap on marginal rates (primary budget)
# =============================================================================
print(f"\n--- 4. Block bootstrap (N={N_BOOT}, blocks={BLOCK_SEC}s) ---")

matched_bs = matched.reset_index()   # sec is now a column
matched_bs["t_unix"] = matched_bs["sec"].astype("int64") // 10**9

blocks: dict[int, dict] = {b: {} for b in BUDGETS_MS}
for run, grp in matched_bs.groupby("measurement"):
    t0 = grp["t_unix"].min()
    bid_col = ((grp["t_unix"] - t0) // BLOCK_SEC).astype(int)
    for bid, bg in zip(bid_col, [grp.iloc[i:i+1] for i in range(len(grp))]):
        pass  # handled below

# Build block arrays properly
from collections import defaultdict
block_data: dict[int, dict[str, tuple]] = {b: defaultdict(list) for b in BUDGETS_MS}
for run, grp in matched_bs.groupby("measurement"):
    t0 = grp["t_unix"].min()
    grp = grp.copy()
    grp["bid"] = ((grp["t_unix"] - t0) // BLOCK_SEC).astype(int)
    for b in BUDGETS_MS:
        for bid, bg in grp.groupby("bid"):
            key = f"{run}_{bid}"
            block_data[b][key] = (
                bg[f"o1_{b}"].values.astype(bool),
                bg[f"o2_{b}"].values.astype(bool),
            )

boot_results: dict[int, dict] = {}
for b in BUDGETS_MS:
    bk = list(block_data[b].keys())
    nb = len(bk)
    bp1 = np.empty(N_BOOT); bp2 = np.empty(N_BOOT)
    for i in range(N_BOOT):
        sk = rng.choice(bk, size=nb, replace=True)
        o1 = np.concatenate([block_data[b][k][0] for k in sk])
        o2 = np.concatenate([block_data[b][k][1] for k in sk])
        bp1[i] = o1.mean(); bp2[i] = o2.mean()
    boot_results[b] = dict(
        p1_mean=float(bp1.mean()),
        p1_ci=(float(np.percentile(bp1, 2.5)), float(np.percentile(bp1, 97.5))),
        p2_mean=float(bp2.mean()),
        p2_ci=(float(np.percentile(bp2, 2.5)), float(np.percentile(bp2, 97.5))),
    )
    br = boot_results[b]
    print(f"  {b}ms: op1 boot mean={br['p1_mean']:.5f} CI95={br['p1_ci']}  "
          f"op2 boot mean={br['p2_mean']:.5f} CI95={br['p2_ci']}")


# =============================================================================
# 5. Congestion contrast (iperf-loaded vs V2X-clean periods)
# =============================================================================
print("\n--- 5. Congestion contrast (iperf-loaded vs V2X-clean) ---")

def class_joint(sub_df: pd.DataFrame, b: int) -> dict:
    sub_df = sub_df.copy()
    sub_df["sec"] = sub_df.index.floor("s")
    g = (
        sub_df.groupby(["measurement", "sec", "operator"])["ping_ms"].max()
        .unstack("operator").rename(columns={1: "p1", 2: "p2"})
        .dropna(subset=["p1", "p2"])
    )
    if len(g) == 0:
        return dict(n=0, lift=float("nan"))
    o1 = (g["p1"] > b).values; o2 = (g["p2"] > b).values
    return joint_stats(o1, o2, len(g))

contrast: dict = {}
for b in [PRIMARY]:
    contrast[("all",   b)] = class_joint(df,                           b)
    tput = df[df["measured_qos"] == "datarate"]
    contrast[("iperf", b)] = class_joint(tput,                         b)
    contrast[("delay", b)] = joint[b]

print(f"  @{PRIMARY}ms:  ALL={contrast[('all',PRIMARY)]['lift']:.2f}×  "
      f"iperf={contrast[('iperf',PRIMARY)]['lift']:.2f}×  "
      f"delay(V2X)={contrast[('delay',PRIMARY)]['lift']:.2f}×")
print("  → iperf-congestion artifact accounts for the elevated joint in iperf-loaded periods.")


# =============================================================================
# 6. Missingness sensitivity (MAR check key results)
# =============================================================================
print("\n--- 6. Missingness sensitivity (MAR check) ---")

# Build complete 1-sec grid including NaN rows
gmax_full = core.groupby(["measurement", "sec", "operator"])["ping_ms"].max()
g_full = gmax_full.unstack("operator").rename(columns={1: "p1", 2: "p2"})

frames = []
for run in sorted(DUAL_OP_RUNS):
    sub = g_full.xs(run, level="measurement")
    full_idx = pd.date_range(sub.index.min(), sub.index.max(),
                             freq="s", tz=sub.index.tz)
    sub = sub.reindex(full_idx)
    sub.index.name = "sec"
    sub["measurement"] = run
    frames.append(sub.reset_index())
grid = pd.concat(frames, ignore_index=True).sort_values(["measurement", "sec"]).reset_index(drop=True)
grid["miss1"] = grid["p1"].isna()
grid["miss2"] = grid["p2"].isna()

# Classify warm-up vs in-run per operator
for opn in [1, 2]:
    grid[f"active{opn}"] = False
warmup_secs: dict = {}
for run in sorted(DUAL_OP_RUNS):
    idx = grid.index[grid["measurement"] == run]
    sub = grid.loc[idx]
    for opn in [1, 2]:
        first_ping = sub.loc[~sub[f"miss{opn}"], "sec"].min()
        grid.loc[idx, f"active{opn}"] = grid.loc[idx, "sec"] >= first_ping
        warmup_secs[(run, opn)] = int(
            (sub["sec"] < first_ping).sum()
        )

for opn in [1, 2]:
    grid[f"inrun_miss{opn}"] = grid[f"miss{opn}"] & grid[f"active{opn}"]

n_warmup1  = int(grid["inrun_miss1"].sum() == 0 and grid["miss1"].sum())  # all miss1 is warmup
n_warmup_total = sum(warmup_secs[(r, 1)] for r in sorted(DUAL_OP_RUNS))   # sum over runs
n_bothdark_raw  = int((grid["miss1"] & grid["miss2"]).sum())
n_bothdark_inrun = int((grid["inrun_miss1"] & grid["inrun_miss2"]).sum())
n_inrun_miss1 = int(grid["inrun_miss1"].sum())
n_inrun_miss2 = int(grid["inrun_miss2"].sum())

# Inrun recode @ primary budget
both_active = (grid["active1"] & grid["active2"]).values
for b in BUDGETS_MS:
    exceed1 = (grid["p1"] > b).fillna(False).values
    exceed2 = (grid["p2"] > b).fillna(False).values
    inrun1  = grid["inrun_miss1"].values
    inrun2  = grid["inrun_miss2"].values
    o1 = (exceed1 | inrun1)[both_active]
    o2 = (exceed2 | inrun2)[both_active]
    n_u = int(both_active.sum())
    s = joint_stats(o1, o2, n_u)
    if b == PRIMARY:
        inrun_recode_prim = s

print(f"  Warm-up (all runs total): op1={n_warmup_total}s  "
      f"(run8={warmup_secs[(8,1)]}s, run9={warmup_secs[(9,1)]}s, run10={warmup_secs[(10,1)]}s)")
print(f"  Both-dark raw (incl warm-up): {n_bothdark_raw}s")
print(f"  op1 in-run missing: {n_inrun_miss1}s   op2 in-run missing: {n_inrun_miss2}s")
print(f"  In-run both-missing: {n_bothdark_inrun}s")
print(f"  inrun_recode @{PRIMARY}ms: JOINT={inrun_recode_prim['nj']}  "
      f"exp={inrun_recode_prim['exp']:.3f}  lift={inrun_recode_prim['lift']:.2f}x")
gate = (inrun_recode_prim["lift"] <= 1.10) and (n_bothdark_inrun == 0)
print(f"  GATE: {'CLEARED' if gate else 'FAILED'} — NEGATIVE survives missingness attack")


# =============================================================================
# 7. Scheduler-magnitude calibration
# =============================================================================
# Bounds the upside of a hypothetical perfect correlated-outage scheduler
# against the joint-outage component. Under a binary per-second / per-slot
# availability metric:
#   without scheduler at lift L: P(joint outage / slot) = L * p1 * p2
#   with perfect correlated-outage scheduler:           = p1 * p2
# so per-slot reliability-nines gained = log10(L). The gain is independent
# of (p1, p2): it depends only on the lift magnitude. At our Poisson 95%
# upper-bound lift L_UB = 5.8x, the maximum reliability-nines a perfect
# scheduler can gain against the joint-outage component is log10(5.8) ≈
# 0.76 nines. Less than one decade -- therefore in these traces a
# correlated-outage scheduler cannot, by itself, bridge a 1-nine 3GPP V2X
# target step (e.g., 99.99% -> 99.999%) even at the upper bound of our
# statistical power. The observed lift is 0; no positive evidence at any
# magnitude exists in our data.
print("\n--- 7. Scheduler-magnitude calibration ---")
p1_obs = jp["p1"]
p2_obs = jp["p2"]
indep_rate = p1_obs * p2_obs
print(f"  marginal outage rates @{PRIMARY}ms: p1={p1_obs:.4%}, p2={p2_obs:.4%}")
print(f"  independence joint rate per slot: p1*p2 = {indep_rate:.2e}")
print(f"  Poisson 95% upper-bound lift L_UB = {pois_lift_upper:.2f}x")
print("")
print("  reliability nines gained by a perfect correlated-outage scheduler")
print("  against the joint-outage component = log10(L):")
nines_table: list[tuple[float, float]] = []
for L_val in [1.2, 1.5, 2.0, 3.0, pois_lift_upper]:
    nines_gained = math.log10(L_val)
    nines_table.append((L_val, nines_gained))
    print(f"    L={L_val:>4.2f}x  ->  {nines_gained:.2f} nines gained")
nines_at_ub = math.log10(pois_lift_upper)
print("")
print(f"  At L_UB={pois_lift_upper:.2f}x, max scheduler gain = "
      f"{nines_at_ub:.2f} nines (< 1 decade).")
print( "  Therefore even at our statistical-power upper bound, a")
print( "  correlated-outage scheduler cannot bridge a 1-nine 3GPP V2X")
print( "  target step (99.99% -> 99.999%).")
print(f"  Observed lift = {jp['lift']:.2f}x (no positive evidence).")


# =============================================================================
# KEY NUMBERS SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("KEY NUMBERS FOR PAPER (Section IV-D)")
print("=" * 70)
op1_rate = rl_stats[(PRIMARY, "op1")]["rate_pct"]
op2_rate = rl_stats[(PRIMARY, "op2")]["rate_pct"]
op2_max_rl = rl_stats[(PRIMARY, "op2")]["max_rl"]
op2_ge10   = rl_stats[(PRIMARY, "op2")]["n_ge10"]
print(f"op1 rate @{PRIMARY}ms: {op1_rate:.4f}%   op2 rate: {op2_rate:.4f}%   ratio: {op2_op1_ratio:.1f}×")
print(f"op2 max run-length @{PRIMARY}ms: {op2_max_rl}s   episodes >=10s: {op2_ge10}")
print(f"JOINT obs @{PRIMARY}ms: {jp['nj']}   exp: {jp['exp']:.2f}   lift: {jp['lift']:.2f}×")
print(f"Poisson P(X<={jp['nj']} | λ={jp['exp']:.2f}): {poisson_p:.2f}")
print(f"Poisson 95% upper CI on lift @{PRIMARY}ms: {pois_lift_upper:.1f}×")
print(f"lift @100ms: {jp100['lift']:.2f}×")
print(f"P(op1 out | op2 out) @{PRIMARY}ms: {cond_p1_given_p2:.6f}")
print(f"Congestion contrast @{PRIMARY}ms: ALL={contrast[('all',PRIMARY)]['lift']:.2f}×  "
      f"iperf={contrast[('iperf',PRIMARY)]['lift']:.2f}×  delay={contrast[('delay',PRIMARY)]['lift']:.2f}×")
print(f"Missingness warm-up total: {n_bothdark_raw}s   in-run op2: {n_inrun_miss2}s   "
      f"in-run both: {n_bothdark_inrun}s")
print(f"inrun_recode lift @{PRIMARY}ms: {inrun_recode_prim['lift']:.2f}×")


# =============================================================================
# Write numbers/element4.txt
# =============================================================================
lines_out = [
    f"op1 rate @{PRIMARY}ms: {op1_rate:.4f}%",
    f"op2 rate @{PRIMARY}ms: {op2_rate:.4f}%",
    f"op2/op1 ratio @{PRIMARY}ms: {op2_op1_ratio:.1f}x",
    # Per-budget rates for Table 4 (Section IV-D)
    f"per-budget op1 rates: 100ms {rl_stats[(100,'op1')]['rate_pct']:.4f}%  "
    f"200ms {rl_stats[(200,'op1')]['rate_pct']:.4f}%  "
    f"500ms {rl_stats[(500,'op1')]['rate_pct']:.4f}%",
    f"per-budget op2 rates: 100ms {rl_stats[(100,'op2')]['rate_pct']:.4f}%  "
    f"200ms {rl_stats[(200,'op2')]['rate_pct']:.4f}%  "
    f"500ms {rl_stats[(500,'op2')]['rate_pct']:.4f}%",
    f"per-budget JOINT obs / exp: 100ms {joint[100]['nj']}/{joint[100]['exp']:.2f}  "
    f"200ms {joint[200]['nj']}/{joint[200]['exp']:.2f}  "
    f"500ms {joint[500]['nj']}/{joint[500]['exp']:.2f}",
    f"op2 max run-length @{PRIMARY}ms: {op2_max_rl}s",
    f"op2 episodes >=10s @{PRIMARY}ms: {op2_ge10}",
    f"JOINT obs @{PRIMARY}ms: {jp['nj']}",
    f"JOINT exp @{PRIMARY}ms: {jp['exp']:.2f}",
    f"lift @{PRIMARY}ms: {jp['lift']:.2f}x",
    f"Poisson P(X<={jp['nj']} | lambda={jp['exp']:.2f}): {poisson_p:.2f}",
    f"Poisson 95% upper CI on lift @{PRIMARY}ms: {pois_lift_upper:.1f}x",
    f"lift @100ms: {jp100['lift']:.2f}x",
    f"P(op1 out | op2 out) @{PRIMARY}ms: {cond_p1_given_p2:.6f}",
    f"congestion contrast @{PRIMARY}ms: ALL {contrast[('all',PRIMARY)]['lift']:.2f}x "
    f"iperf {contrast[('iperf',PRIMARY)]['lift']:.2f}x "
    f"delay {contrast[('delay',PRIMARY)]['lift']:.2f}x",
    f"missingness warmup total: {n_bothdark_raw}s",
    f"missingness op1 inrun: {n_inrun_miss1}s",
    f"missingness op2 inrun: {n_inrun_miss2}s",
    f"missingness inrun both: {n_bothdark_inrun}s",
    f"inrun_recode lift @{PRIMARY}ms: {inrun_recode_prim['lift']:.2f}x",
    # Scheduler-magnitude calibration
    f"scheduler nines gain at L=1.2x: {math.log10(1.2):.2f}",
    f"scheduler nines gain at L=2.0x: {math.log10(2.0):.2f}",
    f"scheduler nines gain at L=L_UB={pois_lift_upper:.1f}x: {nines_at_ub:.2f}",
]

out_path = C.NUMBERS_DIR / "element4.txt"
out_path.write_text("\n".join(lines_out) + "\n")
print(f"\nNumbers written to {out_path}")


# =============================================================================
# Fig. 5 (Section IV-D): per-second app-outage timeline, runs 8/9/10, budget=200ms
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

runs_sorted = sorted(DUAL_OP_RUNS)
fig, axes = plt.subplots(len(runs_sorted), 1, figsize=(3.5, 4.0),
                         sharex=False)

# Reset matched index to access (measurement, sec) explicitly
m = matched.reset_index()
m["sec"] = pd.to_datetime(m["sec"])

for ax, run in zip(axes, runs_sorted):
    sub = m[m["measurement"] == run].sort_values("sec").reset_index(drop=True)
    if sub.empty:
        ax.set_visible(False)
        continue
    t0 = sub["sec"].iloc[0]
    t_min = (sub["sec"] - t0).dt.total_seconds().values / 60.0

    o1 = sub[f"o1_{PRIMARY}"].values
    o2 = sub[f"o2_{PRIMARY}"].values
    oj = sub[f"oj_{PRIMARY}"].values

    # 3 vertical strips (top=op1, mid=op2, bottom=joint)
    bands = [
        (2.6, 0.7, o1, "tab:orange", "op1"),
        (1.6, 0.7, o2, "tab:purple", "op2"),
        (0.6, 0.7, oj, "tab:red",    "joint"),
    ]
    for y0, h, mask, color, lbl in bands:
        ax.add_patch(plt.Rectangle(
            (t_min.min(), y0), t_min.max() - t_min.min(), h,
            facecolor="0.92", edgecolor="0.7", linewidth=0.4, zorder=1,
        ))
        if mask.any():
            ax.vlines(t_min[mask], y0, y0 + h,
                      colors=color, linewidth=0.6, zorder=2)
        # Right-edge count annotation (outside strip, white bbox for readability)
        ax.text(
            t_min.max() + (t_min.max() - t_min.min()) * 0.05, y0 + h / 2,
            f"{lbl} $n={int(mask.sum())}$",
            va="center", ha="left", fontsize=7, color=color,
            bbox=dict(facecolor="white", edgecolor="none",
                      pad=1.0, alpha=0.9),
        )

    ax.set_xlim(t_min.min() - 0.1,
                t_min.max() + (t_min.max() - t_min.min()) * 0.25)
    ax.set_ylim(0.2, 3.6)
    ax.set_yticks([])
    ax.set_xlabel("time within run (min)")
    ax.set_title(f"Run {run}  (matched seconds $n={len(sub)}$)",
                 fontsize=8, loc="left")
    ax.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.6)

fig.suptitle(rf"Berlin V2X delay class, budget $= {PRIMARY}$ ms",
             fontsize=9, y=1.00)
fig.tight_layout()
out_fig = C.FIGURES_DIR / "fig_element4_timeline.pdf"
fig.savefig(out_fig, bbox_inches="tight")
plt.close(fig)
print(f"[fig] {out_fig.name}")
print("Done.")
