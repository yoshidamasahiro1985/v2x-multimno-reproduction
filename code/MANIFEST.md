# Manifest — claim → producing script

Every figure, table, and **number cited in the prose** must map to exactly
one script + command that regenerates it (reproducibility-packaging Rule 3).
Filled incrementally in Phase 3 as each result enters the manuscript.

## Figures

| Figure | Caption (short) | Script | Command |
|---|---|---|---|
| 1 / fig_element1_radio_scatter.pdf | §IV-A: 2-panel hex-bin (a) Vienna LTE SINR vs DL throughput, deg-3 fit overlay R²=0.13; (b) DoNext 5G-NSA SS-RSRP vs DL throughput, \|r\|=0.21 | element1_radio_throughput.py | `python3 element1_radio_throughput.py` |
| 2 / fig_element2_forest_grid.pdf | §IV-B: (a) Forest plot of MH ratios (A&C all-geo/same/diff + A&B/B&C controls) w/ 95% CIs; (b) earth-grid descent 220→28m | element2_dependence.py | `python3 element2_dependence.py` |
| 3 / fig_element2_mechanism.pdf | §IV-B: Earth-grid descent overlay of DEF_A / RSRP-only / SINR-only A&C joint-outage R; three curves tracking → coverage-driven mechanism, not RAN-scheduler interference | element2_dependence.py | `python3 element2_dependence.py` |
| 4 / fig_element3_gain_vs_or.pdf | §IV-C: Vienna A&C availability gain (ppm) vs Mantel-Haenszel rate ratio R, with markers at R=1.00/1.30/1.47/1.97 + erosion% labels | element3_redundancy.py | `python3 element3_redundancy.py` |
| 5 / fig_element4_timeline.pdf | §IV-D: Per-second outage timeline for Berlin V2X runs 8/9/10 @budget=200ms; op1/op2/joint strips; joint = 0 across all 3 runs | element4_berlin.py | `python3 element4_berlin.py` |

## Tables

| Table | Content | Script | Command |
|---|---|---|---|
| 1 | §III: Datasets overview (roles); sample counts trace to dataset_stats.py | dataset_stats.py | `python3 dataset_stats.py` |
| 2 | §IV-B: sub-1 GHz fallback fractions (A/B/C × overall/JO slots) — mechanism evidence | element2_dependence.py | `python3 element2_dependence.py` |
| 3 | §IV-C: rate-ratio R family gain table (independence/1.30/1.47/1.97 → gain ppm + erosion%); Berlin app-layer row (indep-predicted + observed) | element3_redundancy.py | `python3 element3_redundancy.py` |
| 4 | §IV-D: Per-budget joint app-outage stats (100/200/500ms × matched-op1-rate/op2-rate/JOINT obs/exp/lift) | element4_berlin.py | `python3 element4_berlin.py` |

## In-text numbers

Each number stated in the prose, the script that **prints** it (never
hand-transcribed), and the `numbers/` log line it appears in.

| § | Number (claim) | Value | Script | numbers/ log |
|---|---|---|---|---|
| III | Vienna LTE samples (A/B/C %) | 1,183,683 (49.7/24.7/25.6) | dataset_stats.py | numbers/dataset_stats.txt |
| III | Vienna sessions / dates / active h | 12 / 8 / 48.5 | dataset_stats.py | numbers/dataset_stats.txt |
| III | Vienna A&C co-presence dates | 6 | dataset_stats.py | numbers/dataset_stats.txt |
| III | Vienna cadence A / B,C | 3.46 / 2.00 Hz | dataset_stats.py | numbers/dataset_stats.txt |
| III | Vienna RSRP / SINR present | 100% / 99.9% | dataset_stats.py | numbers/dataset_stats.txt |
| III | Berlin samples (op1/op2 %) | 207,434 (49.9/50.1) | dataset_stats.py | numbers/dataset_stats.txt |
| III | Berlin delay / datarate class | 77,848 / 129,586 | dataset_stats.py | numbers/dataset_stats.txt |
| III | Berlin ping_ms null rate | 28.9% | dataset_stats.py | numbers/dataset_stats.txt |
| III | Berlin V2X dual-op runs 8–10 samples | 43,200 (21,600/op) | dataset_stats.py | numbers/dataset_stats.txt |
| III | Berlin convoy op1↔op2 sep (median/p90) | 3.1 / 10.9 m | dataset_stats.py | numbers/dataset_stats.txt |
| III | DoNext mobile samples (A/B/C) | 2,382,602 (48.2/51.8/13 rows) | dataset_stats.py | numbers/dataset_stats.txt |
| III | DoNext RAT split 5G-NSA / 4G | 91.5% / 5.3% | dataset_stats.py | numbers/dataset_stats.txt |
| R1 | Vienna LTE radio→throughput R² (deg-3, n=1,155,265) | 0.13 | element1_radio_throughput.py | numbers/element1.txt |
| R1 | DoNext clean downlink full-buffer tests / median | 104,742 / 94 Mbps | element1_radio_throughput.py | numbers/element1.txt |
| R1 | DoNext 5G-NSA R² (SINR / multivar-4 / 7-feat best) | 0.004 / 0.037 / 0.10 | element1_radio_throughput.py | numbers/element1.txt |
| R1 | DoNext 4G R² (best, multivar-4) | 0.18 | element1_radio_throughput.py | numbers/element1.txt |
| R1 | DoNext 5G-NSA max single |Pearson corr| (ss_rsrp) | 0.21 | element1_radio_throughput.py | numbers/element1.txt |
| R2 | A&C geolocatable co-present slots / JO (all/geo) | 70,531 / 79 / 28 | element2_dependence.py | numbers/element2.txt |
| R2 | A&C all-geo MH R + CI | 1.90 [1.49, 2.65] | element2_dependence.py | numbers/element2.txt |
| R2 | A&C same-site R (10 JO, 3 cell-pairs) | ~4.0 | element2_dependence.py | numbers/element2.txt |
| R2 | A&C diff-site residual R + CI (18 JO) | 1.47 [1.15, 1.98] | element2_dependence.py | numbers/element2.txt |
| R2 | diff-site Bonferroni (p1 / ×3 / 98.3% LB) | 0.007 / 0.021 / 1.11 | element2_dependence.py | numbers/element2.txt |
| R2 | neg controls A&B / B&C R + CI | 1.07 [0.84,1.12] / 0.85 [0.66,1.36] | element2_dependence.py | numbers/element2.txt |
| R2 | earth-grid descent (220m→28m), Kendall τ | 1.97→1.25, +1 | element2_dependence.py | numbers/element2.txt |
| R2 | earth-grid DEF_A @ 110 m R+CI / n_JO | 1.64 [1.36, 1.83] / 79 | element2_dependence.py | numbers/element2.txt |
| R2 | earth-grid RSRP-only @ 110 m R+CI / n_JO | 1.81 [1.33, 2.18] / 42 | element2_dependence.py | numbers/element2.txt |
| R2 | earth-grid SINR-only @ 110 m R+CI / n_JO | 1.70 [1.09, 2.55] / 13 | element2_dependence.py | numbers/element2.txt |
| R2 | coord validation: same-site max / diff-site min / margin / safety | 40 / 309 / 259 m / 6.4× | element2_dependence.py | numbers/element2.txt |
| R2 | sub-1 GHz @JO vs overall (A / C / B) | 34.2/12.7/0% vs 14.2/3.4/0% | element2_dependence.py | numbers/element2.txt |
| R2 / Tab II | sub-1 GHz Lift (at-JO / overall, A / C / B) | 2.41× / 3.77× / N/A | element2_dependence.py | numbers/element2.txt |
| R3 | Vienna A / C marginal outage rates (DEF_A) | 0.85% / 1.71% | element3_redundancy.py | numbers/element3.txt |
| R3 | Vienna single-best availability | 99.15% | element3_redundancy.py | numbers/element3.txt |
| R3 | Vienna A&C gain under independence | 8,351 ppm | element3_redundancy.py | numbers/element3.txt |
| R3 | Vienna A&C gain @R=1.30 + erosion | 8,308 ppm / 0.52% | element3_redundancy.py | numbers/element3.txt |
| R3 | gain erosion @R=1.97 (CI upper) | 1.68% | element3_redundancy.py | numbers/element3.txt |
| R3 | retains at R=1.97 CI upper | 98.3% | element3_redundancy.py | numbers/element3.txt |
| R3 | break-even Vienna (R=1.30) | 0.83% | element3_redundancy.py | numbers/element3.txt |
| R3 | Berlin gain (independence-predicted) | 2290.2 ppm | element3_redundancy.py | numbers/element3.txt |
| R3 | Berlin gain (JOINT=0, observed) | 2338.6 ppm | element3_redundancy.py | numbers/element3.txt |
| R3 | Berlin break-even (observed) | 0.23% | element3_redundancy.py | numbers/element3.txt |
| R4 | op1 / op2 outage rate @200ms | 0.23% / 2.07% | element4_berlin.py | numbers/element4.txt |
| R4 | op2/op1 ratio @200ms | ~9× | element4_berlin.py | numbers/element4.txt |
| R4 | Berlin JOINT @≥200ms (exp / Poisson P) | 0 / 0.52 / 0.60 | element4_berlin.py | numbers/element4.txt |
| R4 | lift @100ms | 0.53× | element4_berlin.py | numbers/element4.txt |
| R4 | P(op1 out \| op2 out) @200ms | 0.00 | element4_berlin.py | numbers/element4.txt |
| R4 | op2 max run-length / episodes ≥10s @200ms | 33s / 6 | element4_berlin.py | numbers/element4.txt |
| R4 | Poisson 95% upper CI on lift | 5.8× | element4_berlin.py | numbers/element4.txt |
| R4 | congestion contrast (all/iperf/delay) | 1.45/1.14/0.00× | element4_berlin.py | numbers/element4.txt |
| R4 | missingness warm-up total / inrun op2 / inrun both | 110s / 116s / 0 | element4_berlin.py | numbers/element4.txt |
| R4 | inrun_recode lift @200ms | 0.00× | element4_berlin.py | numbers/element4.txt |
| R4 | per-budget op1 rates (100/200/500ms) | 0.61/0.23/0.06% | element4_berlin.py | numbers/element4.txt |
| R4 | per-budget op2 rates (100/200/500ms) | 2.93/2.07/1.62% | element4_berlin.py | numbers/element4.txt |
| R4 | per-budget JOINT obs/exp (100/200/500ms) | 1/1.88, 0/0.52, 0/0.10 | element4_berlin.py | numbers/element4.txt |
| R3 | Vienna gain @R=1.47 + erosion | 8,283 ppm / 0.82% | element3_redundancy.py | numbers/element3.txt |
| R2 | Cochran-Q diff-site homogeneity (K/Q/df/p/I²) | 8 / 13.6 / 7 / 0.058 / 49% | element2_dependence.py | numbers/element2.txt |
| R2 | Cochran-Q same-site (K/Q/df/p/I²) | 3 / 29.0 / 2 / <0.001 / 93% | element2_dependence.py | numbers/element2.txt |
| R2 | Cochran-Q all-geo (K/Q/df/p/I²) | 11 / 51.4 / 10 / <0.001 / 81% | element2_dependence.py | numbers/element2.txt |
| §VI | scheduler nines gain at L=L_UB=5.8× | 0.76 nines | element4_berlin.py | numbers/element4.txt |
| §VI | scheduler nines gain at L=1.5× / 2.0× (illustrative) | 0.18 / 0.30 nines | element4_berlin.py | numbers/element4.txt |

> The values above are the **recorded targets** from the analysis (SESSION.md
> / session logs). In Phase 3 each must be re-emitted by its script and the
> printed value checked against the manuscript before the row is final.
