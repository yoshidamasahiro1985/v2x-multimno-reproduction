# Reproduction package — multi-MNO V2X reliability measurement study

This package regenerates every figure, table, and in-text number in the
manuscript *"A Cross-Layer Measurement Study of Multi-Operator Cellular
Redundancy for V2X Reliability"*. See `MANIFEST.md` for the claim →
script mapping.

## Requirements

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Python 3.14; pinned `pandas`, `numpy`, `pyarrow`, `matplotlib` (see
`requirements.txt`). All scripts use a fixed seed (42) for deterministic
output. The figures are emitted via matplotlib's `pgf` backend with
`pdflatex`, so a working LaTeX installation is required for the figure
step (`numbers/` logs are produced regardless).

## Data (not redistributed — download from the canonical sources)

This package does **not** include raw data. Download the three datasets
from their original repositories (cited in the paper) and point
`V2X_DATA_DIR` at the directory containing them:

```
export V2X_DATA_DIR=/path/to/datasets
```

| Dataset | Used for | License / source |
|---|---|---|
| Vienna 4G/5G smartphone | elements 1, 2, 3 | CC BY 4.0 — DOI: [10.48550/arXiv.2603.02638](https://doi.org/10.48550/arXiv.2603.02638) |
| Berlin V2X | element 4 | IEEE DataPort licence — DOI: [10.1109/VTC2023-Spring57618.2023.10200750](https://doi.org/10.1109/VTC2023-Spring57618.2023.10200750) |
| DoNext | element 1 | Open access (consult the dataset record for the precise licence) — DOI: [10.1109/TMLCN.2025.3564239](https://doi.org/10.1109/TMLCN.2025.3564239) |

Expected layout under `V2X_DATA_DIR` is documented at the top of
`config.py`. Run `python3 -c "import config; config.check_data()"` to
verify placement.

## Run

```
bash make_all.sh          # regenerate everything
# or individually:
python3 dataset_stats.py
python3 element1_radio_throughput.py
python3 element2_dependence.py
python3 element3_redundancy.py
python3 element4_berlin.py
```

Outputs:

- Figures → `../figures/` (one directory above this package; populated
  with `fig_element*.pdf`)
- In-text numbers → `numbers/` (one printed log per element)

## Layout

```
config.py            # paths (via V2X_DATA_DIR), seed, outage thresholds
requirements.txt     # pinned dependency versions
make_all.sh          # driver: runs every script in order
dataset_stats.py     # §III dataset-characterization numbers
element1_*.py        # §IV-A radio quality → throughput regressions
element2_*.py        # §IV-B radio-layer multi-MNO outage dependence
element3_*.py        # §IV-C redundancy utility under measured dependence
element4_*.py        # §IV-D Berlin V2X application-layer independence
numbers/             # printed in-text statistics (one file per element)
MANIFEST.md          # figure / table / number → producing script map
../figures/          # output target for the manuscript figures
```

## Outage definition

The radio-layer outage definition (`DEF_A`) is encoded at the top of
`config.py` as the rule `RSRP < -110 dBm or SINR < -6 dB` and is applied
uniformly across elements 2 and 3. Adjust these floors there if you wish
to sweep the definition.
