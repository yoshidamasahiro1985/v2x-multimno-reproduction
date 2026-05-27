"""Central configuration for the reproduction package.

All data paths and output locations live here (reproducibility-packaging
Rule 3: no hardcoded absolute paths in element scripts). Element scripts
import from this module only.

Data location
-------------
Set the environment variable ``V2X_DATA_DIR`` to the directory that
contains the three datasets (see README.md for download instructions and
DOIs). If unset, it defaults to ``<repo>/data/external`` (the private
research repo's symlink target). The datasets themselves are NOT
redistributed with this package.

Expected layout under ``V2X_DATA_DIR``::

    vienna-4g5g/Vienna-full-4G-5G-dataset-v1.0.3/phone/phone_data_lte.parquet
    berlin-v2x/data/cellular_dataframe.parquet
    berlin-v2x/data/sidelink_dataframe.parquet
    donext/mobile/cell_data.csv
    donext/mobile/iperf_data.csv
    donext/mobile/latency_data.csv
    donext/mobile/neighboring_data.csv
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

# --- Reproducibility ---------------------------------------------------
SEED = 42


def rng() -> np.random.Generator:
    """Fresh seeded generator. Call once per script for determinism."""
    return np.random.default_rng(SEED)


# --- Paths -------------------------------------------------------------
_THIS = Path(__file__).resolve()
PKG_DIR = _THIS.parent                       # paper/code/
PAPER_DIR = PKG_DIR.parent                   # paper/
REPO_ROOT = PAPER_DIR.parent                 # repo root

DATA_DIR = Path(os.environ.get("V2X_DATA_DIR", REPO_ROOT / "data" / "external"))

# Dataset files
VIENNA_LTE = (
    DATA_DIR / "vienna-4g5g" / "Vienna-full-4G-5G-dataset-v1.0.3"
    / "phone" / "phone_data_lte.parquet"
)
VIENNA_CELLINFO = (
    DATA_DIR / "vienna-4g5g" / "Vienna-full-4G-5G-dataset-v1.0.3"
    / "estimated_cell_info" / "cell_info_final_lte.parquet"
)
BERLIN_CELLULAR = DATA_DIR / "berlin-v2x" / "data" / "cellular_dataframe.parquet"
BERLIN_SIDELINK = DATA_DIR / "berlin-v2x" / "data" / "sidelink_dataframe.parquet"
DONEXT_DIR = DATA_DIR / "donext" / "mobile"
DONEXT_CELL = DONEXT_DIR / "cell_data.csv"
DONEXT_IPERF = DONEXT_DIR / "iperf_data.csv"
DONEXT_LATENCY = DONEXT_DIR / "latency_data.csv"
DONEXT_NEIGHBORING = DONEXT_DIR / "neighboring_data.csv"

# Outputs (figures shared with the LaTeX manuscript; numbers = in-text stats)
FIGURES_DIR = PAPER_DIR / "figures"
NUMBERS_DIR = PKG_DIR / "numbers"

# --- Outage definition (DEF_A) -----------------------------------------
RSRP_FLOOR_DBM = -110.0
SINR_FLOOR_DB = -6.0


def ensure_dirs() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    NUMBERS_DIR.mkdir(parents=True, exist_ok=True)


def check_data() -> None:
    """Fail fast with a clear message if datasets are not where expected."""
    missing = [
        str(p)
        for p in (VIENNA_LTE, BERLIN_CELLULAR, DONEXT_CELL)
        if not p.exists()
    ]
    if missing:
        raise FileNotFoundError(
            "Dataset(s) not found:\n  "
            + "\n  ".join(missing)
            + f"\nSet V2X_DATA_DIR (currently: {DATA_DIR}). See README.md."
        )
