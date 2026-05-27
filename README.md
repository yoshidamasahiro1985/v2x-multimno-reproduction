# Reproduction Package — When Does Multi-Operator Cellular Redundancy Help V2X?

Masahiro Yoshida · Graduate School of Global Informatics, Chuo University

---

## What is this?

Companion code for the manuscript

> **When Does Multi-Operator Cellular Redundancy Help V2X? A Cross-Layer Measurement Study**

This repository contains the analysis scripts and pre-generated figures needed to
reproduce every result, table, and in-text number reported in the manuscript.

Raw datasets are **not** redistributed here — they are publicly available under their
respective licences (see below and `code/README.md`).

## Repository layout

```
code/          Analysis scripts and configuration (Python 3)
               └── README.md   ← start here for setup + run instructions
figures/       Pre-generated PDF figures (byte-identical to manuscript)
```

## Quick start

See **[`code/README.md`](code/README.md)** for the full setup and run instructions.

The short version:

```bash
cd code
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
V2X_DATA_DIR=/path/to/v2x-data/external bash make_all.sh
```

Running `make_all.sh` re-generates all figures to `../figures/` and prints
all numbers to `numbers/*.txt`.

## Datasets

| Dataset | DOI / source |
|---|---|
| Vienna 4G/5G drive test | [10.48550/arXiv.2603.02638](https://doi.org/10.48550/arXiv.2603.02638) (CC BY 4.0) |
| Berlin V2X convoy | [10.1109/VTC2023-Spring57618.2023.10200750](https://doi.org/10.1109/VTC2023-Spring57618.2023.10200750) (IEEE DataPort) |
| DoNext 5G-NSA drive test | [10.1109/TMLCN.2025.3564239](https://doi.org/10.1109/TMLCN.2025.3564239) (open access) |

## Citation

If you use this code, please cite this reproduction package via its Zenodo DOI:

> **DOI**: [10.5281/zenodo.20404338](https://doi.org/10.5281/zenodo.20404338)

Once the accompanying manuscript has a public DOI, please cite that too. See
`code/CITATION.cff` for the machine-readable citation.

BibTeX entry for the manuscript (update with DOI on acceptance):

```bibtex
@article{yoshida2026v2x,
  title   = {When Does Multi-Operator Cellular Redundancy Help V2X?
             A Cross-Layer Measurement Study},
  author  = {Yoshida, Masahiro},
  journal = {IEEE Access},
  year    = {2026},
  note    = {Submitted}
}
```

## Licence

Code: MIT — see `code/LICENSE`.
