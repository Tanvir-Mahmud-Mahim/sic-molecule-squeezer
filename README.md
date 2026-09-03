# sic-molecule-squeezer

**End-to-end open-source design of a 4H-SiC photonic-molecule source of
multimode squeezed light from soliton-crystal microcombs.**

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache--2.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21995634.svg)](https://doi.org/10.5281/zenodo.21995634)
[![Paper](https://img.shields.io/badge/Paper-10.1364%2FOE.612248-blue)](https://doi.org/10.1364/OE.612248)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org)

This repository contains the complete simulation pipeline for the article

> T. M. Mahim, M. M. Rahman, and A. S. M. Mohsin, "Overcoming the 3 dB squeezing extraction limit in silicon carbide microcombs with a photonic molecule," Optics Express 34(18), 34822-34834 (2026). https://doi.org/10.1364/OE.612248 (open access)

The pipeline runs from public material data to detected quantum noise:

```
Sellmeier data (CC0)  ->  FEM waveguide modes  ->  ring dispersion D1..D4
                      ->  Lugiato-Lefever soliton crystal (co-moving frame)
                      ->  linearized multimode Heisenberg-Langevin model
                      ->  squeezing spectra, supermodes, entanglement
```

The raw simulation **data and the validation testbench** are archived
separately on Zenodo: **doi:10.5281/zenodo.21995634**.

---

## Repository layout

| Path | Contents |
|---|---|
| `src/materials.py` | Sellmeier models (4H-SiC Wang 2013; SiO2 Malitson 1965; refractiveindex.info, CC0) |
| `src/fem_modes.py` | femwell/scikit-fem full-vector FEM mode solver wrapper |
| `src/sweep_fem.py` | geometry sweep (widths x thicknesses) |
| `src/fem_final.py` | high-accuracy run for the selected geometry, effective area, mode field, mesh convergence |
| `src/dispersion.py` | ring resonance grid and D1, D2, D3, D4 extraction |
| `src/lle.py` | modal Lugiato-Lefever solver, soliton-crystal ansatz, drift measurement |
| `src/exp_lle.py` | crystal continuation across detuning; D3 scaling family |
| `src/exp_d3boundary.py` | D3 stability-boundary scan (step 0.25) and breathing-mechanism diagnosis |
| `src/quantum.py` | linearized multimode quantum model (input-output, covariance, log-negativity) |
| `src/exp_quantum.py` | production quantum experiments (bare device, molecule sweep, entanglement, D3 study) |
| `src/exp_addendum.py` | ideal-detection and molecule squeezing across the detuning family |
| `src/aux_ring.py` | two-ring resonance alignment: mismatch, heater offset, kappa_P(mu), radius-error tolerance |
| `src/exp_tworing.py` | full (non-adiabatic) two-ring quantum model of the photonic molecule |
| `src/validate_quantum.py` | analytic testbench (vacuum, Bogoliubov, 3 dB bound, uncertainty) |
| `src/convergence_checks.py` | numerical convergence testbench (LLE grid/step, mode truncation) |
| `src/make_numbers.py` | regenerates every number quoted in the paper from the data |
| `src/fig_*.py`, `src/figstyle.py` | publication figure generation |
| `figures/` | rendered figures (PNG previews) |

## Installation

Requires Python 3.10+ and a C/C++ toolchain for `gmsh`.

```bash
git clone https://github.com/Tanvir-Mahmud-Mahim/sic-molecule-squeezer.git
cd sic-molecule-squeezer
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

On Debian/Ubuntu, `gmsh` additionally needs system libraries:

```bash
sudo apt-get install libglu1-mesa libxrender1 libxcursor1 libxft2 libxinerama1
```

## Quick start

Reproduce the full study (about 1-2 h on a modern laptop):

```bash
cd src
python sweep_fem.py         # FEM geometry sweep          (~5 min)
python fem_final.py         # selected geometry, order-2  (~5 min)
python exp_lle.py           # soliton-crystal families    (~40 min)
python exp_d3boundary.py    # D3 stability boundary       (~15 min)
python exp_quantum.py       # quantum experiments         (~15 min)
python exp_addendum.py      # detuning-family squeezing   (~10 min)
python aux_ring.py          # two-ring alignment           (~2 min)
python exp_tworing.py       # full two-ring quantum model  (~10 min)
python validate_quantum.py  # analytic testbench          (~1 min)
python convergence_checks.py# convergence testbench       (~15 min)
python make_numbers.py      # regenerate paper numbers
for f in fig_abstract fig_device fig_comb fig_quantum fig_molecule fig_tworing fig_d3; do
    python $f.py            # figures -> ../figures/
done
```

To skip the compute and analyze the published outputs directly, download
the Zenodo archive (doi:10.5281/zenodo.21995634) and unpack its `data/`
folder into `src/`.

## Testbench

Two scripts verify the physics and the numerics:

- `validate_quantum.py` checks the quantum solver against closed-form
  results: exact vacuum output (machine precision), Bogoliubov drift
  eigenvalues of the pumped resonator, the 3 dB critical-coupling bound,
  and the Heisenberg uncertainty product.
- `convergence_checks.py` verifies LLE grid/time-step convergence and the
  fluctuation-mode truncation dependence of the quoted squeezing.

Expected outputs are archived with the Zenodo record under
`testbench/expected_output/`.

## Key results

| Quantity | Value |
|---|---|
| Selected geometry | 1.85 um x 500 nm 4H-SiC, oxide clad, R = 50.06 um |
| Dispersion | D2/2pi = 6.89 MHz, D3/2pi = -63.4 kHz (FSR 350 GHz) |
| Detected squeezing | 7.9 dB at 8.3 mW pump (full two-ring model, kappa_P = 10 kappa, kappa_aux/2pi = 8 GHz; adiabatic limit 8.5 dB) |
| Squeezing bandwidth | 1.81 GHz (full width at half-peak noise reduction) |
| Entanglement | odd-mode lattice, E_N up to 0.13 at the drop port |
| D3 tolerance | crystal survives 3.25x design D3, lost at 3.5x (~222 kHz) via a breathing instability |

## Citing

If you use this code, please cite the article and the software/data
archive (see `CITATION.cff`):

```bibtex
@article{Mahim2026SiCSqueezer,
  author  = {Mahim, Tanvir M. and Rahman, M. Mosaddequr and Mohsin, Abu S. M.},
  title   = {Overcoming the 3 dB squeezing extraction limit in silicon carbide microcombs with a photonic molecule},
  journal = {Optics Express},
  volume  = {34},
  number  = {18},
  pages   = {34822--34834},
  year    = {2026},
  doi     = {10.1364/OE.612248},
  note    = {Code: https://github.com/Tanvir-Mahmud-Mahim/sic-molecule-squeezer;
             Data: doi:10.5281/zenodo.21995634}
}
```

## License

Code: [Apache-2.0](LICENSE) (includes an express patent grant from the
authors; see the NOTICE file). The material-dispersion coefficients derive
from the refractiveindex.info database (CC0). Archived data and testbench
outputs on Zenodo are CC-BY-4.0.
