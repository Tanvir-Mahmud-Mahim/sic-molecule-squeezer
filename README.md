# Overcoming the 3 dB squeezing extraction limit in silicon carbide microcombs with a photonic molecule

End-to-end open-source simulation study: FEM waveguide dispersion -> Lugiato-Lefever
soliton-crystal mean field -> linearized multimode quantum model -> detected squeezing
and entanglement at a photonic-molecule drop port.

Code: https://github.com/Tanvir-Mahmud-Mahim/sic-molecule-squeezer (Apache-2.0)
Data + testbench: https://doi.org/10.5281/zenodo.21471674 (CC-BY-4.0)

## Contents

- `paper/` - Optics Express manuscript (Optica universal template)
  - `manuscript.tex` - single self-contained top-level file (all numbers and
    links inlined as macros; manual DOI-hyperlinked bibliography, 33 references)
  - `manuscript.pdf` (10 pages, abstract 145 words)
  - `fig0_abstract.pdf` (graphical abstract), `fig1_device.pdf` ... `fig5_d3.pdf`
  - Table 1: quantitative figure-of-merit comparison with the state of the art
- `supplement/` - supplemental document (Optica supplemental template)
  - `supplement.tex` - single self-contained top-level file
  - `supplement.pdf` (4 pages), `figS1.pdf`
- `sim/` - complete simulation pipeline (Python 3, open source only)
  - `materials.py` - Sellmeier models (refractiveindex.info, CC0)
  - `fem_modes.py`, `sweep_fem.py`, `fem_final.py` - femwell/scikit-fem FEM mode solving
  - `dispersion.py` - ring resonance grid, D1..D4 extraction
  - `lle.py`, `exp_lle.py` - modal LLE, co-moving-frame soliton crystals
  - `quantum.py` - linearized Heisenberg-Langevin multimode model
  - `exp_quantum.py`, `exp_addendum.py` - production quantum experiments
  - `validate_quantum.py` - analytic validation tests
  - `convergence_checks.py` - numerical convergence studies
  - `make_numbers.py` - regenerates all quoted numbers from data
  - `fig_*.py`, `figstyle.py` - figure generation
  - `*.json`, `*.npz`, `*.log` - all raw simulation data and run logs
- `figures/` - all figures as PDF and PNG
- `upload/` - published packages and the guide used to publish them
  - `github-repo/sic-molecule-squeezer/` - code package as uploaded to GitHub
    (README, Apache-2.0 LICENSE, NOTICE, CITATION.cff, .zenodo.json)
  - `zenodo-dataset/sic-molecule-squeezer-data-v1.0/` - data + testbench package
    as published on Zenodo (CC-BY-4.0)
  - `UPLOAD_GUIDE.md` - how both were published and wired into the paper
- `opex_submission_main.zip` / `opex_submission_supplement.zip` - flat,
  single-.tex submission packages for Optics Express (verified to compile standalone)
- `SiC_squeezing_project_final.zip` - this whole project, archived
- `UPLOAD_RECORD.md` - record of how the GitHub and Zenodo uploads were performed

## Key results

- Selected geometry: 1.85 um x 500 nm oxide-clad 4H-SiC, R = 50.06 um, FSR 350 GHz
- D2/2pi = 6.89 MHz, D3/2pi = -63.4 kHz (FEM, converged to <1e-7 in neff)
- Stationary 2-FSR soliton crystal for zeta0 = 5.75-11.0 at 8.3 mW pump
- D3-induced repetition-rate shift: -1.3 to -2.2 MHz (linear in D3; in-situ diagnostic)
- Bare critically coupled device: squeezing pinned at the 3 dB coupling limit
- Photonic molecule at fixed pump: optimum Purcell rate ~10 kappa,
  8.5 dB detected squeezing over 1.74 GHz; 9.6 dB near annihilation
- Drop-port entangled odd-mode lattice, log-negativity up to 0.12
- Squeezing robust to D3 (spread < 0.01 dB) until crystal loss at ~3.5x design D3

## Reproducing

```
pip install femwell gmsh shapely scikit-fem matplotlib scipy numpy
cd sim
python sweep_fem.py && python fem_final.py
python exp_lle.py && python exp_quantum.py && python exp_addendum.py
python validate_quantum.py && python convergence_checks.py
python make_numbers.py && python fig_abstract.py && python fig_device.py \
  && python fig_comb.py && python fig_quantum.py && python fig_molecule.py \
  && python fig_d3.py
cd ../paper && pdflatex manuscript && pdflatex manuscript
```

Material data: refractiveindex.info database (CC0): 4H-SiC (Wang 2013,
doi:10.1002/lpor.201300068), SiO2 (Malitson 1965, doi:10.1364/JOSA.55.001205).

## Reference verification

All 33 manuscript references (the 7 supplement references are a subset) were
checked against Crossref and publisher records: every DOI resolves to the cited
paper, and every bibliographic detail (authors, journal, volume, pages, year)
matches. In-text factual claims about cited work were verified against the
sources. No errors found.
