# sic-molecule-squeezer: Simulation Code

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)
[![Paper DOI](https://img.shields.io/badge/Paper-10.1364%2FOE.612248-blue)](https://doi.org/10.1364/OE.612248)
[![Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21995634.svg)](https://doi.org/10.5281/zenodo.21995634)

This repository holds the code for our article **"Overcoming the 3 dB
squeezing extraction limit in silicon carbide microcombs with a photonic
molecule"**, which I wrote with M. Mosaddequr Rahman and Abu S. M. Mohsin
(Department of Electrical and Electronic Engineering, BRAC University).
It was published in *Optics Express* **34**(18), 34822-34834 (2026).

- Article (open access): https://doi.org/10.1364/OE.612248
- Supplemental document: [`supplement.pdf`](supplement.pdf) in this repository
- Repository: https://github.com/Tanvir-Mahmud-Mahim/sic-molecule-squeezer
- Archived data and testbench: https://doi.org/10.5281/zenodo.21995634

With these scripts you can run the whole study. It goes from published
material data all the way to the quantum noise that a detector sees:

```
material data  ->  waveguide modes (FEM)  ->  ring dispersion D1..D4
               ->  soliton crystal (Lugiato-Lefever equation)
               ->  linearized multimode quantum model
               ->  squeezing spectra, supermodes, entanglement
```

One figure, Supplementary Fig. S1, has no script here. I also found two
statements in the published supplement that do not match what the code
and the archived data give. They are the first two rows of Table S1, and
the breathing frequency in Section 6. I list them in
[Section 10](DETAILS.md#10-notes-on-the-calculations) under "Known gaps" and
"Differences from the published supplement".

---

## Contents

1. [The idea in one minute](#1-the-idea-in-one-minute)
2. [What is in this repository](#2-what-is-in-this-repository)
3. [Installation](#3-installation)
4. [Quick start: three ways to use the code](#4-quick-start-three-ways-to-use-the-code)
5. [The scripts, step by step](DETAILS.md#5-the-scripts-step-by-step)
6. [Which script makes which figure](DETAILS.md#6-which-script-makes-which-figure)
7. [The Python modules](DETAILS.md#7-the-python-modules)
8. [Where the numbers come from](DETAILS.md#8-where-the-numbers-come-from)
9. [Built-in checks](DETAILS.md#9-built-in-checks)
10. [Notes on the calculations](DETAILS.md#10-notes-on-the-calculations)
11. [Version history](DETAILS.md#11-version-history)
12. [How to cite](#12-how-to-cite)
13. [License and contact](#13-license-and-contact)

Sections 5 to 11 are in [DETAILS.md](DETAILS.md). So are the
[extra notes for Sections 1 to 4](DETAILS.md#extra-notes-for-sections-1-to-4).

---

## 1. The idea in one minute

**Squeezed light** is light whose random noise is pushed below the noise
of empty space (the "vacuum" or "shot-noise" level). This happens in one
of its two wave properties, called quadratures. Squeezing is quoted in
decibels (dB) below the vacuum level. More dB means quieter light.

A **microcomb** is a tiny ring of optical material, here silicon carbide
(4H-SiC). A single laser turns it into many evenly spaced colors, like the
teeth of a comb. The spacing of the teeth is the **free spectral range
(FSR)**. In our design, two light pulses circulate in the ring, evenly
spaced (a "2-FSR soliton crystal"). Then only every second comb tooth is
bright (the even modes). The odd modes in between carry no laser light,
but the bright teeth feed them with squeezed vacuum.

Here is the problem. Take a ring coupled to one output waveguide at the
usual balanced ("critical") setting. Half of the light is lost inside the
ring, and only half reaches the detector. This caps the squeezing you can
detect at 3 dB, however strong the source.

The fix we study is a **photonic molecule**: a second, smaller ring (half
the radius, so twice the FSR) placed next to the main ring. Its resonances
line up only with the odd modes. So it pulls the squeezed odd modes out
quickly to a separate "drop" waveguide, before they are lost. The bright
even modes stay in the main ring.

Here are the main results as recorded in the repository (the previous
README, the supplement, and the values drawn in `fig_abstract.py`):

| Quantity | Value |
|---|---|
| Detected squeezing | 7.9 dB at 8.3 mW pump (full two-ring model, κP = 10κ, κaux/2π = 8 GHz); 8.5 dB in the simplified (adiabatic) model |
| Squeezing bandwidth | 1.81 GHz (full width at half of the peak noise reduction) |
| Entanglement | odd-mode lattice, log-negativity up to 0.13 at the drop port |

(κ is the loss rate of a main-ring mode. κP is the extra extraction rate
that the small ring gives the odd modes. κaux is the loss rate of the small
ring into its drop waveguide. I explain these symbols again in
[Section 10](DETAILS.md#10-notes-on-the-calculations).)

The full results table and a list of what the code computes are in
[DETAILS.md](DETAILS.md#more-on-section-1-what-the-code-computes-and-the-main-results).

---

## 2. What is in this repository

The main parts are:

- `src/`: all Python scripts and modules.
- `figures/`: PNG previews of the article figures.
- `supplement.pdf`: supplemental document of the article (6 pages).
- `requirements.txt`: Python packages to install.
- [CHANGELOG.md](CHANGELOG.md): what changed, newest first.
- `CITATION.cff`: citation details (drives the "Cite this repository" button).
- `LICENSE` and `NOTICE`: the Apache-2.0 license; copyright holders and
  credit for the CC0 material data.

You run all scripts from inside `src/`. They read and write their results
(`.json` and `.npz` files) in `src/` itself. The figure scripts write to
`figures/` (both PDF and PNG; only the PNG previews are stored on GitHub).

The full annotated file tree and more notes on the result files are in
[DETAILS.md](DETAILS.md#more-on-section-2-the-full-file-tree).

---

## 3. Installation

You need **Python 3.10 or newer**. The repository states this, and the
current releases of `scikit-fem` and `shapely` also require Python 3.10 or
newer. I checked this guide with Python 3.11.

```
git clone https://github.com/Tanvir-Mahmud-Mahim/sic-molecule-squeezer.git
cd sic-molecule-squeezer
python -m venv venv
source venv/bin/activate          # on Windows: venv\Scripts\activate
pip install -r requirements.txt
```

`requirements.txt` asks for `numpy>=1.24`, `scipy>=1.10`,
`matplotlib>=3.7`, `shapely>=2.0`, `scikit-fem>=8.0`, `gmsh>=4.11` and
`femwell>=0.1`. The versions are not pinned. The exact versions used for
the article are not recorded either (the supplement says "femwell (v0.x)").

My notes from checking this guide are in
[DETAILS.md](DETAILS.md#more-on-section-3-what-i-noticed-while-checking-the-installation).
They cover disk space for `femwell`, which steps need which packages, and
the system libraries that `gmsh` needs on Debian/Ubuntu.

---

## 4. Quick start: three ways to use the code

Run all commands from the `src/` folder:

```
cd src
```

### Way A: check that the code works (a few seconds)

```
python validate_quantum.py
python materials.py
python fem_modes.py
```

`validate_quantum.py` needs no input files. It prints the results of the
exact-result checks of the quantum model (explained in
[Section 9](DETAILS.md#9-built-in-checks)). It does not print PASS or FAIL, so
compare its output with the expected values given there. `materials.py`
prints the refractive indices at 1.30, 1.55 and 1.80 µm. `fem_modes.py`
solves one waveguide mode (1.85 µm by 0.5 µm at 1.55 µm). It should print
an effective index of about 2.30724 with 4612 mesh elements. This matches
the 40 nm, first-order row of Table S1 in the supplement.

### Way B: work from the archived results

The previous README gives this route. Download the Zenodo archive
(https://doi.org/10.5281/zenodo.21995634) and unpack its `data/` folder
into `src/`, so that the `.json` and `.npz` files sit next to the scripts.
Then you can run any analysis or figure script from Section 5 and
Section 6, for example:

```
python fig_comb.py
```

What I found in the archive is in
[DETAILS.md](DETAILS.md#way-b-the-archive-contents).

### Way C: recompute everything from scratch

Run the steps of [Section 5](DETAILS.md#5-the-scripts-step-by-step) in order. The
previous README estimated about 1 to 2 hours in total on a modern laptop.
My measured run times are in [DETAILS.md](DETAILS.md#way-c-measured-run-times).

On Linux or macOS you can write the same order as:

```
for s in sweep_fem fem_final exp_lle exp_d3boundary exp_quantum \
         exp_addendum aux_ring exp_tworing exp_tworing_family exp_d3mechanism \
         validate_quantum convergence_checks make_numbers; do
    python $s.py
done
for f in fig_abstract fig_device fig_comb fig_quantum fig_molecule fig_tworing fig_d3; do
    python $f.py
done
```

---

## More details (Sections 5 to 11)

The full notes are in [DETAILS.md](DETAILS.md). Here is what each section holds:

- [5. The scripts, step by step](DETAILS.md#5-the-scripts-step-by-step):
  a table of every script, with what it needs, what it writes and its run time.
- [6. Which script makes which figure](DETAILS.md#6-which-script-makes-which-figure):
  which script draws each figure file, and what each figure shows.
- [7. The Python modules](DETAILS.md#7-the-python-modules):
  what each module contains, and which scripts import other scripts.
- [8. Where the numbers come from](DETAILS.md#8-where-the-numbers-come-from):
  sources of the material data, values set in the code, and numerical settings.
- [9. Built-in checks](DETAILS.md#9-built-in-checks):
  the checks the scripts print, and my reproduction checks of 30 September 2026.
- [10. Notes on the calculations](DETAILS.md#10-notes-on-the-calculations):
  symbols and units, the simplified and full molecule models, approximations,
  hard-coded values, known gaps, and differences from the published supplement.
- [11. Version history](DETAILS.md#11-version-history):
  versions, dates and data DOIs.

---

## 12. How to cite

Please cite the article and the data archive. GitHub also shows a
**"Cite this repository"** button in the right-hand column, which reads
`CITATION.cff`.

> T. M. Mahim, M. M. Rahman, and A. S. M. Mohsin, "Overcoming the 3 dB
> squeezing extraction limit in silicon carbide microcombs with a photonic
> molecule," Optics Express 34(18), 34822-34834 (2026).
> https://doi.org/10.1364/OE.612248
>
> Data and testbench: Zenodo, https://doi.org/10.5281/zenodo.21995634

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

---

## 13. License and contact

Code: Apache License 2.0 (see `LICENSE`; Section 3 of the license contains
the contributors' patent grant). `NOTICE` names the copyright holders and
credits the refractiveindex.info database (CC0) for the material
coefficients. Archived data and testbench outputs on Zenodo: CC BY 4.0.

If you have questions or find a bug, please open an issue on this
repository. You can also contact me, Tanvir M. Mahim, at BRAC University
(tanvir.mahim@bracu.ac.bd).
