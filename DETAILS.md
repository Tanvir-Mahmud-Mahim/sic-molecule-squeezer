# sic-molecule-squeezer: details

This file has the full notes. For the overview, installation and quick start,
see [README.md](README.md).

## Contents

- [Extra notes for Sections 1 to 4](#extra-notes-for-sections-1-to-4)
  - [More on Section 1: what the code computes and the main results](#more-on-section-1-what-the-code-computes-and-the-main-results)
  - [More on Section 2: the full file tree](#more-on-section-2-the-full-file-tree)
  - [More on Section 3: what I noticed while checking the installation](#more-on-section-3-what-i-noticed-while-checking-the-installation)
  - [More on Section 4: the archive and the run times](#more-on-section-4-the-archive-and-the-run-times)
- [5. The scripts, step by step](#5-the-scripts-step-by-step)
- [6. Which script makes which figure](#6-which-script-makes-which-figure)
- [7. The Python modules](#7-the-python-modules)
- [8. Where the numbers come from](#8-where-the-numbers-come-from)
- [9. Built-in checks](#9-built-in-checks)
- [10. Notes on the calculations](#10-notes-on-the-calculations)
- [11. Version history](#11-version-history)

---

## Extra notes for Sections 1 to 4

### More on Section 1: what the code computes and the main results

The code computes:

- the optical modes and dispersion of the silicon carbide waveguide and ring,
- the soliton-crystal state of the main ring,
- the squeezing, the squeezed "supermodes" (the best combinations of many
  comb lines) and the entanglement between comb lines, for the plain ring
  and for the photonic molecule,
- a full model of the two coupled rings, including misalignment of their
  resonances and fabrication errors in the small ring's radius,
- checks against exact results and numerical convergence checks.

Here are the main results as recorded in the repository (the previous
README, the supplement, and the values drawn in `fig_abstract.py`):

| Quantity | Value |
|---|---|
| Selected waveguide | 4H-SiC, 1.85 µm wide, 500 nm thick, oxide cladding; ring radius R = 50.06 µm |
| Dispersion | D2/2π = 6.89 MHz, D3/2π = -63.4 kHz (FSR 350 GHz) |
| Detected squeezing | 7.9 dB at 8.3 mW pump (full two-ring model, κP = 10κ, κaux/2π = 8 GHz); 8.5 dB in the simplified (adiabatic) model |
| Squeezing bandwidth | 1.81 GHz (full width at half of the peak noise reduction) |
| Entanglement | odd-mode lattice, log-negativity up to 0.13 at the drop port |
| Tolerance to third-order dispersion D3 | crystal survives at 3.25 times the design D3 and is lost at 3.5 times (about 222 kHz) |

(κ is the loss rate of a main-ring mode. κP is the extra extraction rate
that the small ring gives the odd modes. κaux is the loss rate of the small
ring into its drop waveguide. I explain these symbols again in
[Section 10](#10-notes-on-the-calculations).)

### More on Section 2: the full file tree

```
sic-molecule-squeezer/
|-- README.md              this guide
|-- CHANGELOG.md           what changed, newest first
|-- CITATION.cff           citation details (drives the "Cite this repository" button)
|-- .zenodo.json           metadata for the Zenodo software record (article title outdated; see Section 10)
|-- LICENSE                Apache-2.0 license
|-- NOTICE                 copyright holders; credit for the CC0 material data
|-- requirements.txt       Python packages to install
|-- supplement.pdf         supplemental document of the article (6 pages)
|-- figures/               PNG previews of the article figures
|   |-- fig0_abstract.png  graphical abstract
|   |-- fig1_device.png    device concept and waveguide design
|   |-- fig2_comb.png      the soliton crystal
|   |-- fig3_quantum.png   squeezing of the plain ring
|   |-- fig4_molecule.png  squeezing with the photonic molecule
|   |-- fig5_d3.png        third-order dispersion and entanglement map
|   |-- fig6_tworing.png   full two-ring model
|   `-- figS1.png          supplementary figure S1 (no script in this repository)
`-- src/
    |-- materials.py           refractive index of 4H-SiC and silica
    |-- fem_modes.py           waveguide mode solver (femwell / scikit-fem)
    |-- sweep_fem.py           step 1: scan waveguide widths and thicknesses
    |-- fem_final.py           step 2: accurate run for the chosen waveguide
    |-- dispersion.py          ring resonances and dispersion D1..D4
    |-- lle.py                 Lugiato-Lefever solver and soliton starting shape
    |-- exp_lle.py             step 3: soliton crystals across detuning; D3 family
    |-- exp_d3boundary.py      step 4: where the crystal is lost as D3 grows
    |-- quantum.py             linearized multimode quantum model
    |-- exp_quantum.py         step 5: squeezing, supermodes, entanglement, D3 study
    |-- exp_addendum.py        step 6: ideal detection and molecule across detuning
    |-- aux_ring.py            step 7: alignment of the two rings, radius errors
    |-- exp_tworing.py         step 8: full two-ring quantum model
    |-- exp_tworing_family.py  step 9: full two-ring model across detuning
    |-- exp_d3mechanism.py     step 10: breathing state just past the D3 boundary
    |-- validate_quantum.py    checks of quantum.py against exact results
    |-- convergence_checks.py  numerical convergence checks
    |-- make_numbers.py        writes the article's numbers as LaTeX macros (numbers/numbers.tex)
    |-- figstyle.py            shared figure style
    `-- fig_*.py               one script per figure (see Section 6)
```

`make_numbers.py` writes `numbers/numbers.tex` in the repository folder
(it creates `numbers/` when needed). The result files are not stored on
GitHub (`.gitignore` excludes `src/*.npz`, `src/*.json` and `numbers/`).
You can get the archived copy from Zenodo. The FEM scripts also leave a
temporary mesh file `mesh.msh` in `src/` (also ignored by git).

### More on Section 3: what I noticed while checking the installation

Here is what I noticed while checking this guide (30 September 2026):

- `femwell` 0.1.12 depends on `meshwell`, which pulls in large packages
  (for example `cadquery` and `vtk`). So allow plenty of disk space. The
  code here only imports `femwell.mesh` and `femwell.maxwell.waveguide`.
- Only the first two steps (the FEM mode solving) need `femwell`,
  `scikit-fem`, `gmsh` and `shapely`. Everything after that needs only
  `numpy`, `scipy` and `matplotlib`.
- `gmsh` needs some system graphics libraries, even when nothing is drawn
  on screen. On Debian/Ubuntu you can install them with:

  ```
  sudo apt-get install libglu1-mesa libxrender1 libxcursor1 libxft2 libxinerama1 libopengl0
  ```

  (The earlier README listed all of these except `libopengl0`. On the
  Ubuntu 24.04 machine I used for checking, `gmsh` also needed
  `libOpenGL.so.0`, which `libopengl0` provides.)

### More on Section 4: the archive and the run times

#### Way B: the archive contents

I checked the unpacked archive on 30 September 2026. Its `data/` folder
has the 15 result files the scripts write (including
`q_tworing_family.npz` and `d3_mech.npz`). Its `testbench/` folder has a
frozen copy of the scripts (`testbench/src/`) and reference logs
(`testbench/expected_output/`). The archive does not contain the two
scripts added in the fix of 30 September 2026 (`exp_tworing_family.py`,
`exp_d3mechanism.py`).

#### Way C: measured run times

When I checked this guide, steps 1 to 8 plus `convergence_checks.py` took
about 70 minutes together on a shared 2-core machine. Each figure script
took a few seconds. Steps 9 and 10 (added later) took about 9 minutes more
(see the table in Section 5).

---

## 5. The scripts, step by step

"Repo estimate" is the run time given in the previous README (for "a
modern laptop"). "Measured" is the time I measured while checking this
guide (30 September 2026). I used a shared 2-core machine that was also
running other jobs, so a machine to yourself may well be faster.

| Step | Command (in `src/`) | What it does | Needs | Writes | Repo estimate | Measured |
|---|---|---|---|---|---|---|
| 1 | `python sweep_fem.py` | Effective index of the fundamental TE mode (electric field lying in the chip plane) at 13 wavelengths from 1.30 to 1.85 µm, for widths 1.40, 1.60, 1.85, 2.10 and 2.40 µm and thicknesses 0.50 and 0.60 µm (first-order elements) | - | `fem_sweep.json` | ~5 min | 14 min |
| 2 | `python fem_final.py` | Chosen waveguide (1.85 µm x 0.50 µm): mesh-convergence table at 1.55 µm, second-order run at 17 wavelengths, effective mode area, mode-field map | - | `fem_final.json`, `mode_field.npz` | ~5 min | 11.5 min |
| 3 | `python exp_lle.py` | Sets the ring radius for a 350 GHz FSR and extracts D1..D4. Finds the two-soliton crystal and follows it in pump detuning from 6.0 upward (toward 12.0) and from 5.75 downward (toward 2.0) in steps of 0.25, stopping in each direction once the crystal is lost. Then runs a D3-scaling family (0, 1, 2, 2.5, 3 times the design D3) at detuning 6.5 | `fem_final.json` | `lle_family.npz`, `lle_d3family.npz` | ~40 min | 16 min |
| 4 | `python exp_d3boundary.py` | Scans the D3 scale from 3.0 to 4.0 in steps of 0.25 to find where the crystal is lost, and tracks how the odd modes grow just past that point | `fem_final.json` | `d3_boundary.npz` | ~15 min | 5.7 min |
| 5 | `python exp_quantum.py` | Plain ring: squeezing across all stationary crystals. At detuning 6.5: spectra, supermodes, photonic-molecule scan (κP = 0.5 to 50 κ), entanglement matrix at κP = 20κ, and the D3 family at κP = 20κ | `lle_family.npz`, `lle_d3family.npz` | `q_sweep.npz`, `q_rep.npz`, `q_entanglement.npz`, `q_d3.npz` | ~15 min | 1.4 min |
| 6 | `python exp_addendum.py` | For every stationary crystal: squeezing with ideal detection (all loss monitored) and with the molecule at κP = 10κ and 20κ | `lle_family.npz` | `q_addendum.npz` | ~10 min | 2.6 min |
| 7 | `python aux_ring.py` | Resonances of the small ring (radius for a 700 GHz FSR), one heater offset that best aligns them with the odd main-ring modes, the leftover mismatch, and the effect of 5 nm and 20 nm radius errors | `fem_final.json`, `lle_family.npz` | `aux_ring.npz` | ~2 min | 28 s |
| 8 | `python exp_tworing.py` | Full two-ring quantum model: comparison with the simplified model for κP = 5 to 30 κ and κaux/2π = 2, 4, 8 GHz; convergence as κaux grows; design point (κP = 10κ, 8 GHz) with bandwidth, supermodes, entanglement and radius errors | `lle_family.npz`, `aux_ring.npz` | `q_tworing.npz` | ~10 min | 15 min |
| 9 | `python exp_tworing_family.py` | Full two-ring model at the design point (κP = 10κ, κaux/2π = 8 GHz, perfect alignment) for every stationary crystal (ζ0 = 5.75 to 11.0): peak squeezing at the drop port | `lle_family.npz`, `aux_ring.npz` | `q_tworing_family.npz` | - | 7.3 min |
| 10 | `python exp_d3mechanism.py` | Converges the crystal at 3.25 times D3 (ζ0 = 6.5), switches to 3.5 times D3 and follows the resulting breathing state for 150 time units (2/κ). Records the peak field and comb lines μ = 2 and 88, their mean and oscillation, the dominant breathing frequency, and where Dint changes sign | `fem_final.json` | `d3_mech.npz` | - | 1.7 min |
| - | `python validate_quantum.py` | Checks of the quantum model against exact results (Section 9) | - | printed only | ~1 min | 1.4 s |
| - | `python convergence_checks.py` | Soliton solver: grid size and time step; quantum model: number of modes kept (M = 20, 30, 40) | `fem_final.json`, `lle_family.npz` | printed only | ~15 min | 3.7 min |
| - | `python make_numbers.py` | Collects every quoted number into LaTeX macros (also printed as JSON) | `fem_final.json` and the outputs of steps 3 to 10 | `../numbers/numbers.tex` (folder created if needed) | - | 1.4 s |
| - | `python fig_*.py` | Draws the figures (Section 6) | see Section 6 | `../figures/*.pdf`, `../figures/*.png` | - | see Section 6 |

No script uses random numbers.

---

## 6. Which script makes which figure

The table lists the file each script writes, and what the script's own
description says it shows. Be careful: the figure numbers in the file
names are not always the figure numbers in the published article. For
example, the supplement calls the third-order-dispersion figure
(`fig5_d3`) "Fig. 7" of the main text. Please check the article for the
final numbering.

| Figure file | Content | Data from | Drawn by | Measured |
|---|---|---|---|---|
| `fig0_abstract` | Graphical abstract: device sketch, soliton crystal, bare versus molecule squeezing spectra, headline numbers (the three numbers are typed into the script) | steps 3, 5, 8 | `fig_abstract.py` | 2 s |
| `fig1_device` | (a) 3D sketch of the two rings, (b) which modes the small ring extracts, (c) mode field in the waveguide cross-section, (d) integrated dispersion, (e) D2 versus width for both thicknesses | steps 1, 2 | `fig_device.py` | 8 s |
| `fig2_comb` | Soliton crystal: (a) intensity around the ring, (b) comb spectrum, (c) where the crystal exists versus detuning, (d) shift of the comb spacing caused by D3 | step 3 | `fig_comb.py` | 3 s |
| `fig3_quantum` | Plain ring: (a) squeezing and anti-squeezing spectra, (b) the two strongest supermodes, (c) squeezing versus detuning with the 3 dB bound and ideal detection | steps 5, 6 | `fig_quantum.py` | 2 s |
| `fig4_molecule` | Molecule: (a) spectra for several κP, (b) peak squeezing versus κP, (c) bandwidth versus κP, (d) squeezing versus detuning | steps 5, 6 | `fig_molecule.py` | 3 s |
| `fig5_d3` | (a) comb-spacing shift versus D3 with the loss point, (b) spectra for all D3 scales, (c) supermode asymmetry, (d) entanglement map between odd modes | steps 3, 5 | `fig_d3.py` | 3 s |
| `fig6_tworing` | Full two-ring model: (a) mismatch after heater alignment, (b) resulting extraction rate per mode, (c) spectra, simplified versus full model, (d) peak squeezing versus κaux and radius-error points | steps 7, 8 | `fig_tworing.py` | 12 s |
| `figS1` | Supplementary Fig. S1: supermodes at the drop port (κP = 20κ) and spectra of bare ring and molecule | - | no script in this repository | - |

---

## 7. The Python modules

| File | What it contains |
|---|---|
| `materials.py` | Sellmeier formulas (refractive index versus wavelength) for 4H-SiC (ordinary and extraordinary ray) and fused silica |
| `fem_modes.py` | Builds the cross-section mesh (SiC core in a 6.4 µm x 4.9 µm silica window) and solves for the fundamental TE mode with femwell; `make_mesh`, `solve_neff` |
| `dispersion.py` | `RingDispersion`: fits the effective index with a degree-6 polynomial, finds ring resonances from β(ω)·2πR = 2πm, and fits D1, D2, D3, D4 over mode numbers \|μ\| ≤ 40; also the radius for a given FSR |
| `lle.py` | `LLE`: split-step solver of the Lugiato-Lefever equation in the mode basis; `measure_drift` (speed of a drifting pattern); `cw_background`; `soliton_crystal_ansatz` (starting shape) |
| `quantum.py` | `CombQuantum`: linearized quantum model around the comb; output noise spectra, squeezing spectra, supermodes, intracavity covariance; `log_negativity` (entanglement measure for two modes) |
| `figstyle.py` | Shared figure style (Okabe-Ito color-blind-safe colors, DejaVu Sans font, 5.9 inch full width) and the `save` helper |

Three of the `exp_*.py` scripts also serve as modules for other scripts:

- `exp_lle.py` (`converge_state`, `is_crystal`, `zeta_of`) is used by
  `exp_d3boundary.py`, `exp_d3mechanism.py` and `convergence_checks.py`;
- `exp_quantum.py` (`build`, `max_squeezing`, `centered_psi`, `zeta_vec`)
  is used by `exp_addendum.py`, `exp_tworing.py` and `convergence_checks.py`;
- `exp_tworing.py` (`build_tworing`, `peak`) is used by
  `exp_tworing_family.py`.

All three load files as soon as you import them. Importing `exp_lle.py`
loads `fem_final.json`. Importing `exp_quantum.py` loads `lle_family.npz`.
Importing `exp_tworing.py` loads `lle_family.npz` and `aux_ring.npz`. So
those files must already exist.

---

## 8. Where the numbers come from

**Material data** (`materials.py`, from the refractiveindex.info database,
CC0 public domain):

- 4H-SiC: S. Wang et al., Laser Photonics Rev. 7, 831 (2013),
  doi:10.1002/lpor.201300068 (ordinary ray used for the TE mode; range
  0.4047 to 5 µm).
- Silica: I. H. Malitson, J. Opt. Soc. Am. 55, 1205 (1965),
  doi:10.1364/JOSA.55.001205 (range 0.21 to 6.7 µm).

At 1.55 µm these give n(SiC, ordinary) = 2.5644 and n(SiO2) = 1.4440
(printed by `python materials.py`).

**Values set in `exp_lle.py`**:

| Value | Used | Source as recorded |
|---|---|---|
| Kerr coefficient n2 | 6.9e-19 m²/W | code comment "4H-SiC (Guidry et al.)"; the supplement cites Wang et al. and Guidry et al. |
| Linear index in the Kerr coupling n0 | 2.5644 | the 4H-SiC ordinary index at 1.55 µm |
| Pump wavelength | 1.55 µm | set in the code |
| FSR | 350 GHz | set in the code; gives R = 50.06 µm |
| Loaded quality factor Q | 1.5 million | set in the code; no source recorded. It gives κ/2π ≈ 128.9 MHz |
| Pump strength | f = 3 (f² = 9) | set in the code; the supplement quotes 8.3 mW |
| Coupling | critical (η = 1/2) | set in the code |
| Mode-field area | from `fem_final.py` | computed |

**Design values of the photonic molecule** (`exp_quantum.py`,
`aux_ring.py`, `exp_tworing.py`):

- small-ring radius for a 700 GHz FSR;
- κaux/2π = 2 GHz in `aux_ring.py`, scanned 2, 4, 8, 16 and 64 GHz in
  `exp_tworing.py`, with 8 GHz as the design point;
- extraction rate κP from 0.5κ to 50κ, with 10κ as the design point;
- radius errors of 5 nm and 20 nm.

**Numerical settings**: the soliton grid is N = 192 modes, with time step
0.002. The quantum model keeps fluctuation modes |μ| ≤ 30 (M = 30) and comb lines
|μ| ≤ 60. The representative detuning is ζ0 = 6.5.

---

## 9. Built-in checks

**`validate_quantum.py`** (quantum model against exact results):

| Printed line | What it compares | Result when I checked it |
|---|---|---|
| `T1 vacuum` | With no comb, the output noise matrix must equal the vacuum (identity) | max difference 2.2e-16 |
| `T2` | Largest real part of the numerical eigenvalues for a continuously pumped ring versus the exact (Bogoliubov) formula | -0.5 exact, -0.4999999999999985 numerical |
| `T3 critical` | Critically coupled ring: squeezing must stay below 3.01 dB; the product of smallest and largest noise must be at least 1 (uncertainty principle) | 2.84 dB; product 1.004 |
| `T3 overcoupled 10x` | Output coupling 10 times stronger: squeezing may pass 3 dB | 4.49 dB; product 1.002 |

The supplement (Section 8) gives the pass criteria:

- vacuum to better than 1e-15;
- Bogoliubov eigenvalues matching the closed form to 1e-13;
- no more than 3.01 dB at critical coupling;
- the uncertainty product at least 1 "on all frequency grids".

The run above is consistent with these on what the script actually
tests. T1 is below 1e-15. T2 compares only the largest real part of the
eigenvalues with the formula (difference 1.5e-15), not every eigenvalue.
T3 stays below 3.01 dB. The uncertainty product is at least 1 on the one
frequency grid the script uses (61 points from 0 to 3 κ), not on several
grids. The supplement's list of analytic checks (Section 5E) also
includes a "noise-flow identity" between detected and perfectly detected
spectra. No script in this repository contains it.

**`fem_final.py`** prints a mesh-convergence table for the chosen
waveguide at 1.55 µm (compare it with Table S1 of the supplement). I
checked it with gmsh 4.15.2, scikit-fem 12.0.2 and femwell 0.1.12. The
rows for 40 nm and 30 nm with first-order elements matched Table S1 to
all 7 printed digits. So did all three second-order rows (for example
2.3081259 for the production setting). The first two rows did not: 60 nm
gave 2.3064460 (2958 elements), and 50 nm gave 2.3068171 (3426 elements).
Table S1 instead lists 2.3068171 (3426 elements) for 60 nm and 2.3067888
(3426 elements) for 50 nm. The archived `fem_final.json` agrees with the
code, not with the table. See "Differences from the published supplement"
in [Section 10](#10-notes-on-the-calculations).

**`convergence_checks.py`** prints the soliton result for grid sizes 192
and 256 and time steps 0.004, 0.002 and 0.001 (compare with Table S2 of
the supplement). It also prints the peak squeezing for M = 20, 30 and 40
kept modes (compare with Table S4).

**`exp_tworing.py`** starts with a vacuum check of the two-ring model
(the printed value should be close to 0).

**Reproduction check (30 September 2026).** For this guide, I ran steps 1
to 8 and `convergence_checks.py` from scratch (Python 3.11, numpy 2.4.4,
scipy 1.17.1, matplotlib 3.10.9). Among other things, they reproduced:

- R = 50.06 µm, D2/2π = 6.890 MHz, D3/2π = -63.37 kHz, D4/2π = -148.32 Hz,
  κ/2π = 128.9 MHz, and 0.93 mW at f = 1, so 8.3 mW at f² = 9
  (printed by `exp_lle.py`);
- Table S3 of the supplement exactly (-6.84, -7.35, -7.89, -8.18,
  -8.41 dB; adiabatic -8.49 dB), and at the design point -7.886 dB,
  bandwidth 1.81 GHz, largest log-negativity 0.133 (printed by
  `exp_tworing.py`);
- Table S5 exactly (-7.274 to -7.284 dB) (printed by `exp_quantum.py`);
- Tables S2 and S4 exactly (printed by `convergence_checks.py`);
- figure previews: `fig0_abstract`, `fig1_device`, `fig3_quantum` and
  `fig4_molecule` came out pixel-identical to the PNG files in `figures/`.
  `fig2_comb` and `fig5_d3` differed in fewer than 0.3% of pixels (in
  `fig2_comb`, the single "lost" point at ζ0 = 3.75 sits at a different
  height);
- the D3 boundary as described in supplement Section 6: crystal at 3.25
  times D3 (residual 6.2e-4), lost at 3.5 times (residual 0.22), odd/even
  contrast -20.6 dB at 4 times (printed by `exp_d3boundary.py`). Its last
  line, the odd-mode growth rate, prints `nan`. This is because the
  odd/even contrast never enters the -160 to -40 dB window the fit uses.
  That agrees with the supplement's statement that the state stays
  symmetric and breathes instead.

**Expected outputs**: the Zenodo archive has reference logs in
`testbench/expected_output/`:

- `exp_lle.log`, `exp_d3boundary.log`, `exp_quantum.log`,
  `exp_addendum.log`, `exp_tworing.log` and `convergence_checks.log`;
- `d3_mechanism.log`: the script that wrote it is in neither this
  repository nor the archive, but `exp_d3mechanism.py` reproduces it (see
  below);
- `d3_gridcheck.log`: no script in this repository (see "Known gaps").

**Check of the fix (30 September 2026).** For this check, I copied the
archived input files from Zenodo v1.1 into `src/`. Then:

- `exp_tworing.py` wrote a `q_tworing.npz` with the same 35 arrays, in the
  same order, as the archived file. Its printed output equals
  `testbench/expected_output/exp_tworing.log` line for line.
  All 35 arrays, including the new `smin_full2_10`, are bit for bit
  identical to the archived ones. (This run took 41 min, against 15 min
  in the earlier check, because the shared machine was heavily loaded.)

The scripts that originally wrote `q_tworing_family.npz` and `d3_mech.npz`
are in neither this repository nor the Zenodo archive.
`exp_tworing_family.py` and `exp_d3mechanism.py` were written for this
fix. Their settings are reconstructions, chosen so that the output matches
the archived files. For `exp_tworing_family.py`, the reconstructed setting
is the choice of which states to include (co-moving residual below 1e-2 and crystal test
passed). For `exp_d3mechanism.py`, the reconstructed settings are the run length (150 time units),
the time step (0.002), the sampling interval (0.25), the second comb line
(μ = 88), and the choice to take the statistics from the second half of the record.
The comparisons below show that these settings reproduce the archived
files.

- `exp_tworing_family.py` reproduced the archived `q_tworing_family.npz`.
  `zeta0`, `kaux` and `kp` were identical, and `full_db` was within
  4.6e-14 dB of the archived values (against a tolerance of 1e-12 dB).
  This run used one BLAS thread (`OPENBLAS_NUM_THREADS=1`), which changes
  rounding in the last digits.
- `exp_d3mechanism.py` reproduced the archived `d3_mech.npz` bit for bit
  (all 5 arrays identical). Its printed output equals
  `testbench/expected_output/d3_mechanism.log` line for line.
- `make_numbers.py` then wrote `numbers/numbers.tex`. Its 73 macros are
  identical to those it writes from the archived files alone. Nothing was
  written outside the repository.
- `fig_tworing.py` ran and drew `figures/fig6_tworing.png` pixel for
  pixel (and byte for byte) identical to the PNG stored in the
  repository. `fig_abstract.py`, which also reads `q_tworing.npz`, still
  gives an identical `fig0_abstract.png`.
- Afterward, I reran the repository's checks. `validate_quantum.py`
  printed the values in the table above. The output of
  `convergence_checks.py` equals `testbench/expected_output/convergence_checks.log`
  line for line.

---

## 10. Notes on the calculations

**Symbols and units.**

- μ is the comb-line (mode) number counted from the pumped line (μ = 0).
- κ is the loaded loss rate (linewidth) of a main-ring mode. Quantum-model
  rates and sideband frequencies ω are given in units of κ.
- D1 is the FSR (as an angular frequency). D2, D3 and D4 describe how the
  resonance spacing changes with μ (dispersion). The integrated dispersion
  is Dint(μ) = D2 μ²/2 + D3 μ³/6 + D4 μ⁴/24.
- The pump detuning is ζ0 = 2δ0/κ, and each mode has ζμ = 2(δ0 + Dint(μ))/κ.
- Soliton time is measured in units of 2/κ. The field is scaled so that
  the Kerr (intensity-dependent) frequency shift equals (κ/2)|ψ|².
- The FEM uses micrometers. Dispersion values are in rad/s inside the
  code and are divided by 2π for printing.
- Squeezing is 10·log10 of the noise relative to vacuum. Spectra are
  averaged over positive and negative sideband frequency, as a homodyne
  detector measures them.
- η is the escape efficiency, the fraction of the odd-mode light that
  reaches the detector. It is 0.5 for the plain critically coupled ring,
  and κP/(κ + κP) with the molecule (`exp_quantum.py`).

**Drifting crystal.** Third-order dispersion makes the crystal drift
around the ring. `exp_lle.py` absorbs the drift speed v into the detuning
(a term -vμ), re-measures the drift and repeats. It does at most 4 rounds
and stops early when the change is below 1e-8. v equals the shift of the
comb spacing, δf_rep = v(κ/2)/2π.

**Which crystals count.** In `exp_lle.py`, a state counts as a 2-FSR
crystal when two things hold. The strongest odd line is more than 30 dB
below the strongest even line. And the field still has a clear peak
(maximum |ψ| above 1.5). States whose remaining time change (residual) is
1e-2 or more are treated as "breathing" (pulsing). They are left out of
the quantum analysis.

**Simplified versus full molecule model.** `exp_quantum.py` and
`exp_addendum.py` treat the small ring through a single extra loss rate
κP on the odd modes. This is called adiabatic elimination: the small ring
is assumed much faster than everything else. `exp_tworing.py` keeps one
explicit small-ring mode per odd mode, coupled at rate J with
κP = 4J²/κaux. It detects at the small ring's drop port.

**Approximations recorded in the code or supplement.** Bending of the
ring is neglected in the mode solving (supplement Section 3). The small
ring uses the same fitted propagation constant as the main ring. The
geometry sweep uses first-order elements for speed. The chosen waveguide
uses second-order elements.

**Hard-coded values.** Some numbers are typed into scripts rather than
read from results:

- `exp_lle.py` stores `boundary_scale=3.5`;
- `make_numbers.py` sets `nDthreeLast = "3.25"` and `nDthreeMuDW = "91"`
  (the 91 is the mode where Dint changes sign at 3.5 times D3, which
  `exp_d3mechanism.py` prints);
- `fig_molecule.py` converts bandwidth with 0.1289 GHz per κ;
- `fig_abstract.py` prints "7.9 dB", "1.81 GHz" and "E_N = 0.13".

**Known gaps** (still open after the fix of 30 September 2026):

1. `figures/figS1.png` has no script in this repository.
2. The D3 boundary on grids of 256 and 320 modes (supplement Section 6) is
   not scripted here. The Zenodo archive has its reference log
   (`testbench/expected_output/d3_gridcheck.log`) but no script for it.
3. `.zenodo.json` still gives an older article title in its description
   ("Photonic-molecule extraction of multimode squeezed light from a
   dispersion-engineered 4H-silicon-carbide soliton-crystal microcomb"),
   not the published title. It also does not list the paper DOI.
4. `CITATION.cff` describes the tagged code (version 1.0.0, 21 July 2026),
   but the code on `main` includes later, untagged changes. Also, the data
   DOI it lists (10.5281/zenodo.21995634, called "Zenodo v1.1" in commit
   47a689f) is not the one the v1.0.0 tag cited (10.5281/zenodo.21471674).
   See [Section 11](#11-version-history).

These gaps were closed on 30 September 2026 (see CHANGELOG.md):

- `fig_tworing.py` stopped with `KeyError: 'smin_full2_10 is not a file in the archive'`
  because `exp_tworing.py` did not save that curve;
- `make_numbers.py` read `q_tworing_family.npz` and `d3_mech.npz`, which
  no script wrote. It also wrote into `../paper/` and `../supplement/`,
  which are not in the repository.

**Differences from the published supplement.** `supplement.pdf` is the
published supplemental document, and it is kept unchanged here. I found
two of its statements that do not match the code and the archived data:

1. *Table S1, first two rows.* With first-order elements, `fem_final.py`
   gives 2.3064460 (2958 elements) at 60 nm and 2.3068171 (3426
   elements) at 50 nm. The archived production output `fem_final.json`
   (Zenodo v1.1) holds the same values (2.306446026374873 with 2958
   elements; 2.3068170818047125 with 3426 elements). A rerun on
   30 September 2026 (gmsh 4.15.2, scikit-fem 12.0.2, femwell 0.1.12)
   reproduced them to about 1e-15. Table S1 instead lists 2.3068171 (3426
   elements) for 60 nm, which is the 50 nm result. For 50 nm it lists
   2.3067888 (3426 elements). The value 2.3067888 is in no archived file.
   None of these settings gave it either: first order with core resolution
   70, 65, 60, 55, 50, 45, 40 or 35 nm; the 50 nm mesh with a 50.06 µm
   bend radius; second order at 60 nm. Table S1 is identical in every
   version of `supplement.pdf` in the git history (first added in commit
   b4f81a4, 22 July 2026). `fem_final.py` has not changed since it was
   added (commit 1ce6f1a, 21 July 2026). So the two rows look like a
   mistake in the table, not in the code. No other number depends on
   them. The other five rows match, and the scripts use only the
   second-order results (the production setting and `nNeffConv` in
   `make_numbers.py`). The code was not changed.
2. *Breathing frequency (supplement Section 6: "58 MHz").*
   `exp_d3mechanism.py` samples the peak field every 0.25 time units for
   150 units. The time unit is 2/κ (`lle.py`). The dominant frequency
   `fdom` = 0.9067 is read from `numpy.fft.rfftfreq`, so it is in cycles
   per time unit. That means one oscillation every 1/0.9067 = 1.10 time
   units. The archived series also shows this directly. In the second
   half, upward zero crossings of the peak field are on average 4.4
   samples (that is, 1.10 time units) apart. `make_numbers.py`
   (`nBreathMHz`) and the archived `d3_mechanism.log` convert it as
   fdom·(κ/2)/2π = 58.5 MHz. This treats `fdom` as an angular frequency.
   Read as cycles per time unit, the same data give fdom·κ/2 ≈ 367 MHz
   (about 2.85 κ/2π), 2π times higher. The conversion was left as it is,
   so that `make_numbers.py` and the new script reproduce the published
   value and the archived log. We, the authors, should confirm which one
   is intended. The peak-field oscillation of about ±17% in the same
   paragraph agrees with `d3_mech.npz` (0.6665 around a mean of 3.979).
   The paragraph also says the odd/even contrast "stays below -260 dB".
   The archived `d3_boundary.npz` track of the same setup (3.25 crystal
   switched to 3.5, up to t = 100) ranges from -281 to -214 dB. So it is
   far below the even modes, but not always below -260 dB.

---

## 11. Version history

| Version | Date | What |
|---|---|---|
| fixes (on `main`, not tagged) | 30 Sep 2026 | `exp_tworing.py` saves the 2 GHz curve that `fig_tworing.py` needs; new `exp_tworing_family.py` and `exp_d3mechanism.py` write the two files `make_numbers.py` needs; `make_numbers.py` writes to `numbers/` instead of `../paper/` and `../supplement/` |
| documentation update (on `main`, not tagged) | 30 Sep 2026 | This guide rewritten; CHANGELOG added |
| untagged changes on `main` | 18 Aug to 4 Sep 2026 | Revision 1 of the article: full two-ring model, small-ring alignment, D3 boundary scan, updated figures and supplement, data DOI changed to 10.5281/zenodo.21995634; then the published citation |
| **v1.0.0** (git tag) | 21 Jul 2026 | First release of the code, figures and metadata |

You can find the details in [CHANGELOG.md](CHANGELOG.md).

**Data DOIs.** At the `v1.0.0` tag, the README and `CITATION.cff` cited
the Zenodo archive 10.5281/zenodo.21471674. Commit 47a689f (18 Aug 2026)
replaced it everywhere with 10.5281/zenodo.21995634, which the commit
message calls "Zenodo v1.1". This guide and `CITATION.cff` cite
10.5281/zenodo.21995634. That is the version that goes with the current
code and the published article. `CITATION.cff` keeps version 1.0.0 and
date 21 July 2026, the last tagged release. No later version number has
been assigned in the repository.
