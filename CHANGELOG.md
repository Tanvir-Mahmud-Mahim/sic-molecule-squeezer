# Changelog

All notable changes to this repository are listed here, newest first.
Entries before the documentation update are taken from the git history and
the `v1.0.0` tag.

## Documentation update (30 September 2026, branch `docs`, not tagged)

Documentation only; no code, data or figure was changed.

- README rewritten as a step-by-step guide: plain-language summary, file
  tree, installation notes (including the system libraries `gmsh` needs),
  three ways to use the code, a table of scripts with inputs, outputs and
  run times, a table of which script draws which figure file, parameter
  sources, built-in checks, notes on units and conventions, and a list of
  known gaps found while checking.
- Added this CHANGELOG.
- `CITATION.cff`: data-DOI descriptions now say which archive version
  each DOI is (10.5281/zenodo.21995634 is "Zenodo v1.1" per commit
  47a689f; 10.5281/zenodo.21471674 was cited at tag v1.0.0). Version
  (1.0.0) and date unchanged.

## Untagged changes on `main` after v1.0.0

### 4 September 2026

- Article published: `README.md` and `CITATION.cff` now cite
  Optics Express 34(18), 34822-34834 (2026), doi:10.1364/OE.612248, and a
  paper badge was added.

### 18 August 2026 (Revision 1 of the article)

- New scripts: `aux_ring.py` (alignment of the two rings, radius-error
  tolerance), `exp_tworing.py` (full two-ring quantum model),
  `exp_d3boundary.py` (third-order-dispersion stability boundary, step
  0.25) and `fig_tworing.py` (new figure `fig6_tworing`).
- `exp_lle.py` and `make_numbers.py` updated for the new results.
- `supplement.pdf` updated, including a corrected reference author order
  (Gouzien et al., Phys. Rev. Res. 5, 023178).
- Dataset links in `README.md`, `CITATION.cff`, `.zenodo.json` and the
  supplement changed from 10.5281/zenodo.21471674 to
  10.5281/zenodo.21995634 (described in the commit message as Zenodo v1.1).
- Figure polish (annotation placement; graphical abstract now shows the
  full-model values) in `fig_abstract.py`, `fig_d3.py`, `fig_quantum.py`,
  `fig_tworing.py`, with updated PNG previews.
- README key-results table updated to the full two-ring model values
  (commit 391f79e).
- `supplement.pdf` re-uploaded once more after the reference fix (commit
  4df24ea; the commit message does not say what changed).

### 22 to 31 July 2026

- `supplement.pdf` added (22 July).
- README edited twice (31 July).

## v1.0.0 (21 July 2026)

First release (git tag `v1.0.0`).

- Repository root: README, Apache-2.0 LICENSE, NOTICE, `CITATION.cff`,
  `.zenodo.json`, `requirements.txt`, `.gitignore`.
- PNG previews of the figures `fig0_abstract` to `fig5_d3` and `figS1`.
- `src/`: material models, FEM mode solving, ring dispersion,
  Lugiato-Lefever solver and experiments, linearized quantum model and
  experiments, validation and convergence scripts, `make_numbers.py`, and
  the figure scripts.
