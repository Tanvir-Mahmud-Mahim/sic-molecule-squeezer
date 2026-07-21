"""Material dispersion models.

Sellmeier data from the refractiveindex.info database (CC0, public domain):
- 4H-SiC ordinary/extraordinary: S. Wang et al., Laser Photonics Rev. 7, 831 (2013),
  doi:10.1002/lpor.201300068  (formula 2 / formula 4, 0.4047-5 um)
- SiO2: I. H. Malitson, J. Opt. Soc. Am. 55, 1205 (1965),
  doi:10.1364/JOSA.55.001205  (formula 1, 0.21-6.7 um)
Wavelengths in micrometres.
"""
import numpy as np


def n_sic_o(lam_um):
    """4H-SiC ordinary ray (in-plane E field; TE modes on c-axis-normal SiCOI).
    refractiveindex.info formula 2 with coefficients from Wang 2013."""
    l2 = np.asarray(lam_um, dtype=float) ** 2
    n2 = (1.0 + 0.0
          + 0.20075 * l2 / (l2 - (-12.07224))
          + 5.54861 * l2 / (l2 - 0.02641)
          + 35.65066 * l2 / (l2 - 1268.24708))
    return np.sqrt(n2)


def n_sic_e(lam_um):
    """4H-SiC extraordinary ray (E parallel to c axis). Formula 4, Wang 2013."""
    lam = np.asarray(lam_um, dtype=float)
    l2 = lam ** 2
    n2 = (6.79485
          + 0.15558 * lam ** 0 / (l2 - 0.03535 ** 1)
          + 0.0
          + (-0.02296) * lam ** 2)
    return np.sqrt(n2)


def n_sio2(lam_um):
    """Fused silica, Malitson 1965 Sellmeier (formula 1)."""
    l2 = np.asarray(lam_um, dtype=float) ** 2
    n2 = (1.0
          + 0.6961663 * l2 / (l2 - 0.0684043 ** 2)
          + 0.4079426 * l2 / (l2 - 0.1162414 ** 2)
          + 0.8974794 * l2 / (l2 - 9.896161 ** 2))
    return np.sqrt(n2)


if __name__ == "__main__":
    for lam in (1.30, 1.55, 1.80):
        print(f"lam={lam} um: n_SiC_o={n_sic_o(lam):.4f} "
              f"n_SiC_e={n_sic_e(lam):.4f} n_SiO2={n_sio2(lam):.4f}")
