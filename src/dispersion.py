"""Ring-resonator dispersion from FEM neff(lambda) tables.

Fits neff(omega) with a polynomial, builds the resonance grid
beta(omega_m) * 2*pi*R = 2*pi*m, and extracts D1, D2, D3, Dint.
"""
import numpy as np
from scipy.optimize import brentq

C = 299792458.0  # m/s


class RingDispersion:
    def __init__(self, lams_um, neffs, deg=6):
        lams = np.asarray(lams_um) * 1e-6
        self.omega_grid = 2 * np.pi * C / lams
        self.neff_grid = np.asarray(neffs)
        # fit neff in normalized omega for conditioning
        self.om0n = self.omega_grid.mean()
        self.oms = self.omega_grid.std()
        x = (self.omega_grid - self.om0n) / self.oms
        self.pfit = np.polyfit(x, self.neff_grid, deg)
        self.resid = np.max(np.abs(np.polyval(self.pfit, x) - self.neff_grid))

    def neff(self, omega):
        return np.polyval(self.pfit, (omega - self.om0n) / self.oms)

    def beta(self, omega):
        return self.neff(omega) * omega / C

    def ng(self, omega, do=1e9):
        return C * (self.beta(omega + do) - self.beta(omega - do)) / (2 * do)

    def resonances(self, radius_um, lam0_um=1.55, mu_max=80):
        """Return mu array, omega_mu, and D1, D2, D3, D4 (angular, rad/s)."""
        L = 2 * np.pi * radius_um * 1e-6
        om_target = 2 * np.pi * C / (lam0_um * 1e-6)
        m0 = int(round(self.beta(om_target) * L / (2 * np.pi)))
        om_lo, om_hi = self.omega_grid.min(), self.omega_grid.max()

        def f(om, m):
            return self.beta(om) * L - 2 * np.pi * m

        mus, oms = [], []
        for mu in range(-mu_max, mu_max + 1):
            m = m0 + mu
            try:
                om = brentq(f, om_lo, om_hi, args=(m,), xtol=1e-3)
            except ValueError:
                continue
            mus.append(mu)
            oms.append(om)
        mus = np.array(mus)
        oms = np.array(oms)
        # fit local expansion omega_mu = om0 + D1 mu + D2 mu^2/2 + D3 mu^3/6 + D4 mu^4/24
        sel = np.abs(mus) <= min(40, mu_max)
        pf = np.polyfit(mus[sel], oms[sel], 4)
        D1 = pf[3]
        D2 = 2 * pf[2]
        D3 = 6 * pf[1]
        D4 = 24 * pf[0]
        om0 = pf[4]
        Dint = oms - (om0 + D1 * mus)
        return dict(mu=mus, omega=oms, om0=om0, D1=D1, D2=D2, D3=D3, D4=D4,
                    Dint=Dint, m0=m0)

    def radius_for_fsr(self, fsr_hz, lam0_um=1.55):
        om0 = 2 * np.pi * C / (lam0_um * 1e-6)
        ng = self.ng(om0)
        return C / (2 * np.pi * ng * fsr_hz) * 1e6  # um
