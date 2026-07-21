"""Lugiato-Lefever simulation in the modal basis with arbitrary Dint.

Normalized units: time in 2/kappa, field psi with g0|A|^2 = (kappa/2)|psi|^2,
detuning zeta_mu = 2*(delta0 + Dint(mu))/kappa, pump
f^2 = 8*g0*eta_c*P_in/(kappa^2*hbar*omega0).
d psi_mu/dt = -(1 + i zeta_mu) psi_mu + i FFT[|psi|^2 psi]_mu + f delta_mu0
"""
import numpy as np


class LLE:
    def __init__(self, N, zeta_mu, f):
        self.N = N
        self.mu = np.fft.fftfreq(N, d=1.0 / N).astype(int)
        self.zeta = np.asarray(zeta_mu, dtype=float)  # fft order
        self.f = f

    def _lin_half(self, psi_mu, dt):
        L = -(1.0 + 1j * self.zeta)
        eL = np.exp(L * dt / 2)
        pump = np.zeros_like(psi_mu)
        pump[0] = self.f * self.N ** 0.5   # uniform pump in theta domain
        return eL * psi_mu + (eL - 1.0) / L * pump

    def step(self, psi_mu, dt):
        psi_mu = self._lin_half(psi_mu, dt)
        psi = np.fft.ifft(psi_mu) * self.N ** 0.5
        psi = psi * np.exp(1j * np.abs(psi) ** 2 * dt)
        psi_mu = np.fft.fft(psi) / self.N ** 0.5
        return self._lin_half(psi_mu, dt)

    def run(self, psi_mu, T, dt=0.005):
        for _ in range(int(T / dt)):
            psi_mu = self.step(psi_mu, dt)
        return psi_mu

    def residual(self, psi_mu):
        """|d psi/dt| for stationarity check."""
        psi = np.fft.ifft(psi_mu) * self.N ** 0.5
        nl = np.fft.fft(np.abs(psi) ** 2 * psi) / self.N ** 0.5
        pump = np.zeros_like(psi_mu)
        pump[0] = self.f * self.N ** 0.5
        d = -(1.0 + 1j * self.zeta) * psi_mu + 1j * nl + pump
        return np.max(np.abs(d)) / self.N ** 0.5


def measure_drift(sim, psi_mu, dt=0.002, T=1.0):
    """Drift velocity v (rad per normalized time) from modal phase slopes:
    for a rigidly drifting state psi_mu(t) = psi_mu(0) exp(-i mu v t)."""
    psi2 = psi_mu.copy()
    for _ in range(int(T / dt)):
        psi2 = sim.step(psi2, dt)
    mu = sim.mu
    sel = (np.abs(psi_mu) > 1e-3 * np.max(np.abs(psi_mu))) & (mu != 0) \
        & (np.abs(mu) < 60)
    dphi = np.angle(psi2[sel] * np.conj(psi_mu[sel]))
    slope = np.sum(mu[sel] * dphi) / np.sum(mu[sel] ** 2)
    return -slope / T


def cw_background(zeta0, f, branch="lower"):
    rho2 = 0.0
    for _ in range(500):
        rho2 = f ** 2 / (1 + (zeta0 - rho2) ** 2)
    return f / (1 + 1j * (zeta0 - rho2))


def soliton_crystal_ansatz(N, zeta0, d2, f, n_solitons=2):
    """Approximate n-soliton crystal on the lower-branch cw background."""
    theta = 2 * np.pi * np.arange(N) / N - np.pi
    psi = np.full(N, cw_background(zeta0, f))
    arg = np.sqrt(8 * zeta0) / (np.pi * f)
    phi0 = np.arccos(min(arg, 1.0))
    width = np.sqrt(d2 / (2 * zeta0)) if d2 > 0 else 0.01
    for k in range(n_solitons):
        c = -np.pi + (2 * k + 1) * np.pi / n_solitons
        d = np.angle(np.exp(1j * (theta - c)))
        psi += np.sqrt(2 * zeta0) / np.cosh(d / width) * np.exp(1j * phi0)
    return np.fft.fft(psi) / N ** 0.5
