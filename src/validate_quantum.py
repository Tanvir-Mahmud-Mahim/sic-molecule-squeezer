"""Validation of the linearized quantum module against analytic results.

Test 1: vacuum (no comb) -> V(omega) = identity.
Test 2: cw-pumped Kerr resonator: drift-matrix eigenvalues for the (mu,-mu)
        pair must match the analytic Bogoliubov result
        lambda = -kt/2 +/- sqrt(|F|^2 - dbar^2), F = (1/2)psi0^2,
        dbar = Delta_mu + |psi0|^2 (kappa units, symmetric zeta).
Test 3: output squeezing of the critically coupled cw Kerr OPO below
        threshold is bounded by 3 dB (escape efficiency 1/2), and grows
        towards the bound at threshold; overcoupled case exceeds it.
"""
import numpy as np
from quantum import CombQuantum
from lle import cw_background, LLE

M = 6
n = 2 * M + 1

# --- Test 1: vacuum
psi = np.zeros(2 * M + 1, complex)
cq = CombQuantum(psi, np.zeros(n), 0.5 * np.ones(n), 0.5 * np.ones(n), M)
V = cq.quad_covariance_sym(0.7, list(range(n)))
print("T1 vacuum: max|V - I| =", np.max(np.abs(V - np.eye(2 * n))))

# --- Test 2: cw Kerr Bogoliubov eigenvalues
zeta0, f = 1.2, 1.1
psi0 = cw_background(zeta0, f)
rho2 = np.abs(psi0) ** 2
zeta = np.full(n, zeta0)   # flat dispersion for the test
psic = np.zeros(n, complex)
psic[M] = psi0
cq = CombQuantum(psic, zeta, 0.5 * np.ones(n), 0.5 * np.ones(n), M)
ev = np.linalg.eigvals(cq.Mdrift)
mu_test = 3
Delta = 0.5 * zeta0
dbar = -Delta + 2 * rho2   # effective detuning of fluctuation (sign conv.)
gain = 0.5 * rho2
lam_an = -0.5 + np.sqrt(complex(gain ** 2 - (Delta - rho2) ** 2))
print("T2 analytic Re(lam) =", np.real(lam_an),
      " max Re(numeric) =", np.max(np.real(ev)))

# --- Test 3: squeezing bounds vs coupling for near-threshold cw state
# threshold for pair (mu,-mu): gain = |dbar| happens at ... just scan pump
for eta, kc in (("critical", 0.5), ("overcoupled 10x", 10 * 0.5)):
    ki = 0.5
    kt = ki + kc
    # rescale rates so total kappa = kt (units still fine: kappa ref = 1)
    cqa = CombQuantum(psic, zeta * (1.0 / kt), np.full(n, ki) / kt,
                      np.full(n, kc) / kt, M)
    om = np.linspace(0, 3, 61)
    smin, smax = cqa.squeezing_spectrum(om, list(range(n)))
    print(f"T3 {eta}: min V = {smin.min():.4f} "
          f"({10*np.log10(smin.min()):.2f} dB), "
          f"uncertainty check minV*maxV>=1: {(smin*smax).min():.3f}")
