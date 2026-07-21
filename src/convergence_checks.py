"""Numerical convergence checks reported in the supplementary document:
(i) LLE grid size N and time step dt; (ii) quantum mode truncation M;
(iii) frequency-grid resolution for peak squeezing."""
import numpy as np
import exp_lle as X
from exp_quantum import build, max_squeezing, ODD
import exp_quantum as Q

z = 6.5
print("--- LLE dt/N convergence (comb line mu=2 amplitude, peak |psi|)")
ref = None
for (Ngrid, dt) in [(192, 0.004), (192, 0.002), (192, 0.001), (256, 0.002)]:
    X.N = Ngrid
    psi, r, v, _ = X.converge_state(z, T=400, dt=dt)
    th = np.abs(np.fft.ifft(psi)) * Ngrid ** 0.5
    a2 = np.abs(psi[2]) / Ngrid ** 0.5
    print(f"N={Ngrid} dt={dt}: |phi_2|={a2:.6f} max|psi|={th.max():.6f} "
          f"resid={r:.1e} v={v:.6e}", flush=True)
X.N = 192

print("--- quantum truncation M (bare, critical coupling)")
dat = np.load("lle_family.npz")
iz = int(np.argmin(np.abs(dat["zeta0"] - z)))
psi = dat["psi"][iz]
vrep = float(dat["v"][iz])
zz = float(dat["zeta0"][iz])
for Mtr in (20, 30, 40):
    Q.M = Mtr
    Q.ODD = [i for i, m in enumerate(range(-Mtr, Mtr + 1)) if m % 2 != 0]
    cq = Q.build(psi, zz, v=vrep)
    oms, smin, smax, om_at, s_at = Q.max_squeezing(cq, om_max=8, n_om=161)
    print(f"M={Mtr}: min={10*np.log10(s_at):.4f} dB at omega={om_at:.3f}",
          flush=True)
    # molecule case
    cqm = Q.build(psi, zz, kc_odd=20, v=vrep)
    _, _, _, om2, s2 = Q.max_squeezing(cqm, om_max=60, n_om=481)
    print(f"M={Mtr} molecule kP=20: min={10*np.log10(s2):.4f} dB "
          f"at omega={om2:.3f}", flush=True)
