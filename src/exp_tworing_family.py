"""Full two-ring model across the detuning family at the design point.

Runs the explicit two-ring model of exp_tworing.py (perfect alignment,
detection at the auxiliary drop port) at the design point
kappa_P = 10 kappa, kappa_aux/2pi = 8 GHz, for every stationary 2-FSR
crystal of lle_family.npz (co-moving residual < 1e-2 and crystal test
passed), and records the peak drop-port squeezing of each state.
At zeta0 = 6.5 this is the design-point value of exp_tworing.py.

Needs lle_family.npz and aux_ring.npz. Outputs q_tworing_family.npz
(used by make_numbers.py).
"""
import numpy as np
from exp_tworing import build_tworing, peak, DETECT_AUX, dat, aux, kappa

KP = 10.0                                   # kappa_P / kappa
KAUX = 4 * (float(aux["kaux"]) / kappa)     # 8 GHz design value, kappa units

if __name__ == "__main__":
    stat = (dat["resid"] < 1e-2) & dat["crystal"]
    zeta0s = dat["zeta0"][stat]
    psis = dat["psi"][stat]
    vs = dat["v"][stat]
    J = np.sqrt(KP * KAUX) / 2
    print(f"kappa_P = {KP:.0f} kappa, kappa_aux = {KAUX:.2f} kappa "
          f"({KAUX*kappa/2/np.pi/1e9:.0f} GHz), {len(zeta0s)} states")
    full_db = []
    for z, psi, v in zip(zeta0s, psis, vs):
        cq = build_tworing(psi, float(z), J=J, kaux=KAUX, v=float(v))
        _, _, om, s = peak(cq, DETECT_AUX)
        full_db.append(10 * np.log10(s))
        print(f"zeta0={z:.2f}: full model {full_db[-1]:+.3f} dB "
              f"at omega={om:.2f}", flush=True)
    np.savez("q_tworing_family.npz", zeta0=zeta0s, full_db=full_db,
             kaux=KAUX, kp=KP)
    print("DONE")
