"""Addendum: (a) ideal-detection (eta=1) squeezing across the family;
(b) molecule squeezing (kP=10, 20) across the detuning family."""
import numpy as np
import exp_quantum as Q
from exp_quantum import build, max_squeezing, ODD
from quantum import CombQuantum

dat = np.load("lle_family.npz")
rows = []
for z, ok, psi, v, r in zip(dat["zeta0"], dat["crystal"], dat["psi"],
                            dat["v"], dat["resid"]):
    if not (ok and r < 1e-2):
        continue
    row = dict(zeta0=float(z))
    # ideal: all loss monitored (eta -> 1), same total linewidth
    n = 2 * Q.M + 1
    cqi = CombQuantum(Q.centered_psi(psi), Q.zeta_vec(float(z), v=float(v)),
                      np.full(n, 1e-6), np.full(n, 1.0), Q.M)
    _, _, _, om_i, s_i = max_squeezing(cqi, om_max=10, n_om=161)
    row["ideal_db"] = 10 * np.log10(s_i)
    for kp in (10.0, 20.0):
        cqm = build(psi, float(z), kc_odd=kp, v=float(v))
        _, _, _, om_m, s_m = max_squeezing(cqm, om_max=40, n_om=161)
        row[f"mol{int(kp)}_db"] = 10 * np.log10(s_m)
    rows.append(row)
    print(row, flush=True)
np.savez("q_addendum.npz",
         zeta0=[r["zeta0"] for r in rows],
         ideal_db=[r["ideal_db"] for r in rows],
         mol10_db=[r["mol10_db"] for r in rows],
         mol20_db=[r["mol20_db"] for r in rows])
print("DONE")
