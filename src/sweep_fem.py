"""Sweep waveguide geometries: neff(lambda) tables for dispersion engineering."""
import numpy as np, json, time
from fem_modes import make_mesh, solve_neff

lams = np.linspace(1.30, 1.85, 13)
heights = [0.50, 0.60]
widths = [1.40, 1.60, 1.85, 2.10, 2.40]

results = {}
t00 = time.time()
for h in heights:
    for w in widths:
        mesh = make_mesh(w, h, res_core=0.035)
        neffs = []
        for lam in lams:
            neff, _ = solve_neff(mesh, w, h, lam, order=1)
            neffs.append(neff)
        key = f"w{w:.2f}_h{h:.2f}"
        results[key] = dict(width=w, height=h, lams=list(lams), neff=neffs)
        print(f"{key}: neff(1.55)~{np.interp(1.55, lams, neffs):.5f} "
              f"[{time.time()-t00:.0f}s]", flush=True)

with open("fem_sweep.json", "w") as f:
    json.dump(results, f, indent=1)
print("DONE", time.time() - t00)
