"""High-accuracy FEM run for the selected geometry (w=1.85, h=0.50 um):
order-2 elements, denser wavelength grid, effective area, field export,
and a mesh-convergence table."""
import numpy as np, json, time
from fem_modes import make_mesh, solve_neff

W, H = 1.85, 0.50
lams = np.linspace(1.30, 1.85, 17)

out = {}
t0 = time.time()

# --- convergence table at 1.55 um
conv = []
for res, order in [(0.06, 1), (0.05, 1), (0.04, 1), (0.03, 1),
                   (0.05, 2), (0.04, 2), (0.03, 2)]:
    mesh = make_mesh(W, H, res_core=res)
    neff, mode = solve_neff(mesh, W, H, 1.55, order=order)
    conv.append(dict(res=res, order=order, nelem=int(mesh.t.shape[1]),
                     neff=neff))
    print(f"conv res={res} order={order} neff={neff:.7f} [{time.time()-t0:.0f}s]",
          flush=True)
out["convergence"] = conv

# --- final sweep, order 2
mesh = make_mesh(W, H, res_core=0.04)
neffs = []
for lam in lams:
    neff, mode = solve_neff(mesh, W, H, lam, order=2)
    neffs.append(neff)
    print(f"lam={lam:.3f} neff={neff:.7f} [{time.time()-t0:.0f}s]", flush=True)
out["final"] = dict(width=W, height=H, lams=list(lams), neff=neffs)

# --- effective area and field at 1.55
neff, mode = solve_neff(mesh, W, H, 1.55, order=2)
try:
    aeff = float(np.real(mode.calculate_effective_area()))
except Exception as e:
    aeff = None
    print("aeff failed:", e)
out["aeff_um2"] = aeff
print("Aeff =", aeff, "um^2")

# export |E| on a grid for the figure
from skfem import Basis, ElementTriP0
xs = np.linspace(-2.2, 2.2, 221)
ys = np.linspace(-1.4, 1.9, 166)
Xg, Yg = np.meshgrid(xs, ys)
pts = np.vstack([Xg.ravel(), Yg.ravel()])
basis = mode.basis
(et, et_basis), (ez, ez_basis) = mode.basis.split(mode.E)
from skfem import ElementVector
probes = et_basis.probes(pts)
Et = probes @ et
Emag2 = np.abs(Et[:len(pts[0])])**2 + np.abs(Et[len(pts[0]):])**2
np.savez("mode_field.npz", x=xs, y=ys,
         E2=Emag2.reshape(len(ys), len(xs)))
with open("fem_final.json", "w") as f:
    json.dump(out, f, indent=1)
print("DONE", time.time() - t0)
