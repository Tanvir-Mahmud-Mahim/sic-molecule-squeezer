"""Full-vector FEM mode solving of 4H-SiCOI waveguides with femwell (scikit-fem).

Computes neff(lambda) for the fundamental TE mode of an oxide-clad
4H-SiC-on-insulator waveguide, optionally with bend radius.
Units: micrometres.
"""
import numpy as np
from collections import OrderedDict
from shapely.geometry import box
from skfem import Basis, ElementTriP0
from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import compute_modes
from skfem.io.meshio import from_meshio

from materials import n_sic_o, n_sio2

WINDOW_W = 3.2   # half-width of computational window
WINDOW_TOP = 2.2
WINDOW_BOT = 2.2


def make_mesh(width, height, res_core=0.04, res_clad=0.35):
    core = box(-width / 2, 0, width / 2, height)
    clad = box(-WINDOW_W, -WINDOW_BOT, WINDOW_W, height + WINDOW_TOP)
    polygons = OrderedDict(core=core, clad=clad)
    resolutions = dict(
        core={"resolution": res_core, "distance": 0.6},
    )
    mesh = from_meshio(
        mesh_from_OrderedDict(polygons, resolutions,
                              default_resolution_max=res_clad,
                              filename="mesh.msh"))
    return mesh


def solve_neff(mesh, width, height, lam, radius=np.inf, order=1):
    basis0 = Basis(mesh, ElementTriP0())
    eps = basis0.zeros()
    for subdomain, n_func in (("core", n_sic_o), ("clad", n_sio2)):
        eps[basis0.get_dofs(elements=subdomain)] = n_func(lam) ** 2
    modes = compute_modes(basis0, eps, wavelength=lam, num_modes=2,
                          order=order, radius=radius,
                          n_guess=float(n_sic_o(lam)))
    # pick highest-neff TE-polarized mode
    best = None
    for m in modes:
        if m.te_fraction > 0.6:
            best = m
            break
    if best is None:
        best = modes[0]
    return float(np.real(best.n_eff)), best


if __name__ == "__main__":
    import time
    mesh = make_mesh(1.85, 0.5)
    t0 = time.time()
    neff, mode = solve_neff(mesh, 1.85, 0.5, 1.55)
    print(f"neff(1.55um, w=1.85, h=0.5) = {neff:.6f}  "
          f"te_frac={mode.te_fraction:.3f}  t={time.time()-t0:.1f}s  "
          f"nelem={mesh.t.shape[1]}")
