"""Experiment 1: 2-FSR soliton-crystal states of the FEM-designed resonator.

Builds zeta_mu from the FEM dispersion (D2, D3, D4), seeds a two-soliton
crystal, converges it, and continues it adiabatically across the pump
detuning range to map the existence window. Also computes a single-soliton
reference and a D3-scaling family for the higher-order-dispersion study.
"""
import numpy as np, json
from lle import LLE, soliton_crystal_ansatz
from dispersion import RingDispersion

C = 299792458.0
HBAR = 1.054571817e-34

# ---- physical parameters
fem = json.load(open("fem_final.json"))
d = fem["final"]
rd = RingDispersion(d["lams"], d["neff"])
R_UM = rd.radius_for_fsr(350e9)
res = rd.resonances(R_UM, mu_max=80)
D1, D2, D3, D4 = res["D1"], res["D2"], res["D3"], res["D4"]

lam0 = 1.55e-6
om0 = 2 * np.pi * C / lam0
Q_LOAD = 1.5e6
kappa = om0 / Q_LOAD
n2 = 6.9e-19            # m^2/W, 4H-SiC (Guidry et al.)
n0 = 2.5644
Aeff = fem["aeff_um2"] * 1e-12
L = 2 * np.pi * R_UM * 1e-6
Veff = Aeff * L
g0 = HBAR * om0 ** 2 * C * n2 / (n0 ** 2 * Veff)
P1mW = kappa ** 2 * HBAR * om0 / (8 * g0 * 0.5) * 1e3  # power at f=1, mW

print(f"R={R_UM:.2f}um D1/2pi={D1/2/np.pi/1e9:.2f}GHz "
      f"D2/2pi={D2/2/np.pi/1e6:.3f}MHz D3/2pi={D3/2/np.pi/1e3:.2f}kHz "
      f"D4/2pi={D4/2/np.pi:.2f}Hz")
print(f"kappa/2pi={kappa/2/np.pi/1e6:.1f}MHz g0/2pi={g0/2/np.pi:.2f}Hz "
      f"Aeff={fem['aeff_um2']:.3f}um2 P(f=1)={P1mW:.2f}mW")

N = 192
F_PUMP = 3.0


def zeta_of(mu, d3scale=1.0, d4scale=1.0):
    Dint = D2 * mu ** 2 / 2 + d3scale * D3 * mu ** 3 / 6 \
        + d4scale * D4 * mu ** 4 / 24
    return 2 * Dint / kappa


def converge_state(zeta0, psi_mu=None, d3scale=1.0, T=300, dt=0.002,
                   v0=0.0, n_iter=4):
    """Converge to a state that is stationary in a co-moving frame.
    Iteratively absorbs the drift velocity v into a linear detuning term
    (repetition-rate shift of the comb grid). Returns psi, co-moving
    residual, total v, sim (with effective zeta)."""
    from lle import measure_drift
    mu = np.fft.fftfreq(N, d=1.0 / N).astype(int)
    v = v0
    if psi_mu is None:
        d2n = 2 * D2 / kappa
        psi_mu = soliton_crystal_ansatz(N, zeta0, d2n, F_PUMP, n_solitons=2)
    for it in range(n_iter):
        zeta = zeta0 + zeta_of(mu, d3scale) - v * mu
        sim = LLE(N, zeta, F_PUMP)
        psi_mu = sim.run(psi_mu, T if it == 0 else T / 2, dt)
        dv = measure_drift(sim, psi_mu, dt=dt)
        v += dv
        if abs(dv) < 1e-8:
            break
    zeta = zeta0 + zeta_of(mu, d3scale) - v * mu
    sim = LLE(N, zeta, F_PUMP)
    return psi_mu, sim.residual(psi_mu), v, sim


def is_crystal(psi_mu, tol_db=-30):
    """2-FSR crystal: odd-mode power far below even-mode power."""
    mu = np.fft.fftfreq(N, d=1.0 / N).astype(int)
    p = np.abs(psi_mu) ** 2
    even = p[(mu % 2 == 0) & (np.abs(mu) > 0) & (np.abs(mu) < 40)].max()
    odd = p[(mu % 2 != 0) & (np.abs(mu) < 40)].max()
    return 10 * np.log10(odd / even) < tol_db, 10 * np.log10(odd / even)


def is_soliton_state(psi_mu, thresh=1.5):
    """Any localized structure present (peak well above background)?"""
    th = np.abs(np.fft.ifft(psi_mu)) * N ** 0.5
    return th.max() > thresh


def continuation(zeta_start, zeta_stop, dz, d3scale=1.0, label=""):
    states = []
    psi, r, v, _ = converge_state(zeta_start, d3scale=d3scale, T=400)
    z = zeta_start
    while (dz > 0 and z <= zeta_stop) or (dz < 0 and z >= zeta_stop):
        ok, ratio = is_crystal(psi)
        alive = is_soliton_state(psi)
        states.append(dict(zeta0=z, psi=psi.copy(), resid=r, v=v,
                           crystal=ok and alive, odd_even_db=ratio))
        print(f"{label} zeta0={z:.2f} resid={r:.2e} v={v:.3e} "
              f"crystal={ok and alive} odd/even={ratio:.1f}dB", flush=True)
        if not (ok and alive) and len(states) > 1:
            break
        z += dz
        psi, r, v, _ = converge_state(z, psi_mu=psi, d3scale=d3scale,
                                      T=150, v0=v)
    return states


if __name__ == "__main__":
    all_out = {}
    # main family (FEM dispersion as computed)
    up = continuation(6.0, 12.0, 0.25, label="up")
    dn = continuation(5.75, 2.0, -0.25, label="dn")
    fam = sorted(dn[::-1] + up, key=lambda s: s["zeta0"])
    np.savez("lle_family.npz",
             zeta0=[s["zeta0"] for s in fam],
             crystal=[s["crystal"] for s in fam],
             resid=[s["resid"] for s in fam],
             v=[s["v"] for s in fam],
             odd_even_db=[s["odd_even_db"] for s in fam],
             psi=np.array([s["psi"] for s in fam]),
             D=[D1, D2, D3, D4], kappa=kappa, g0=g0, f=F_PUMP, N=N,
             R_um=R_UM, P1mW=P1mW, Aeff_um2=fem["aeff_um2"])
    # D3-scaling family at fixed representative detuning
    z_rep = 6.5
    d3_states = {}
    d3_v = {}
    keys = ["0.0", "1.0", "2.0", "2.5", "3.0"]
    for s3 in [float(k) for k in keys]:
        psi, r, v, _ = converge_state(z_rep, d3scale=s3, T=500)
        ok, ratio = is_crystal(psi)
        d3_states[str(s3)] = psi
        d3_v[str(s3)] = v
        print(f"D3x{s3}: resid={r:.2e} v={v:.3e} crystal={ok} "
              f"odd/even={ratio:.1f}dB", flush=True)
    # the stability boundary itself (crystal lost between 3.25 and 3.5,
    # scan step 0.25) is located by exp_d3boundary.py
    np.savez("lle_d3family.npz", zeta0=z_rep, keys=keys,
             boundary_scale=3.5,
             v=[d3_v[k] for k in keys],
             **{f"psi_{k}": v for k, v in d3_states.items()})
    print("DONE")
