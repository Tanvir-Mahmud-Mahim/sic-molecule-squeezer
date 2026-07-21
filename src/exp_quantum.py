"""Experiment 2: multimode squeezing and entanglement of the below-threshold
odd-mode lattice of the 2-FSR soliton crystal, for the bare resonator and
the photonic-molecule (Purcell-outcoupled) source.

All rates in units of the loaded linewidth kappa.
"""
import numpy as np
from quantum import CombQuantum, log_negativity

dat = np.load("lle_family.npz")
D1, D2, D3, D4 = dat["D"]
kappa = float(dat["kappa"])
N = int(dat["N"])
M = 30          # fluctuation modes mu = -M..M
Mc = 60         # classical comb modes kept

mu_fft = np.fft.fftfreq(N, d=1.0 / N).astype(int)


def centered_psi(psi_fft):
    """Convert fft-ordered LLE modal amplitudes to centered theta-normalized
    amplitudes phi_mu = psi_mu/sqrt(N), for which g0 A_m A_n = (kappa/2)
    phi_m phi_n (see lle.py normalization)."""
    out = np.zeros(2 * Mc + 1, complex)
    for i, m in enumerate(mu_fft):
        if -Mc <= m <= Mc:
            out[m + Mc] = psi_fft[i] / np.sqrt(N)
    return out


def zeta_vec(zeta0, d3scale=1.0, v=0.0):
    """Detunings in the co-moving comb frame (v = repetition-rate shift)."""
    mu = np.arange(-M, M + 1)
    Dint = D2 * mu ** 2 / 2 + d3scale * D3 * mu ** 3 / 6 + D4 * mu ** 4 / 24
    return zeta0 + 2 * Dint / kappa - v * mu


ODD = [i for i, m in enumerate(range(-M, M + 1)) if m % 2 != 0]


def build(psi_fft, zeta0, kc_odd=0.5, d3scale=1.0, v=0.0):
    """kc_odd: collection rate for odd modes (kappa units).
    Bare critical-coupled device: kc_odd=0.5 (bus port monitored, ki=0.5).
    Photonic molecule: odd modes get Purcell rate kc_odd=kP to the drop
    port; the intrinsic+bus loss (=1) is then unmonitored."""
    n = 2 * M + 1
    mu = np.arange(-M, M + 1)
    ki = np.full(n, 0.5)
    kc = np.full(n, 0.5)
    if kc_odd != 0.5:
        odd = mu % 2 != 0
        ki[odd] = 1.0          # full bare loss unmonitored
        kc[odd] = kc_odd       # Purcell drop-port rate
    return CombQuantum(centered_psi(psi_fft), zeta_vec(zeta0, d3scale, v),
                       ki, kc, M)


def max_squeezing(cq, om_max=30.0, n_om=121):
    oms = np.linspace(0, om_max, n_om)
    smin, smax = cq.squeezing_spectrum(oms, ODD)
    i = np.argmin(smin)
    return oms, smin, smax, oms[i], smin[i]


if __name__ == "__main__":
    zeta0s = dat["zeta0"]
    crystal = dat["crystal"]
    psis = dat["psi"]
    vs = dat["v"]
    resids = dat["resid"]
    stationary = resids < 1e-2   # exclude breathing states
    out = {}

    # ---- A: sweep across existence range, bare critical coupling
    rows = []
    for z, ok, psi, v, st in zip(zeta0s, crystal, psis, vs, stationary):
        if not (ok and st):
            continue
        cq = build(psi, z, v=v)
        # stability check
        maxRe = float(np.max(np.real(np.linalg.eigvals(cq.Mdrift))))
        oms, smin, smax, om_at, s_at = max_squeezing(cq, om_max=8, n_om=81)
        # intracavity (Lyapunov) minimum variance
        V = cq.steady_covariance()
        # restrict to odd modes' quadratures
        n = 2 * M + 1
        idx = ODD + [n + i for i in ODD]
        w = np.linalg.eigvalsh(V[np.ix_(idx, idx)])
        rows.append(dict(zeta0=float(z), maxRe=maxRe,
                         out_min_db=10 * np.log10(s_at),
                         out_om=float(om_at),
                         intra_min_db=10 * np.log10(w[0]),
                         intra_max_db=10 * np.log10(w[-1])))
        print(rows[-1], flush=True)
    np.savez("q_sweep.npz", rows=rows)

    # ---- B: representative state: spectra, supermodes, molecule sweep
    iz = int(np.argmin(np.abs(zeta0s - 6.5)))
    z, psi, vrep = float(zeta0s[iz]), psis[iz], float(vs[iz])
    print("representative zeta0 =", z)
    # bare
    cq = build(psi, z, v=vrep)
    oms = np.linspace(0, 40, 161)
    smin_b, smax_b = cq.squeezing_spectrum(oms, ODD)
    # supermodes at optimal freq (bare)
    ib = np.argmin(smin_b)
    wb, Ub = cq.supermodes(oms[ib], ODD, k=6)
    # molecule sweep
    kps = np.array([0.5, 1, 2, 5, 10, 20, 30, 50])
    mol = []
    spectra = {}
    for kp in kps:
        cqm = build(psi, z, kc_odd=kp, v=vrep)
        omsm, smin, smax, om_at, s_at = max_squeezing(cqm, om_max=60, n_om=241)
        # squeezing bandwidth: full width where smin < (1+s_at)/2 roughly
        thr = 1 - (1 - s_at) / 2
        bw = omsm[smin < thr]
        bw = bw.max() - bw.min() if len(bw) > 1 else 0.0
        eta = kp / (1.0 + kp) if kp != 0.5 else 0.5
        mol.append(dict(kp=float(kp), min_db=10 * np.log10(s_at),
                        om_at=float(om_at), bw=float(bw), eta=float(eta)))
        spectra[f"smin_{kp}"] = smin
        spectra[f"smax_{kp}"] = smax
        print(mol[-1], flush=True)
    # supermodes for kp=20 at optimal
    cqm = build(psi, z, kc_odd=20, v=vrep)
    omsm, sminm, smaxm, om_atm, s_atm = max_squeezing(cqm, om_max=60, n_om=241)
    wm, Um = cqm.supermodes(om_atm, ODD, k=6)
    np.savez("q_rep.npz", zeta0=z, oms=oms, smin_bare=smin_b, smax_bare=smax_b,
             omsm=np.linspace(0, 60, 241), sup_bare_w=wb, sup_bare_U=Ub,
             sup_mol_w=wm, sup_mol_U=Um, om_opt_bare=oms[ib],
             om_opt_mol=om_atm, kps=kps,
             mol_min_db=[m["min_db"] for m in mol],
             mol_bw=[m["bw"] for m in mol],
             mol_eta=[m["eta"] for m in mol], **spectra)

    # ---- C: entanglement matrix at drop port (kp=20), zero sideband
    Vout = cqm.quad_covariance_sym(om_atm, ODD)
    d = len(ODD)
    EN = np.zeros((d, d))
    for a in range(d):
        for b in range(a + 1, d):
            EN[a, b] = EN[b, a] = log_negativity(Vout, a, b, d)
    np.savez("q_entanglement.npz", EN=EN,
             odd_mu=[m for m in range(-M, M + 1) if m % 2 != 0])
    print("max EN:", EN.max())

    # ---- D: D3 family
    d3 = np.load("lle_d3family.npz")
    z3 = float(d3["zeta0"])
    v3s = d3["v"]
    d3keys = [str(k) for k in d3["keys"]]
    d3rows = []
    d3spec = {}
    for js3, sk in enumerate(d3keys):
        s3 = float(sk)
        psi3 = d3[f"psi_{sk}"]
        cq3 = build(psi3, z3, kc_odd=20, d3scale=s3, v=float(v3s[js3]))
        oms3, smin, smax, om_at, s_at = max_squeezing(cq3, om_max=60, n_om=241)
        w3, U3 = cq3.supermodes(om_at, ODD, k=4)
        d3rows.append(dict(s3=s3, min_db=10 * np.log10(s_at),
                           om_at=float(om_at),
                           second_db=10 * np.log10(w3[1])))
        d3spec[f"smin_{s3}"] = smin
        d3spec[f"U_{s3}"] = U3
        d3spec[f"w_{s3}"] = w3
        print(d3rows[-1], flush=True)
    np.savez("q_d3.npz", oms=np.linspace(0, 60, 241),
             s3=[r["s3"] for r in d3rows],
             min_db=[r["min_db"] for r in d3rows],
             om_at=[r["om_at"] for r in d3rows],
             second_db=[r["second_db"] for r in d3rows], **d3spec)
    print("DONE")
