"""Full two-ring (non-adiabatic) quantum model of the photonic molecule.

Extends the linearized model of quantum.py with one explicit auxiliary-ring
mode per odd main-ring mode, coupled at rate J, with drop-port linewidth
kappa_aux and per-mode detuning taken from the two-ring resonance mismatch
computed in aux_ring.py. Detection is at the auxiliary drop port.

Runs (representative state zeta0 = 6.5, all rates in units of kappa):
  A. adiabatic reference (uniform kappa_P = 4 J^2 / kappa_aux);
  B. full two-ring model at perfect alignment (Delta_mu = 0): isolates the
     non-adiabatic error of the eliminated auxiliary ring (Reviewer 4.1);
  C. full two-ring model with the real heater-aligned mismatch Delta_mu for
     auxiliary radius errors of 0, 5, and 20 nm (Reviewer 3.1);
  D. kappa_aux doubled at fixed kappa_P: convergence toward the adiabatic
     limit.

Outputs q_tworing.npz.
"""
import numpy as np
from quantum import CombQuantum
import exp_quantum as Q
from exp_quantum import build, max_squeezing, ODD

dat = np.load("lle_family.npz")
aux = np.load("aux_ring.npz")
kappa = float(dat["kappa"])
M = Q.M
N_MAIN = 2 * M + 1
MU = np.arange(-M, M + 1)
ODD_MASK = MU % 2 != 0
N_AUX = int(ODD_MASK.sum())
AUX_OF_MAIN = {i: N_MAIN + j for j, i in enumerate(np.where(ODD_MASK)[0])}


def build_tworing(psi_fft, zeta0, J, kaux, dm_kappa=None, v=0.0,
                  kaux_i=0.0):
    """Full two-ring model. J, kaux, kaux_i in kappa units; dm_kappa:
    per-odd-mode resonance mismatch (kappa units) or None for perfect
    alignment. Detection at the auxiliary drop port (rate kaux)."""
    base = CombQuantum(Q.centered_psi(psi_fft),
                      Q.zeta_vec(zeta0, v=v),
                      np.full(N_MAIN, 0.5), np.full(N_MAIN, 0.5), M)
    A_main = base.Mdrift[:N_MAIN, :N_MAIN].copy()
    B_main = base.Mdrift[:N_MAIN, N_MAIN:].copy()
    # main-ring rates: even modes keep bus critical coupling (unmonitored
    # here); odd modes lose the whole bare linewidth to unmonitored baths
    ki = np.full(N_MAIN, 0.5)
    kc = np.full(N_MAIN, 0.5)
    ki[ODD_MASK] = 1.0
    kc[ODD_MASK] = 0.0
    nt = N_MAIN + N_AUX
    ki_t = np.concatenate([ki, np.full(N_AUX, kaux_i)])
    kc_t = np.concatenate([kc, np.full(N_AUX, kaux)])
    # rebuild diagonal decay of main block for the new rates
    zv = Q.zeta_vec(zeta0, v=v)
    Delta_main = 0.5 * np.asarray(zv)
    for i in range(N_MAIN):
        A_main[i, i] = (A_main[i, i]
                        + (0.5 + 0.5) / 2 - (ki[i] + kc[i]) / 2)
    A = np.zeros((nt, nt), complex)
    B = np.zeros((nt, nt), complex)
    A[:N_MAIN, :N_MAIN] = A_main
    B[:N_MAIN, :N_MAIN] = B_main
    # auxiliary modes: decay and detuning
    # mismatch dm = omega_main - omega_aux, so the auxiliary resonance sits
    # at Delta_main - dm in the common rotating frame
    dm = np.zeros(N_AUX) if dm_kappa is None else np.asarray(dm_kappa)
    Delta_aux = Delta_main[ODD_MASK] - dm
    for j in range(N_AUX):
        ja = N_MAIN + j
        A[ja, ja] = -(kaux_i + kaux) / 2 - 1j * Delta_aux[j]
    # beam-splitter coupling H = J (b^dag a + a^dag b) -> db/dt += -iJ a
    for i, ja in AUX_OF_MAIN.items():
        A[i, ja] += -1j * J
        A[ja, i] += -1j * J
    cq = CombQuantum.__new__(CombQuantum)
    cq.M = M
    cq.n = nt
    cq.kappa_i = ki_t
    cq.kappa_c = kc_t
    cq.Mdrift = np.block([[A, B], [np.conj(B), np.conj(A)]])
    return cq


DETECT_AUX = list(range(N_MAIN, N_MAIN + N_AUX))


def peak(cq, detect, om_max=60.0, n_om=241):
    oms = np.linspace(0, om_max, n_om)
    smin, smax = cq.squeezing_spectrum(oms, detect)
    i = int(np.argmin(smin))
    return oms, smin, float(oms[i]), float(smin[i])


def bandwidth(oms, smin):
    """Same criterion as exp_quantum.py: full width over which the linear
    noise reduction 1 - S exceeds half its peak value."""
    s_at = smin.min()
    thr = 1 - (1 - s_at) / 2
    sel = oms[smin < thr]
    return float(sel.max() - sel.min()) if len(sel) > 1 else 0.0


if __name__ == "__main__":
    from quantum import log_negativity
    zeta0s = dat["zeta0"]
    iz = int(np.argmin(np.abs(zeta0s - 6.5)))
    z, psi, v = float(zeta0s[iz]), dat["psi"][iz], float(dat["v"][iz])
    KAUX2GHZ = float(aux["kaux"]) / kappa      # 15.51 kappa  (2 GHz)
    out = {}
    print(f"zeta0 = {z}, kappa_aux(2 GHz) = {KAUX2GHZ:.2f} kappa")

    # sanity: vacuum output with no comb
    cq0 = build_tworing(np.zeros_like(psi), z,
                        J=np.sqrt(10 * KAUX2GHZ) / 2, kaux=KAUX2GHZ, v=0.0)
    V0 = cq0.quad_covariance_sym(1.3, DETECT_AUX)
    print("vacuum check (should be ~0):",
          np.max(np.abs(V0 - np.eye(2 * N_AUX))))

    # ---- A: fixed-pump optimum scan with the full model
    kps = np.array([5.0, 10.0, 15.0, 20.0, 30.0])
    kauxs = np.array([KAUX2GHZ, 2 * KAUX2GHZ, 4 * KAUX2GHZ])  # 2/4/8 GHz
    grid_db = np.zeros((len(kps), len(kauxs)))
    for i, kp in enumerate(kps):
        # adiabatic reference
        cqa = build(psi, z, kc_odd=kp, v=v)
        oms, smin_a, om_a, s_a = peak(cqa, ODD)
        out[f"smin_ad_{kp:.0f}"] = smin_a
        out[f"peak_ad_{kp:.0f}"] = 10 * np.log10(s_a)
        for j, kaux in enumerate(kauxs):
            cqf = build_tworing(psi, z, J=np.sqrt(kp * kaux) / 2,
                                kaux=kaux, v=v)
            _, smin_f, om_f, s_f = peak(cqf, DETECT_AUX)
            grid_db[i, j] = 10 * np.log10(s_f)
            if j == len(kauxs) - 1:
                out[f"smin_full8_{kp:.0f}"] = smin_f
        print(f"kp={kp:4.0f}k  adiab {10*np.log10(s_a):+.2f}  "
              + "  ".join(f"{grid_db[i, j]:+.2f}" for j in range(len(kauxs))),
              flush=True)
    out["oms"] = oms
    out["kps"] = kps
    out["kauxs"] = kauxs
    out["grid_db"] = grid_db

    # ---- B: adiabatic-elimination convergence at kp = 10 (Reviewer 4.1)
    kp = 10.0
    conv_kaux = np.array([KAUX2GHZ, 2 * KAUX2GHZ, 4 * KAUX2GHZ,
                          8 * KAUX2GHZ, 32 * KAUX2GHZ])
    conv_db = []
    s_ad = 10 ** (out["peak_ad_10"] / 10)
    for kaux in conv_kaux:
        cqf = build_tworing(psi, z, J=np.sqrt(kp * kaux) / 2,
                            kaux=kaux, v=v)
        _, _, _, s_f = peak(cqf, DETECT_AUX)
        conv_db.append(10 * np.log10(s_f))
        print(f"conv kaux={kaux:6.1f}k: {conv_db[-1]:+.3f} dB "
              f"(err {conv_db[-1]-out['peak_ad_10']:+.3f})", flush=True)
    out["conv_kaux"] = conv_kaux
    out["conv_db"] = conv_db

    # ---- C: design point kp=10k, kaux = 8 GHz: spectrum, bandwidth,
    #         supermodes, entanglement, radius-error tolerance
    kaux_d = float(kauxs[-1])
    J_d = np.sqrt(kp * kaux_d) / 2
    cqd = build_tworing(psi, z, J=J_d, kaux=kaux_d, v=v)
    oms_d, smin_d, om_d, s_d = peak(cqd, DETECT_AUX)
    bw_d = bandwidth(oms_d, smin_d)
    print(f"design point: {10*np.log10(s_d):+.3f} dB, "
          f"BW = {bw_d:.2f} kappa = {bw_d*kappa/2/np.pi/1e9:.2f} GHz")
    wD, UD = cqd.supermodes(om_d, DETECT_AUX, k=4)
    Vd = cqd.quad_covariance_sym(om_d, DETECT_AUX)
    dEN = len(DETECT_AUX)
    ENd = np.zeros((dEN, dEN))
    for a in range(dEN):
        for b in range(a + 1, dEN):
            ENd[a, b] = ENd[b, a] = log_negativity(Vd, a, b, dEN)
    print(f"design-point max EN = {ENd.max():.3f}, "
          f"SM1 = {10*np.log10(wD[0]):+.2f} dB, "
          f"SM2 = {10*np.log10(wD[1]):+.2f} dB")
    out["design_smin"] = smin_d
    out["design_peak"] = 10 * np.log10(s_d)
    out["design_bw_kappa"] = bw_d
    out["design_EN"] = ENd
    out["design_sup_w"] = wD
    out["design_sup_U"] = UD
    for name in ("0nm", "5nm", "20nm"):
        dm = aux[f"dm_{name}"] / kappa
        cqr = build_tworing(psi, z, J=J_d, kaux=kaux_d, dm_kappa=dm, v=v)
        _, smin_r, om_r, s_r = peak(cqr, DETECT_AUX)
        out[f"design_real_{name}"] = 10 * np.log10(s_r)
        out[f"design_smin_{name}"] = smin_r
        print(f"design real dR={name}: {10*np.log10(s_r):+.3f} dB")
    np.savez("q_tworing.npz", kappa=kappa, **out)
    print("DONE")
