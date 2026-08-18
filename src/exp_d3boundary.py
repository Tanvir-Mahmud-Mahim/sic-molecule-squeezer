"""D3 stability boundary of the 2-FSR crystal, with a documented scan step.

Scans the D3 scaling factor from 3.0 to 4.0 in steps of 0.25 at the
representative detuning zeta0 = 6.5 (same protocol as the D3 family:
fresh two-soliton seed, T = 500 co-moving-frame convergence). For each
scale it records the co-moving residual, drift velocity, odd/even modal
power contrast, and the crystal test result. It also diagnoses the loss
mechanism at the first unstable scale by tracking the growth of odd-mode
power from a converged sub-boundary crystal.

Outputs d3_boundary.npz.
"""
import numpy as np
from exp_lle import converge_state, is_crystal, zeta_of, N, F_PUMP, \
    D2, D3, kappa
from lle import LLE

if __name__ == "__main__":
    scales = np.arange(3.0, 4.01, 0.25)
    rows = []
    psis = {}
    for s3 in scales:
        psi, r, v, _ = converge_state(6.5, d3scale=float(s3), T=500)
        ok, ratio = is_crystal(psi)
        rows.append(dict(s3=float(s3), resid=float(r), v=float(v),
                         odd_even_db=float(ratio),
                         crystal=bool(ok) and r < 1e-2))
        psis[f"psi_{s3:g}"] = psi
        print(f"s3={s3:.2f}: resid={r:.2e} v={v:+.3e} "
              f"odd/even={ratio:.1f} dB crystal={rows[-1]['crystal']}",
              flush=True)
    lost = [r["s3"] for r in rows if not r["crystal"]]
    boundary = min(lost) if lost else None
    print(f"boundary: crystal lost at s3 = {boundary} (step 0.25)")

    # ---- mechanism: seed the *converged* 3.25 crystal at the unstable
    # scale and track odd-mode power growth (C2 symmetry-breaking rate)
    ok_rows = [r for r in rows if r["crystal"]]
    s_last = max(r["s3"] for r in ok_rows)
    psi0 = psis[f"psi_{s_last:g}"]
    v0 = [r["v"] for r in rows if r["s3"] == s_last][0]
    mu = np.fft.fftfreq(N, d=1.0 / N).astype(int)
    s_un = boundary
    zeta = 6.5 + zeta_of(mu, s_un) - v0 * mu
    sim = LLE(N, zeta, F_PUMP)
    psi = psi0.copy()
    odd = mu % 2 != 0
    even = (mu % 2 == 0) & (mu != 0)
    t, todd = [], []
    dt = 0.002
    for k in range(40):
        for _ in range(int(2.5 / dt)):
            psi = sim.step(psi, dt)
        podd = np.sum(np.abs(psi[odd]) ** 2)
        peven = np.sum(np.abs(psi[even]) ** 2)
        t.append((k + 1) * 2.5)
        todd.append(10 * np.log10(podd / peven))
        if k % 8 == 0:
            print(f"t={t[-1]:6.1f}: odd/even = {todd[-1]:.1f} dB", flush=True)
    t = np.array(t)
    todd = np.array(todd)
    # exponential growth rate of the odd-mode field (dB are 10log10 power)
    sel = (todd > -160) & (todd < -40)
    rate = np.polyfit(t[sel], todd[sel], 1)[0] if sel.sum() > 3 else np.nan
    print(f"odd-mode power growth: {rate:.2f} dB per unit time "
          f"(amplitude rate {rate/20*np.log(10):.3f} per 2/kappa)")
    np.savez("d3_boundary.npz",
             s3=[r["s3"] for r in rows],
             resid=[r["resid"] for r in rows],
             v=[r["v"] for r in rows],
             odd_even_db=[r["odd_even_db"] for r in rows],
             crystal=[r["crystal"] for r in rows],
             boundary=boundary, step=0.25,
             growth_t=t, growth_odd_db=todd, growth_rate_db=rate,
             **psis)
    print("DONE")
