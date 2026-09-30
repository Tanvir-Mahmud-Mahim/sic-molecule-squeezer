"""Mechanism of the D3 stability boundary: breathing state past the boundary.

Converges the 2-FSR crystal at the last stable D3 scale (3.25, zeta0 = 6.5,
same protocol as exp_d3boundary.py: fresh two-soliton seed, T = 500), then
switches D3 to the first unstable scale (3.5) in the co-moving frame of the
seed and integrates for 150 time units (2/kappa), sampling every 0.25:
the peak field max|psi|, and the comb lines mu = 2 and mu = 88 (near the
dispersive-wave phase-matching point). Also prints the mode where the
integrated dispersion Dint (including D4) changes sign for both scales.
Mean, oscillation amplitude (half the peak-to-peak range) and the dominant
breathing frequency are taken from the second half of the record.

Needs fem_final.json (through exp_lle.py). Outputs d3_mech.npz (used by
make_numbers.py).
"""
import numpy as np
from exp_lle import converge_state, zeta_of, N, F_PUMP, kappa
from lle import LLE

S_LAST, S_LOST = 3.25, 3.5
T_RUN, T_SAMPLE, DT = 150.0, 0.25, 0.002


def mu_dint_zero(s3):
    """Last positive mode number with Dint(mu) > 0."""
    pos = np.arange(1, N)
    return int(pos[zeta_of(pos, s3) > 0].max())


if __name__ == "__main__":
    psi0, r, v0, _ = converge_state(6.5, d3scale=S_LAST, T=500)
    print(f"seed s3={S_LAST}: resid={r:.2e}", flush=True)
    for s3 in (S_LAST, S_LOST):
        print(f"s3={s3}: Dint zero crossing at mu = {mu_dint_zero(s3)}")

    mu = np.fft.fftfreq(N, d=1.0 / N).astype(int)
    sim = LLE(N, 6.5 + zeta_of(mu, S_LOST) - v0 * mu, F_PUMP)
    psi = psi0.copy()
    n_samp = int(round(T_RUN / T_SAMPLE))
    t, peak, line2, line88 = [], [], [], []
    for k in range(n_samp):
        for _ in range(int(T_SAMPLE / DT)):
            psi = sim.step(psi, DT)
        t.append((k + 1) * T_SAMPLE)
        peak.append(np.max(np.abs(np.fft.ifft(psi)) * N ** 0.5))
        line2.append(np.abs(psi[2]))
        line88.append(np.abs(psi[88]))
    t, peak = np.array(t), np.array(peak)
    line2, line88 = np.array(line2), np.array(line88)

    h = n_samp // 2                      # second half: after the transient

    def stats(x):
        x = x[h:]
        return x.mean(), (x.max() - x.min()) / 2

    m, a = stats(peak)
    print(f"peak field: mean {m:.3f}, osc amplitude {a:.4f}")
    m, a = stats(line2)
    print(f"comb line mu=2: mean {m:.4f}, osc {a:.5f}")
    m, a = stats(line88)
    print(f"line mu=88 (near DW): mean {m:.2e}, osc {a:.2e}")
    x = peak[h:] - peak[h:].mean()
    spec = np.abs(np.fft.rfft(x))
    freqs = np.fft.rfftfreq(len(x), d=T_SAMPLE)   # cycles per unit time
    fdom = float(freqs[1 + np.argmax(spec[1:])])
    # MHz conversion as in make_numbers.py (nBreathMHz) and the archived
    # log; see README, "Notes on the calculations", about its units.
    print(f"dominant breathing frequency: {fdom:.3f} per (2/kappa) = "
          f"{fdom*kappa/2/2/np.pi/1e6:.1f} MHz")
    np.savez("d3_mech.npz", t=t, peak=peak, line2=line2, line88=line88,
             fdom=fdom)
