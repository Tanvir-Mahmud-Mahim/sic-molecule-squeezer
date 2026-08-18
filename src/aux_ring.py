"""Explicit two-ring alignment analysis for the photonic molecule.

Computes the auxiliary-ring resonance grid from the SAME FEM-fitted
propagation constant beta(omega) as the main ring (bending corrections
neglected for both, see Supplement), and from it:

- the resonance mismatch Delta_mu between each odd main-ring mode and its
  nearest auxiliary-ring resonance, after a global thermo-optic alignment
  offset (microheater) chosen to minimize the worst-case mismatch across
  the analyzed lattice |mu| <= 30;
- the resulting mode-dependent Purcell rate
      kappa_P(mu) = kappa_P0 / (1 + (2 Delta_mu / kappa_aux)^2)
  and dispersive shift
      delta_omega(mu) = -J^2 Delta_mu / (Delta_mu^2 + (kappa_aux/2)^2)
  from adiabatic elimination of the detuned auxiliary mode;
- the same quantities for auxiliary-ring radius errors dR (fabrication),
  where the heater still absorbs the global offset but the residual
  FSR error produces a linear walk-off across the lattice.

Outputs aux_ring.npz used by exp_tworing.py and the figures.
"""
import json
import numpy as np
from dispersion import RingDispersion

C = 299792458.0

fem = json.load(open("fem_final.json"))
d = fem["final"]
rd = RingDispersion(d["lams"], d["neff"])

# main ring: exactly the published geometry
R = rd.radius_for_fsr(350e9)              # um
main = rd.resonances(R, mu_max=40)
# auxiliary ring: FSR = 2 x main FSR -> half radius (same n_g at 1.55 um)
Raux = rd.radius_for_fsr(700e9)           # um  (= R/2 to numerical precision)
aux = rd.resonances(Raux, mu_max=25)

dat = np.load("lle_family.npz")
kappa = float(dat["kappa"])               # rad/s

MU_ODD = np.array([m for m in range(-30, 31) if m % 2 != 0])


def mismatch(offset, dfsr=0.0):
    """Delta_mu (rad/s) between odd main modes and nearest aux resonance.
    offset: global aux-grid shift (rad/s), thermo-optic.
    dfsr: aux FSR error (rad/s per aux index), from radius error."""
    out = []
    for mu in MU_ODD:
        om_main = np.interp(mu, main["mu"], main["omega"])
        oms_aux = aux["omega"] + offset + dfsr * aux["mu"]
        out.append(om_main - oms_aux[np.argmin(np.abs(om_main - oms_aux))])
    return np.array(out)


def best_offset(dfsr=0.0):
    """Global offset minimizing the worst-case |Delta_mu| (Chebyshev).
    Coarse search over one auxiliary FSR, then two local refinements."""
    span = float(aux["D1"])
    offs = np.linspace(-span / 2, span / 2, 4001)
    worst = [np.max(np.abs(mismatch(o, dfsr))) for o in offs]
    i = int(np.argmin(worst))
    o_best, w_best = offs[i], worst[i]
    step = offs[1] - offs[0]
    for _ in range(2):
        offs = np.linspace(o_best - 2 * step, o_best + 2 * step, 2001)
        worst = [np.max(np.abs(mismatch(o, dfsr))) for o in offs]
        i = int(np.argmin(worst))
        o_best, w_best = offs[i], worst[i]
        step = offs[1] - offs[0]
    return o_best, w_best


if __name__ == "__main__":
    D1 = float(main["D1"])
    D1aux = float(aux["D1"])
    print(f"R = {R:.3f} um, Raux = {Raux:.4f} um, R/2 = {R/2:.4f} um")
    print(f"main FSR = {D1/2/np.pi/1e9:.4f} GHz, aux FSR = "
          f"{D1aux/2/np.pi/1e9:.4f} GHz, 2*main = {2*D1/2/np.pi/1e9:.4f}")
    print(f"aux D2/2pi = {aux['D2']/2/np.pi/1e6:.2f} MHz per aux index")

    # ---- ideal radius: heater-aligned residual mismatch
    off0, worst0 = best_offset()
    dm0 = mismatch(off0)
    print(f"ideal radius: offset = {off0/2/np.pi/1e9:.4f} GHz, "
          f"max |Delta| = {worst0/2/np.pi/1e6:.2f} MHz")

    # ---- radius errors: dFSR_aux/FSR_aux = -dR/Raux (n_g fixed)
    cases = {"0nm": 0.0}
    for dr_nm in (5.0, 20.0):
        dR_um = dr_nm * 1e-3
        dfsr = -D1aux * (dR_um / Raux)    # rad/s per aux index
        off, worst = best_offset(dfsr)
        cases[f"{dr_nm:.0f}nm"] = dfsr
        print(f"dR = {dr_nm:.0f} nm: dFSR = {dfsr/2/np.pi/1e6:.1f} MHz, "
              f"max |Delta| = {worst/2/np.pi/1e9:.3f} GHz "
              f"(= {worst/kappa:.2f} kappa_aux units below)")

    # ---- kappa_P(mu) and shifts for the design point
    kaux_hz = 2e9                          # kappa_aux / 2pi (Hz)
    kaux = 2 * np.pi * kaux_hz             # rad/s
    out = dict(mu_odd=MU_ODD, kappa=kappa, kaux=kaux,
               D1=D1, D1aux=D1aux, D2aux=float(aux["D2"]),
               R_um=R, Raux_um=Raux)
    for kp_target in (10.0, 20.0):         # in units of kappa
        kP0 = kp_target * kappa
        J = np.sqrt(kP0 * kaux) / 2        # rad/s
        out[f"J_{kp_target:.0f}"] = J
        for name, dr_nm in (("0nm", 0.0), ("5nm", 5.0), ("20nm", 20.0)):
            dfsr = -D1aux * (dr_nm * 1e-3 / Raux)
            off, _ = best_offset(dfsr)
            dm = mismatch(off, dfsr)
            kP_mu = kP0 / (1 + (2 * dm / kaux) ** 2)
            shift = -J ** 2 * dm / (dm ** 2 + (kaux / 2) ** 2)
            out[f"dm_{name}"] = dm
            out[f"kPmu_{kp_target:.0f}_{name}"] = kP_mu / kappa
            out[f"shift_{kp_target:.0f}_{name}"] = shift / kappa
            print(f"kP={kp_target:.0f}k {name}: kP(mu) range "
                  f"{kP_mu.min()/kappa:.2f}-{kP_mu.max()/kappa:.2f} kappa, "
                  f"max|shift| = {np.max(np.abs(shift))/kappa:.3f} kappa")
    np.savez("aux_ring.npz", **out)
    print("DONE")
