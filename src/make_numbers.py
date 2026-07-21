"""Generate numbers.tex: every quantitative claim in the manuscript is a
macro extracted directly from the simulation outputs."""
import numpy as np, json

fem = json.load(open("fem_final.json"))
from dispersion import RingDispersion
d = fem["final"]
rd = RingDispersion(d["lams"], d["neff"])
R = rd.radius_for_fsr(350e9)
res = rd.resonances(R, mu_max=80)
om0 = 2 * np.pi * 299792458 / 1.55e-6

dat = np.load("lle_family.npz")
kappa = float(dat["kappa"])
g0 = float(dat["g0"])
P1 = float(dat["P1mW"])
f = float(dat["f"])
zeta0 = dat["zeta0"]; resid = dat["resid"]; crystal = dat["crystal"]
stat = (resid < 1e-2) & crystal
vs = dat["v"]
vhz = vs * (kappa / 2) / 2 / np.pi / 1e6

rep = np.load("q_rep.npz")
add = np.load("q_addendum.npz")
sw = np.load("q_sweep.npz", allow_pickle=True)["rows"]
ent = np.load("q_entanglement.npz")
d3f = np.load("lle_d3family.npz")
d3q = np.load("q_d3.npz")

conv = fem["convergence"]
neff_best = [c["neff"] for c in conv if c["order"] == 2 and c["res"] == 0.03][0]
neff_04 = [c["neff"] for c in conv if c["order"] == 2 and c["res"] == 0.04][0]

kps = rep["kps"]; mmin = np.array(rep["mol_min_db"])
iopt = int(np.argmin(mmin))
bw = np.array(rep["mol_bw"]) * (kappa / 2 / np.pi) / 1e9  # GHz
i_bare = int(np.where(kps == 0.5)[0][0])

M = {}
M["nWidth"] = "1.85"
M["nHeight"] = "500"
M["nRadius"] = f"{R:.2f}"
M["nFSR"] = "350"
M["nNeff"] = f"{np.interp(1.55, d['lams'], d['neff']):.3f}"
M["nNg"] = f"{rd.ng(om0):.2f}"
M["nAeff"] = f"{fem['aeff_um2']:.2f}"
M["nDtwo"] = f"{res['D2']/2/np.pi/1e6:.2f}"
M["nDthree"] = f"{res['D3']/2/np.pi/1e3:.1f}"
M["nDfour"] = f"{res['D4']/2/np.pi:.0f}"
M["nKappa"] = f"{kappa/2/np.pi/1e6:.0f}"
M["nGzero"] = f"{g0/2/np.pi:.1f}"
M["nPone"] = f"{P1:.2f}"
M["nPump"] = f"{P1*f**2:.1f}"
M["nFsq"] = f"{f**2:.0f}"
M["nNeffConv"] = f"{abs(neff_best-neff_04)*1e7:.0f}"
M["nZetaLow"] = f"{zeta0[stat].min():.2f}"
M["nZetaHigh"] = f"{zeta0[stat].max():.2f}"
M["nZetaMax"] = f"{np.pi**2*f**2/8:.1f}"
M["nRepShiftLow"] = f"{abs(vhz[stat]).min():.1f}"
M["nRepShiftHigh"] = f"{abs(vhz[stat]).max():.1f}"
M["nBareSq"] = f"{-min(r['out_min_db'] for r in sw):.1f}"
M["nBareSqRange"] = (f"{-max(r['out_min_db'] for r in sw):.1f}--"
                     f"{-min(r['out_min_db'] for r in sw):.1f}")
M["nIdealRep"] = f"{-float(add['ideal_db'][np.argmin(np.abs(add['zeta0']-6.5))]):.1f}"
zz = np.array(add['zeta0']); idl = np.array(add['ideal_db'])
M["nIdealMax"] = f"{-idl[zz <= 10.0].min():.1f}"
M["nIdealTwentyZeta"] = f"{zz[idl < -20][0]:.2f}"
M["nMolKp"] = f"{kps[iopt]:.0f}"
M["nMolSq"] = f"{-mmin[iopt]:.1f}"
M["nMolEta"] = f"{rep['mol_eta'][iopt]:.2f}"
M["nMolSqTwenty"] = f"{-mmin[int(np.where(kps==20)[0][0])]:.1f}"
M["nMolSqTwentyMax"] = f"{-np.min(add['mol20_db']):.1f}"
M["nMolTenMax"] = f"{-np.min(add['mol10_db']):.1f}"
M["nBWbare"] = f"{bw[i_bare]:.2f}"
M["nBWopt"] = f"{bw[iopt]:.2f}"
M["nBWtwenty"] = f"{bw[int(np.where(kps==20)[0][0])]:.2f}"
M["nBWxOpt"] = f"{bw[iopt]/bw[i_bare]:.1f}"
M["nBWxTwenty"] = f"{bw[int(np.where(kps==20)[0][0])]/bw[i_bare]:.1f}"
M["nSMoneBare"] = f"{-10*np.log10(rep['sup_bare_w'][0]):.1f}"
M["nSMtwoBare"] = f"{-10*np.log10(rep['sup_bare_w'][1]):.1f}"
M["nSMthreeBare"] = f"{-10*np.log10(rep['sup_bare_w'][2]):.2f}"
M["nSMoneMol"] = f"{-10*np.log10(rep['sup_mol_w'][0]):.1f}"
M["nSMtwoMol"] = f"{-10*np.log10(rep['sup_mol_w'][1]):.1f}"
M["nENmax"] = f"{ent['EN'].max():.2f}"
M["nDthreeBnd"] = f"{float(d3f['boundary_scale']):.1f}"
M["nDthreeBndkHz"] = f"{float(d3f['boundary_scale'])*abs(res['D3'])/2/np.pi/1e3:.0f}"
M["nDthreeSqSpread"] = f"{np.max(np.abs(np.array(d3q['min_db'])-d3q['min_db'][0])):.2f}"
M["nOddSupp"] = f"{abs(np.median(dat['odd_even_db'][stat])):.0f}"
M["nSqPerPump"] = f"{-mmin[iopt]/(P1*f**2):.2f}"

# odd/even suppression, threshold estimate
with open("../paper/numbers.tex", "w") as fh:
    for k, v in M.items():
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
with open("../supplement/numbers.tex", "w") as fh:
    for k, v in M.items():
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print(json.dumps(M, indent=0))
