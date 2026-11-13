"""CU2 test 2026-11-07/08 (processed HDF5 from the lab) against the 217 V/cm CU2 model.

    python scripts/plot_cu2_test_vs_sim.py <CU2_test_processed_results_2026-11-07_08.h5>
    -> docs/state_cleanup/figures/CU2_test_2026-11-08_vs_sim.{png,pdf}

(a) CU2 laser scan with the triplet start (RC P(2) F'=2), read on R(2) F1'=7/2 F'=4 (sees F1=5/2 F=3 only),
    against the corrdet stage of scripts/cu2_power_check.py (triplet start, F=3 left).
(b) the same with RC P(2) F'=1, read on R(2) F1'=7/2 F'=3 (97 % F1=5/2 F=2 per molecule at E = 0), against a singlet
    (F1=5/2 F=2 mF=0) start, F=2 left.
(c) S(on)/S(off) of the four J=2 lines at the final CU2 point (NLTL2 2068 MHz) against the transfer stage
    (every J=2 level alone at (13, 209) MHz, 20 + 30 mW, proxy points).
x axis: CU2 laser, MHz IR from 275858.76318 GHz (Cs-corrected, as in the HDF5 B1 scan) = (delta0 - 14.5 MHz) / 4.
The 11-08 scan is placed on it through the 11-07 NLTL2 -> frequency relation (linear fit of B1), as the HDF5 does.
"""
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import vz_paths  # noqa: E402
OUT = ROOT / "tmp" / "cu2_power"
FIG = ROOT / "docs" / "state_cleanup" / "figures"
import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK, MUTED, GRID = "#262624", "#6b6b68", "#e4e3dc"
C_AVG, C_AX, C_DATA = "#2a78d6", "#9a9a94", "#262624"
mpl.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED,
                     "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK, "axes.labelcolor": INK,
                     "text.color": INK})
SIM_TRIP_AX = vz_paths.find(OUT / "corrdet_loss5_f209_20_30_e217_trip_vx0_d-130_170_4_result.npz")
SIM_TRIP_AVG = vz_paths.find(OUT / "corrdet_loss5_f209_20_30_e217_trip_proxy_d-62_90_4_result.npz")
SIM_SING_AX = vz_paths.find(OUT / "corrdet_loss5_f209_20_30_e217_sing_vx0_d-94_98_4_result.npz")
SIM_TRANSFER = vz_paths.find(OUT / "transfer_loss5_d13_f209_20_30_e217_result.npz")
SIM_TRANSFER_ACC = vz_paths.find(OUT / "transfer_loss5_d13_f209_20_30_e217_acc_result.npz")      # 24-point acceptance grid
D0_FINAL = 13.0                                   # MHz UV, the transfer stage set point (NLTL2 2068 -> -0.4 MHz IR)
# E = 0 detection strengths per molecule (sum over mF, mF', x y z, / (2F+1)), from cu2_test_line_frequencies.strength
LINES = {  # line -> {(F1, F): strength per molecule}
    "R(2) F1'=7/2 F'=4": {(2.5, 3): 0.7347},
    "R(2) F1'=7/2 F'=3": {(2.5, 2): 0.7619, (2.5, 3): 0.0272},
    "R(2) F1'=5/2 F'=3": {(1.5, 2): 0.7467, (2.5, 2): 0.0025, (2.5, 3): 0.0363},
    "R(2) F1'=5/2 F'=2": {(1.5, 1): 0.8000, (1.5, 2): 0.0533, (2.5, 2): 0.0356, (2.5, 3): 0.0018},
}
ONOFF = [("R2F4_triplet_nltl2_2068", "R(2) F1'=7/2 F'=4", "F1=5/2 F=3\n(triplet)"),
         ("R2F3_singlet_nltl2_2068", "R(2) F1'=7/2 F'=3", "F1=5/2 F=2\n(singlet)"),
         ("R2F5o2F3_F1eq3o2F2_nltl2_2068", "R(2) F1'=5/2 F'=3", "F1=3/2 F=2"),
         ("R2F5o2F2_F1eq3o2F1_nltl2_2068", "R(2) F1'=5/2 F'=2", "F1=3/2 F=1")]


def wvx(vx):
    w = np.exp(-0.5 * (np.asarray(vx) / 0.65) ** 2)
    return w / w.sum()


def sim_curve(f, key="F3", level=None):
    z = np.load(f, allow_pickle=True)
    x = (z["delta0_MHz"] - 14.5) / 4
    if level is None:
        y = wvx(z["vx"]) @ z[key]
    else:
        sel = [k for k, (a, b) in enumerate(zip(z["F1"], z["F"])) if (a, b) == level]
        y = wvx(z["vx"]) @ z["levels"][:, :, sel].sum(-1)
    return x, y


def start_weights(S_off):
    """Initial population of the four levels from the CU2-off signals of the four lines (E = 0 strengths, linear in
    population; DET saturation ignored) -> {(F1, F): N}."""
    lv = [(2.5, 3), (2.5, 2), (1.5, 2), (1.5, 1)]
    A = np.array([[LINES[line].get(l, 0.0) for l in lv] for _, line, _ in ONOFF])
    N = np.linalg.solve(A, np.array([S_off[k] for k, _, _ in ONOFF]))
    return dict(zip(lv, N / N.sum()))


def acceptance_weights(vx, y0, svx=0.65, sy=10.35e-3, vcut=1.0, ycut=8.5e-3):
    """W[vx, y0] reproducing plot_cleanup_maps_current.acceptance_mean (Gaussian vx, |vx| <= vcut, trapezoid ends;
    y0 Gaussian over 0..ycut with sigma_y at CU2), with linear rather than log interpolation in y0 so that it applies
    to every final level, zeros included."""
    vx, y0 = np.asarray(vx, float), np.asarray(y0, float)
    wv = np.exp(-0.5 * (vx / svx) ** 2) * (np.abs(vx) <= vcut + 1e-12)
    t = np.ones_like(vx); t[[0, -1]] = 0.5; wv = wv * t; wv /= wv.sum()
    ym = np.linspace(0.0, ycut, 171)
    wy = np.exp(-0.5 * (ym / sy) ** 2); t = np.ones_like(ym); t[[0, -1]] = 0.5; wy = wy * t; wy /= wy.sum()
    M = np.array([np.interp(ym, y0, np.eye(len(y0))[k]) for k in range(len(y0))])    # y0 node -> ym
    return np.outer(wv, M @ wy)


def transfer_prediction(weights, point=None, sim=None, W=None):
    """S(on)/S(off) of each line from the transfer stage: each level starts uniform over its mF, except the lens
    levels F1=5/2 F=3 (mF = 0, +-1) and F=2 (mF = 0). point = (vx index, y0 index) for one grid point only;
    W[vx, y0] = trajectory weights (default: Gaussian vx, y0 points equal), or {start F1: W}."""
    z = np.load(sim or SIM_TRANSFER, allow_pickle=True)
    if point is None:
        if W is None:
            W = np.outer(wvx(z["vx"]), np.full(len(z["y0"]), 1 / len(z["y0"])))
        if isinstance(W, dict):                    # {start F1: W[vx, y0]}: separate distributions per F1
            fin = np.array([np.einsum("vy,vyl->l", W[float(f1)], z["final"][s])
                            for s, f1 in enumerate(z["start_F1"])])
        else:
            fin = np.einsum("vy,svyl->sl", W, z["final"])                          # start (mF >= 0) x final level
    else:
        fin = z["final"][:, point[0], point[1]]
    lab = list(zip(z["F1"], z["F"], z["mF"]))
    st = list(zip(z["start_F1"], z["start_F"], z["start_mF"]))
    mirror = {l: lab.index((l[0], l[1], -l[2])) for l in lab}

    def evolve(start_lab):                       # final level populations from one start level (mF < 0 mirrored)
        if start_lab[2] >= 0:
            return fin[st.index(start_lab)]
        f = fin[st.index((start_lab[0], start_lab[1], -start_lab[2]))]
        return np.array([f[mirror[l]] for l in lab])
    rho0 = np.zeros(len(lab))
    for (f1, F), w in weights.items():
        ms = [0, -1, 1] if (f1, F) == (2.5, 3) else ([0] if (f1, F) == (2.5, 2) else list(range(-int(F), int(F) + 1)))
        for m in ms:
            rho0[lab.index((f1, F, m))] += w / len(ms)
    after = sum(r * evolve(l) for l, r in zip(lab, rho0) if r > 0)
    out = {}
    for k, line, _ in ONOFF:
        s = np.array([LINES[line].get((l[0], l[1]), 0.0) for l in lab])
        out[k] = float(s @ after / (s @ rho0))
    own = {}
    for (f1, F) in weights:                       # each level alone, fraction of it left in itself
        sel = [i for i, l in enumerate(lab) if (l[0], l[1]) == (f1, F)]
        r = np.zeros(len(lab)); r[[i for i in sel if rho0[i] > 0]] = rho0[[i for i in sel if rho0[i] > 0]]
        a = sum(q * evolve(l) for l, q in zip(lab, r) if q > 0)
        own[(f1, F)] = (float(a[sel].sum() / r.sum()), float(a.sum() / r.sum()))
    return out, own


def main(h5path):
    f = h5py.File(h5path, "r")
    b1 = f["cu2_frequency_scans/B1_pm50MHz_on_R2F4"]
    nl1, x1, R1, e1 = (b1[k][()] for k in ("NLTL2_MHz", "dIR_MHz", "R_template", "sem_bootstrap_R_template"))
    slope, icpt = np.polyfit(nl1, x1, 1)
    to_x = lambda nl: slope * np.asarray(nl) + icpt                                      # noqa: E731
    b4 = f["cu2_frequency_scans/B4a_pm24MHz_on_R2F3"]
    x4, R4, e4 = to_x(b4["NLTL2_MHz"][()]), b4["R_template"][()], b4["sd_bootstrap"][()]
    print(f"NLTL2 -> x: slope {slope:.4f} MHz IR / MHz, 2071.5 -> {to_x(2071.5):+.2f}, 2068 -> {to_x(2068):+.2f}")

    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1.35, 1], hspace=0.42, wspace=0.22)
    axa, axb, axc = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    for ax in (axa, axb):
        ax.axhline(0, color=GRID, lw=0.8, zorder=0); ax.axhline(1, color=GRID, lw=0.8, zorder=0)
        for nl, lab in ((2071.5, "11-07 lock"), (2068, "11-08 final")):
            ax.axvline(to_x(nl), color=MUTED, lw=0.7, ls=(0, (2, 2)), zorder=1)
        ax.set_ylabel("S(CU2 on) / S(CU2 off)")
        ax.set_ylim(-0.08, 1.5)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.spines[["top", "right"]].set_visible(False)
    # (a) triplet
    axa.errorbar(x1, R1, e1, fmt="o", ms=3.5, color=C_DATA, ecolor=MUTED, elinewidth=0.8, capsize=0, zorder=4,
                 label="measured 11-07, R(2) F1'=7/2 F'=4 (F1=5/2 F=3), RC P(2) F'=2")
    if SIM_TRIP_AX.exists():
        xs, ys = sim_curve(SIM_TRIP_AX)
        axa.plot(xs, ys, color=C_AX, lw=1.0, zorder=2, label="model, on axis (vx = 0, y0 = 0)")
    if SIM_TRIP_AVG.exists():
        xs, ys = sim_curve(SIM_TRIP_AVG)
        axa.plot(xs, ys, color=C_AVG, lw=2.0, zorder=3,
                 label="model, averaged over vx = 0, ±0.6 m/s (Gaussian, σ 0.65) and y0 = 0, 6 mm")
    axa.set_xlim(-58, 48)
    axa.set_title("(a) triplet start: F1=5/2 F=3 left vs CU2 laser, 20 + 30 mW, f_AOM 209 MHz, 217 V/cm",
                  loc="left", fontsize=9.5)
    axa.legend(frameon=False, loc="upper left", fontsize=8, ncol=1)
    for ax in (axa, axb):
        ax.text(to_x(2071.5) - 0.3, 0.5, "11-07 lock", rotation=90, ha="right", va="center", fontsize=7, color=MUTED)
        ax.text(to_x(2068) + 0.3, 0.5, "11-08 final", rotation=90, ha="left", va="center", fontsize=7, color=MUTED)
    for xt in ((-49, 0.12), (36, -0.02)):
        axa.text(*xt, "other beam's line\n(not in model)", fontsize=7, color=MUTED, ha="center", va="center")
    # (b) singlet
    axb.errorbar(x4, R4, np.minimum(e4, 0.3), fmt="o", ms=3.5, color=C_DATA, ecolor=MUTED, elinewidth=0.8, capsize=0,
                 zorder=4, label="measured 11-08, R(2) F1'=7/2 F'=3 (97 % F1=5/2 F=2), RC P(2) F'=1")
    if SIM_SING_AX.exists():
        xs, ys = sim_curve(SIM_SING_AX, level=(2.5, 2.0))
        xs3, ys3 = sim_curve(SIM_SING_AX, level=(2.5, 3.0))
        yl = (0.7619 * ys + 0.0272 * ys3) / 0.7619
        axb.plot(xs, yl, color=C_AVG, lw=1.6, zorder=3, label="model, singlet start (F=2 mF=0), on axis")
    if SIM_TRIP_AX.exists():
        xs, ys = sim_curve(SIM_TRIP_AX)
        axb.plot(xs, ys, color=C_AX, lw=0.9, ls="--", zorder=2, label="model, triplet start (panel a), on axis")
    axb.set_xlim(-30, 24)
    axb.set_xlabel("CU2 laser − 275858.76318 GHz (MHz IR, Cs-corrected)  =  (δ0 − 14.5 MHz) / 4")
    axb.set_title("(b) singlet start: F1=5/2 F=2 left", loc="left", fontsize=9.5)
    axb.legend(frameon=False, loc="upper left", fontsize=7.5)
    axa.set_xlabel("CU2 laser − 275858.76318 GHz (MHz IR, Cs-corrected)  =  (δ0 − 14.5 MHz) / 4")
    # (c) on/off at the final point
    g = f["cu2_onoff"]
    S_off = {k: float(g[k].attrs["fl_closed_per_shot_Vms"]) for k, _, _ in ONOFF}
    meas = [(float(g[k].attrs["R_template"]), float(g[k].attrs["R_template_err"])) for k, _, _ in ONOFF]
    xx = np.arange(len(ONOFF))
    axc.errorbar(xx - 0.12, [m for m, _ in meas], [e for _, e in meas], fmt="o", ms=5, color=C_DATA, ecolor=MUTED,
                 elinewidth=1.0, capsize=0, label="measured 11-08 (NLTL2 2068)", zorder=4)
    if SIM_TRANSFER_ACC.exists():
        w = start_weights(S_off)
        z = np.load(SIM_TRANSFER_ACC, allow_pickle=True)
        nv, ny = len(z["vx"]), len(z["y0"])
        pred = {s_: transfer_prediction(w, sim=SIM_TRANSFER_ACC, W=acceptance_weights(z["vx"], z["y0"], svx=s_))[0]
                for s_ in (0.65, 1.28, 0.3)}
        pts = [transfer_prediction(w, (i, j), sim=SIM_TRANSFER_ACC)[0] for i in range(nv) for j in range(ny)]
        lo = [min(q[k] for q in pts) for k, _, _ in ONOFF]; hi = [max(q[k] for q in pts) for k, _, _ in ONOFF]
        axc.vlines(xx + 0.12, lo, hi, color=C_AVG, lw=1.0, alpha=0.5, zorder=3)
        axc.plot(xx + 0.12, [pred[0.65][k] for k, _, _ in ONOFF], "D", ms=5, color=C_AVG, zorder=5,
                 label=f"model, δ0 = {D0_FINAL:g} MHz, 20 + 30 mW, acceptance mean, σvx 0.65 m/s\n"
                       "(line: range over the 24 grid points, |vx| ≤ 1 m/s, y0 0–9 mm)")
        axc.plot(xx + 0.24, [pred[1.28][k] for k, _, _ in ONOFF], "s", ms=4, mfc="white", mec=C_AVG, zorder=5,
                 label="model, σvx 1.28 m/s")
        axc.plot(xx + 0.24, [pred[0.3][k] for k, _, _ in ONOFF], "v", ms=4.5, mfc="white", mec="#eb6834", zorder=5,
                 label="model, σvx 0.3 m/s")
        print("start weights", {k: round(v, 3) for k, v in w.items()})
        for (k, line, lab), (m, e), l, h in zip(ONOFF, meas, lo, hi):
            print(f"{lab.splitlines()[0]:12s} measured {m:.4f} +- {e:.4f}  model 0.65 {pred[0.65][k]:.4f}  "
                  f"1.28 {pred[1.28][k]:.4f}  0.3 {pred[0.3][k]:.4f}  grid {l:.4f}-{h:.4f}  "
                  f"({(m - pred[0.65][k]) / e:+.1f} sigma)")
    axc.axhline(0, color=GRID, lw=0.8, zorder=0)
    axc.set_xticks(xx, [lab for _, _, lab in ONOFF], fontsize=8)
    axc.set_ylabel("S(CU2 on) / S(CU2 off)")
    axc.set_yscale("symlog", linthresh=0.01, linscale=1.5)
    axc.set_ylim(-0.012, 12.0)
    axc.set_yticks([-0.01, 0, 0.005, 0.01, 0.03, 0.1, 0.3], ["−0.01", "0", "0.005", "0.01", "0.03", "0.1", "0.3"])
    axc.minorticks_off()
    axc.spines[["top", "right"]].set_visible(False)
    axc.grid(axis="y", color=GRID, lw=0.5)
    axc.set_title("(c) four J=2 levels at the final point, E = 0 detection", loc="left", fontsize=9.5)
    axc.legend(frameon=False, loc="upper left", fontsize=7.5)
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"CU2_test_2026-11-08_vs_sim.{ext}", dpi=200, bbox_inches="tight")
    print("wrote", FIG / "CU2_test_2026-11-08_vs_sim.png")


if __name__ == "__main__":
    main(sys.argv[1])
