"""Slide schematic: CU2 cleanup efficiency seen by each detection line, measured (CU2 test 2026-11-08, NLTL2 2068)
and simulated, with every J=2 level weighted by its probe photon yield.

    python -I scripts/cu2_test_efficiency_schematic.py <CU2_test_processed_results_2026-11-07_08.h5> [power_mW] [notext] [nocl]
    -> docs/state_cleanup/figures/slides/CU2_test_2026-11-08_efficiency_schematic{,_notext}{,_nocl}.{png,pdf}
    (notext: without the explanation block underneath; that text goes on the slide.
     nocl: left panel only - boxes show measured and model without the 90 % CL lower limit, and no 'CU2 off'
     signal shares under the boxes; the right panel keeps its 90 % CL bars)

Measured efficiency: HDF5 attrs clean_eff = 1 - posterior mean of R = S(CU2 on)/S(CU2 off) with the prior R >= 0,
and clean_eff_lower90 (90 % CL lower limit), as in the CU2 on/off traces figure; the model gives 1 - R. Model (scripts/det_plates_level_diagrams.py,
residual_cu2test): DET plates off (E = 0), sigma+ about k at power_mW (default 24 mW, measured 2026-11-09); photons per
molecule for every J=2 level from the probe OBE (Gaussian vx average, sigma 0.65 m/s); start populations fitted
(non-negative least squares) to the four measured CU2-off signals with these yields; CU2 transfer stage (217 V/cm,
(13, 209) MHz, 20 + 30 mW, 24-point acceptance mean). At zero field the 20 J=2 levels lie within 0.36 MHz, so each
line also excites the other levels its excited state couples to; the probe saturates those weak couplings.
The previous comparison (one-photon line strengths, linear) is printed only (scripts/plot_cu2_test_vs_sim.py).
Level energies: zero field, from generate_reduced_X/B_hamiltonian, drawn in order but not to scale.
"""
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import det_plates_level_diagrams as D  # noqa: E402
import plot_cu2_test_vs_sim as P  # noqa: E402

import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

OUT = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
INK, GRAY, LGRAY, GRID = "#111111", "#555555", "#bdbdbd", "#e6e6e6"
COLORS = ["#2a78d6", "#e8673a", "#2e8b57", "#7a3db8"]
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "mathtext.fontset": "dejavusans",
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.7, "axes.edgecolor": GRAY, "axes.labelcolor": INK,
    "xtick.color": GRAY, "ytick.color": GRAY, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": "white",
})
X_LEV = {(2.5, 3): (1.45, "+359 kHz"), (2.5, 2): (1.10, "+323 kHz"), (1.5, 2): (0.35, "+44 kHz"),
         (1.5, 1): (0.0, "0")}
B_LEV = {"7/2,4": (7.35, "+180 MHz"), "7/2,3": (7.0, "0"), "5/2,3": (8.55, "+921 MHz"), "5/2,2": (8.2, "+577 MHz")}
PAIRS = [((2.5, 3), "7/2,4"), ((2.5, 2), "7/2,3"), ((1.5, 2), "5/2,3"), ((1.5, 1), "5/2,2")]
LINE_KEYS = ["F4", "R7o2F3", "F1p5F3", "F1p5F2"]           # det_plates_level_diagrams / cu2_detection_test names
GROUND_NAME = {(2.5, 3): "F₁ = 5/2, F = 3", (2.5, 2): "F₁ = 5/2, F = 2", (1.5, 2): "F₁ = 3/2, F = 2",
               (1.5, 1): "F₁ = 3/2, F = 1"}
SHORT = {(2.5, 3): "5/2,3", (2.5, 2): "5/2,2", (1.5, 2): "3/2,2", (1.5, 1): "3/2,1"}
EXC_NAME = {"7/2,4": "F₁′ = 7/2, F′ = 4", "7/2,3": "F₁′ = 7/2, F′ = 3", "5/2,3": "F₁′ = 5/2, F′ = 3",
            "5/2,2": "F₁′ = 5/2, F′ = 2"}


def numbers(h5path, power_mW):
    with h5py.File(h5path, "r") as f:
        g = f["cu2_onoff"]
        S_off = {k: float(g[k].attrs["fl_closed_per_shot_Vms"]) for k, _, _ in P.ONOFF}
        meas = [(float(g[k].attrs["clean_eff"]), float(g[k].attrs["clean_eff_lower90"]),
                 float(g[k].attrs["R_template"]), float(g[k].attrs["R_template_err"])) for k, _, _ in P.ONOFF]
    _, _, t = D.residual_cu2test(h5path, power_mW)
    lab = t["labels"]
    model = [t["ratio"][k] for k in LINE_KEYS]
    before = [t["shares"][k] for k in LINE_KEYS]                     # CU2-off signal by start group
    after = []                                                        # CU2-on signal by final (F1, F), / CU2-off
    for k in LINE_KEYS:
        s0 = sum(t["Y"][k][i] * t["rho0"][l] for i, l in enumerate(lab))
        after.append({g_: sum(t["Y"][k][i] * t["after"][l] for i, l in enumerate(lab) if (l[0], l[1]) == g_) / s0
                      for g_ in X_LEV})
    old = P.transfer_prediction(P.start_weights(S_off), sim=P.SIM_TRANSFER_ACC,
                                W=P.acceptance_weights(*(np.load(P.SIM_TRANSFER_ACC)[a] for a in ("vx", "y0"))))[0]
    old = [old[k] for k, _, _ in P.ONOFF]
    return meas, model, old, before, after, t


def main(h5path, power_mW=24, notext=False, nocl=False):
    meas, model, old, before, after, t = numbers(h5path, power_mW)
    eff_m = [(100 * ce, 100 * lo) for ce, lo, _, _ in meas]            # posterior mean, 90 % CL lower limit
    eff_t = [100 * (1 - r) for r in model]
    eff_o = [100 * (1 - r) for r in old]
    for k, (m, e), tt, o, b, a in zip(LINE_KEYS, eff_m, eff_t, eff_o, before, after):
        print(f"{k:7s} measured {m:.2f} % (> {e:.2f} %, 90 % CL)  model {tt:.2f} %  old {o:.2f} %  "
              f"CU2-off {[(SHORT[g_], round(v, 3)) for g_, v in b.items() if v > 0.005]}  "
              f"left {[(SHORT[g_], round(v, 4)) for g_, v in a.items() if v > 1e-4]}")

    if notext:                                                         # no explanation block (text goes on the slide)
        fig = plt.figure(figsize=(8.4, 5.0))
        gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1], left=0.02, right=0.985, top=0.86, bottom=0.11,
                              wspace=0.16)
        ax = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1]); axt = None
    else:
        fig = plt.figure(figsize=(8.4, 6.3))
        gs = fig.add_gridspec(2, 2, width_ratios=[1.75, 1], height_ratios=[4.3, 1.25], left=0.02, right=0.985,
                              top=0.875, bottom=0.01, wspace=0.16, hspace=0.22)
        ax = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1]); axt = fig.add_subplot(gs[1, :])
        axt.axis("off")
    ax.axis("off")

    # ---- (a) level scheme: solid arrow = the line's named level, dashed = other levels in its CU2-off signal
    ax.set_xlim(-0.02, 1.47); ax.set_ylim(-0.55, 9.35)
    xs = [0.13, 0.37, 0.61, 0.85]
    for g_, (y, e) in X_LEV.items():
        ax.plot([0.0, 1.0], [y, y], color=INK, lw=1.4, solid_capstyle="butt")
        ax.text(1.02, y, f"{GROUND_NAME[g_]}   {e}", fontsize=6.6, va="center", color=INK)
    for key, (y, e) in B_LEV.items():
        ax.plot([0.0, 1.0], [y, y], color=INK, lw=1.4, solid_capstyle="butt")
        ax.text(1.02, y, f"{EXC_NAME[key]}   {e}", fontsize=6.6, va="center", color=INK)
    ax.text(0.0, 9.15, "B ³Π₁, J′ = 3, P′ = −   (zero field, energy from the lowest level)", fontsize=7.2,
            color=GRAY, va="center")
    ax.text(0.0, -0.42, "X ¹Σ⁺, J = 2   (zero field, energy from the lowest level; not to scale)", fontsize=7.2,
            color=GRAY, va="center")
    for x, (g_, e), col, (m, lo), tt, b in zip(xs, PAIRS, COLORS, eff_m, eff_t, before):
        y1 = B_LEV[e][0]
        for gg, v in b.items():                       # other levels carrying >= 5 % of the CU2-off signal
            if gg != g_ and v >= 0.05:
                y0 = X_LEV[gg][0]
                ax.add_patch(FancyArrowPatch((x - 0.075, y0), (x - 0.035, y1), arrowstyle="-|>", mutation_scale=7,
                                             lw=1.1, color=col, ls=(0, (3, 2)), alpha=0.9, zorder=4))
                ax.plot([x - 0.11, x - 0.04], [y0, y0], color=col, lw=2.2, alpha=0.55, solid_capstyle="butt",
                        zorder=5)
        y0 = X_LEV[g_][0]
        ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>", mutation_scale=9, lw=1.8, color=col,
                                     zorder=4))
        ax.plot([x - 0.05, x + 0.05], [y0, y0], color=col, lw=3.2, solid_capstyle="butt", zorder=5)
        name = f"R(2)\nF₁′={e.split(',')[0]}\nF′={e.split(',')[1]}"
        ax.text(x, 6.05, name, fontsize=6.5, color=col, ha="center", va="top", weight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none"), zorder=6)
        txt = (f"measured\n{m:.1f} %\n\nmodel\n{tt:.1f} %" if nocl else
               f"measured\n{m:.1f} %\n> {lo:.1f} %\n(90 % CL)\n\nmodel\n{tt:.1f} %")
        ax.text(x, 4.5, txt, fontsize=6.3, color=INK, ha="center", va="top", linespacing=1.15,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=col, lw=0.9), zorder=6)
        if not nocl:                                   # 'CU2 off' signal shares under each box
            comp = "\n".join(f"{SHORT[gg]}: {100 * v:.0f} %" for gg, v in sorted(b.items(), key=lambda q: -q[1])
                             if v >= 0.005)
            ax.text(x, 2.35, "CU2 off\n" + comp, fontsize=5.8, color=GRAY, ha="center", va="top",
                    linespacing=1.1, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"), zorder=6)
    ax.text(-0.02, 9.75, "(a) four detection lines and what each one sees with CU2 off", fontsize=8.5,
            weight="bold", color=INK, transform=ax.transData, va="bottom")

    # ---- (b) efficiency per line
    for i, col in enumerate(COLORS):
        m, lo = eff_m[i]
        axb.errorbar(i - 0.14, m, [[m - lo], [100 - m]], fmt="o", ms=5.5, color=col, ecolor=col, elinewidth=1.2,
                     capsize=0, zorder=4)
        axb.plot(i + 0.14, eff_t[i], "D", ms=5.5, mfc="white", mec=col, mew=1.4, zorder=4)
    axb.errorbar([], [], [], fmt="o", ms=5, color=GRAY,
                 label="measured: 1 − posterior mean of R (R ≥ 0);\nbar: 90 % CL lower limit to 100 %")
    axb.plot([], [], "D", ms=5, mfc="white", mec=GRAY, mew=1.3,
             label=f"model: probe photons per level,\n{power_mW:g} mW, start fitted to CU2 off")
    axb.set_xticks(range(4), ["R(2) 7/2,\nF′=4", "R(2) 7/2,\nF′=3", "R(2) 5/2,\nF′=3", "R(2) 5/2,\nF′=2"])
    axb.set_xlabel("detection line (F₁′, F′)")
    axb.set_xlim(-0.5, 3.5); axb.set_ylim(93.6, 100.6)
    axb.set_yticks([94, 95, 96, 97, 98, 99, 100])
    axb.axhline(100, color=LGRAY, lw=0.8, zorder=0)
    axb.grid(axis="y", color=GRID, lw=0.6, zorder=0)
    axb.set_ylabel("CU2 cleanup efficiency seen by the line (%)")
    axb.legend(frameon=False, fontsize=6.4, loc="lower left")
    axb.text(-0.2, 1.045, "(b) measured vs model", transform=axb.transAxes, fontsize=8.5, weight="bold", color=INK)

    # ---- explanation
    w = t["weights"]
    over = [t["S_fit"][k] / t["S_meas"][k] - 1 for k in ("F1p5F3", "F1p5F2")]
    inside = all(lo <= tt <= 100 for (m, lo), tt in zip(eff_m, eff_t))
    lines = [
        (INK, "CU2 test 2026-11-08: CU2 217 V/cm, 20 + 30 mW, f$_\\mathrm{AOM}$ 209 MHz, NLTL2 2068 (δ₀ = 13 MHz); "
              f"SPB at profile, CU1 off; DET plates off (E = 0), σ+ about k, {power_mW:g} mW"),
        (INK, "efficiency = 1 − S(CU2 on) / S(CU2 off) of the line; measured: 1 − posterior mean of R with R ≥ 0 and the "
              "90 % CL lower limit"),
        (INK, f"model: every J=2 level weighted by its probe photon yield ({power_mW:g} mW, v$_z$ {T_VZ:g} m/s, "
              "vx-averaged); start fitted to the four CU2-off signals: "
              + ", ".join(f"{SHORT[g_]} {100 * v:.0f} %" for g_, v in w.items())),
        (COLORS[2], "→ at zero field the J=2 levels lie within 0.36 MHz and the probe saturates weak couplings, so the "
                    "two F₁ = 3/2 lines see mostly F₁ = 5/2 molecules (dashed arrows, 'CU2 off' shares)"),
        (COLORS[2], f"→ with CU2 on, R(2) F₁′=5/2, F′=3 sees mainly F₁ = 3/2, F = 2 that CU2 pumped in from the lens levels: "
                    f"{100 * after[2][(1.5, 2)]:.1f} % of its {100 - eff_t[2]:.1f} % (F = 3 left: {100 * after[2][(2.5, 3)]:.1f} %)"),
        (COLORS[0], "→ all four model values lie inside the measured 90 % CL ranges" if inside else
                    "→ not every model value lies inside the measured 90 % CL range"),
        (GRAY, f"    caveats: the lens levels alone over-predict the CU2-off signals of the F₁ = 3/2 lines by "
               f"{100 * over[0]:.0f} % and {100 * over[1]:.0f} %, so these weak-coupling yields are uncertain;"
               f" CU2 model at v$_z$ = {t['vz']:g} m/s"),
    ]
    if axt is not None:
        for j, (c, s) in enumerate(lines):
            axt.text(0.0, 0.98 - 0.135 * j, s, fontsize=6.6, color=c, transform=axt.transAxes, va="top")
    fig.suptitle("CU2 test: cleanup efficiency seen by each detection line, measured and simulated", x=0.02,
                 ha="left", fontsize=10, weight="bold", color=INK, y=0.995)
    OUT.mkdir(parents=True, exist_ok=True)
    name = "CU2_test_2026-11-08_efficiency_schematic" + ("_notext" if notext else "") + ("_nocl" if nocl else "")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=220, bbox_inches="tight", pad_inches=0.06)
    print("wrote", OUT / f"{name}.png")


T_VZ = D.T.VZ

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a not in ("notext", "nocl")]
    main(args[0], float(args[1]) if len(args) > 1 else 24, notext="notext" in sys.argv[1:],
         nocl="nocl" in sys.argv[1:])
