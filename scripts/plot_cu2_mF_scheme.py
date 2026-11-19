"""CU2 transitions at 217 V/cm with the mF' sublevels of B J'=3 F1'=7/2 F'=3 and F'=4: which ones the lasers drive.

    python scripts/plot_cu2_mF_scheme.py compute   # builds the 217 V/cm CU2 model, caches sublevel positions (~1 min)
    python scripts/plot_cu2_mF_scheme.py plot      # slide figure
    -> docs/state_cleanup/figures/slides/CU2_mF_scheme_217.{png,pdf}

Positions and strengths come from scripts/cu2_level_scheme.compute (same definitions) run on the CU2 model of
scripts/cu2_laser_aom_scan.build at E = 217 V/cm: for every excited sublevel addressed by a beam, one entry per
transition from the 20 X J=2 sublevels with its resonance on that beam's detuning axis and its |M_eg|^2 summed
over both polarizations (X and Z). The bar is drawn at the coupling-weighted mean resonance (ground Stark spread
<= 1.7 MHz). Beam-1 (F'=4) positions are put on the F'=3-beam axis by adding Delta_main.

Relative excitation rate per excited sublevel, weak field, at the CU2-test carriers delta0 = 13 MHz,
f_AOM = 209 MHz (delta1 = delta0 + f_AOM - Delta_main):
    sum over J=2 sublevels and EOM orders n of |M|^2 J_n(beta)^2 V(pos - line_n),
V = Voigt profile, Lorentzian FWHM Gamma/2pi = 1.61 MHz and Gaussian Doppler sigma = sigma_vx / lambda =
0.65 m/s * 3.68 MHz/(m/s) = 2.39 MHz (transverse velocity spread of the lensed J=2 beam). Power broadening is
not included, so weakly driven levels are if anything underestimated. Rates are relative to the most strongly
driven sublevel of the same beam.
At 217 V/cm the Lambda doublet is Stark mixed, so each mF' appears twice (lower: P'=- parent, upper: P'=+).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
CACHE = ROOT / "tmp" / "cu2_mF_scheme_217" / "level_scheme.npz"
FIG = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
E_V_CM = 217.0
D0, F_AOM = 13.0, 209.0
GAMMA = 1.6076                     # MHz, natural FWHM
SIGMA_D = 0.65 * 3.68              # MHz, Doppler sigma
STRONG, WEAK = 0.20, 0.03          # relative-rate classes


def compute():
    from scripts import cu2_laser_aom_scan as base
    base.E_FIELD = np.array([0.0, 0.0, E_V_CM])          # read by base.build()
    from scripts import cu2_level_scheme as ls
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    ls.CACHE = CACHE
    ls.compute()


def rates(rows, pairs, k, delta, f_eom, beta):
    from scipy.special import jv, voigt_profile
    n = np.arange(-3, 4)
    pw = jv(n, beta) ** 2
    lines = delta + n * f_eom
    r = np.zeros(len(rows))
    for e, pos, w, *_ in pairs:
        e = int(e)
        if rows[e, 0] == k:
            r[e] += w * np.sum(pw * voigt_profile(pos - lines, SIGMA_D, GAMMA / 2))
    sel = rows[:, 0] == k
    return r / r[sel].max()


def plot():
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    from scipy.special import jv
    from scripts.cu2_level_scheme import BETA, F_EOM

    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
                         "mathtext.fontset": "stix"})
    RED, BLUE = "#d1495b", "#1f6fb4"
    z = np.load(CACHE)
    rows, pairs, dmain = z["rows"], z["pairs"], float(z["delta_main"])
    d1 = D0 + F_AOM - dmain
    print(f"Delta_main = {dmain:.3f} MHz, delta1 = {d1:.2f} MHz")

    fig = plt.figure(figsize=(12.5, 8.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[7.4, 9.4], height_ratios=[5.0, 1.25], left=0.075, right=0.93,
                          top=0.855, bottom=0.02, wspace=0.17, hspace=0.12)
    orders = np.arange(-2, 3)
    power = jv(orders, BETA) ** 2
    for col, (k, delta, color, F) in enumerate(((0, D0, RED, 3), (1, d1, BLUE, 4))):
        ax = fig.add_subplot(gs[0, col])
        offset = 0.0 if k == 0 else dmain
        rel = rates(rows, pairs, k, delta, F_EOM, BETA)
        sel = np.flatnonzero(rows[:, 0] == k)
        mmax = int(np.max(np.abs(rows[sel, 1])))
        for n, pw in zip(orders, power):                  # carrier and EOM sidebands of this beam
            y = offset + delta + n * F_EOM
            ax.axhline(y, color=color, lw=0.5 + 5.0 * pw, alpha=0.5, zorder=1)
            ax.text(mmax + 0.75, y, "carrier" if n == 0 else f"{n:+d}", color=color, fontsize=10.5,
                    ha="left", va="center")
        counts = [0, 0]
        for i in sel:
            mF, P, pos = rows[i, 1], rows[i, 2], rows[i, 3] + offset
            r = rel[i]
            if r >= STRONG:
                c, lw = "k", 3.6
                counts[0] += 1
            elif r >= WEAK:
                c, lw = "0.5", 2.4
                counts[1] += 1
            else:
                c, lw = "0.82", 1.4
            ax.plot([mF - 0.36, mF + 0.36], [pos, pos], color=c, lw=lw, solid_capstyle="butt", zorder=3)
            if r >= WEAK:
                ax.text(mF, pos + 1.3, f"{100 * r:.0f}%", fontsize=9, color=c, ha="center", va="bottom", zorder=4)
        for P in (-1.0, 1.0):                             # label the two Lambda-doublet branches at mF' = -mmax
            i = sel[(rows[sel, 1] == -mmax) & (rows[sel, 2] == P)][0]
            ax.text(-mmax - 0.42, rows[i, 3] + offset, "$P'=-$" if P < 0 else "$P'=+$", fontsize=10,
                    color="0.4", ha="right", va="center")
        print(f"F'={F}: strong {counts[0]}, weak {counts[1]}, of {len(sel)}",
              {f"{rows[i,1]:+.0f}{'+' if rows[i,2] > 0 else '-'}": round(100 * rel[i], 1) for i in sel if rows[i, 1] >= 0})
        ys = rows[sel, 3] + offset
        ax.set_ylim(min(ys.min(), offset + delta - 2 * F_EOM) - 5, max(ys.max(), offset + delta + 2 * F_EOM) + 7)
        ax.set_xlim(-mmax - 1.5, mmax + 0.7)
        ax.set_xticks(np.arange(-mmax, mmax + 1))
        ax.tick_params(labelsize=11)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_xlabel("$m_F'$", fontsize=13, labelpad=1)
        ax.set_title(f"$B\\,|J'=3,\\,F_1'=7/2,\\,F'={F}\\rangle$, {'red' if k == 0 else 'blue'} beam", fontsize=13.5)
        if col == 0:
            ax.set_ylabel("transition frequency above the $F'=3$ main line (MHz)", fontsize=12)

    # ground level, beams and AOM (schematic)
    axg = fig.add_subplot(gs[1, :])
    axg.axis("off")
    axg.set_xlim(0, 1)
    axg.set_ylim(0, 1)
    axg.plot([0.22, 0.78], [0.22, 0.22], color="k", lw=3.2, solid_capstyle="butt")
    axg.text(0.5, 0.12, "$X\\,|J=2^+\\rangle$, all 20 sublevels", ha="center", va="top", fontsize=13)
    for xs, xe, color, name, ha in ((0.34, 0.27, RED, f"$F'=3$ beam, $\\delta_0$ = +{D0:g} MHz", "right"),
                                    (0.66, 0.73, BLUE, f"$F'=4$ beam, $\\delta_1$ = +{d1:.1f} MHz", "left")):
        axg.add_patch(FancyArrowPatch((xs, 0.25), (xe, 0.80), arrowstyle="-|>", mutation_scale=22, lw=3.0,
                                      color=color, transform=axg.transAxes))
        axg.text(xe + (-0.02 if ha == "right" else 0.02), 0.62, name, color=color, fontsize=12.5, ha=ha, va="center")
    axg.add_patch(FancyArrowPatch((0.37, 0.55), (0.63, 0.55), arrowstyle="-|>", mutation_scale=16, lw=1.6,
                                  color="0.2", connectionstyle="arc3,rad=-0.15", transform=axg.transAxes))
    axg.text(0.5, 0.70, f"AOM  +$f_{{AOM}}$ ({F_AOM:g} MHz)", ha="center", va="bottom", fontsize=11.5,
             fontweight="bold", color="0.2")

    fig.text(0.075, 0.985, f"CU2 at {E_V_CM:g} V/cm: which $m_F'$ sublevels the two beams talk to",
             ha="left", va="top", fontsize=16, fontweight="bold")
    fig.text(0.075, 0.948, f"%: excitation rate relative to the strongest sublevel of the same beam (from all $J=2$ "
             f"sublevels, X + Z polarization, weak field, Doppler $\\sigma_{{v_x}}$ = 0.65 m/s).  "
             f"Black ≥ {100 * STRONG:g} %, grey {100 * WEAK:g}–{100 * STRONG:g} %, light < {100 * WEAK:g} %.",
             ha="left", va="top", fontsize=10.5, color="0.3")
    fig.text(0.075, 0.922, f"Coloured lines: carrier and EOM sidebands ($f_{{EOM}}$ = {F_EOM:g} MHz, $\\beta$ = {BETA:g}), "
             "width ∝ $J_n(\\beta)^2$.  Each $m_F'$ appears twice: Stark-mixed $\\Lambda$ doublet ($P'=\\pm$ parents).",
             ha="left", va="top", fontsize=10.5, color="0.3")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"CU2_mF_scheme_217.{ext}", dpi=170)
    print(FIG / "CU2_mF_scheme_217.png")


if __name__ == "__main__":
    {"compute": compute, "plot": plot}[sys.argv[1]]()
