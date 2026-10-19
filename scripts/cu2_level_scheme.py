"""Level diagrams of CU2 at 200 V/cm with the laser frequency combs of the R1 and R2 optima.

    python scripts/cu2_level_scheme.py compute   # builds the model, caches line positions (~1 min)
    python scripts/cu2_level_scheme.py paper     # report figure, CeNTREX-paper style
    python scripts/cu2_level_scheme.py slides    # schematic for slides (short-IR slide style)

For every excited sublevel addressed by a beam, the position is the coupling-weighted mean resonance
of its transitions from the 20 X, J=2 sublevels, on that beam's own detuning axis (MHz from its main
resonance); the strength is sum_g sum_pol |M_eg|^2 relative to the strongest sublevel of that beam.
A sublevel counts as driven when its weak-field excitation rate (excitation_rates) is >= 5 % of the
most strongly driven sublevel of the same beam.
Writes docs/state_cleanup/figures/CU2_200V_level_diagram_R1_R2.{pdf,png}            (paper)
and    docs/state_cleanup/figures/slides/CU2_200V_level_scheme_R1_R2_slides.{pdf,png} (slides).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cu2_laser_aom_scan import BETA, CHECKPOINT_DIR, EOM_FREQUENCY, build  # noqa: E402

CACHE = CHECKPOINT_DIR / "level_scheme.npz"
OUT = ROOT / "docs" / "state_cleanup" / "figures"
F_EOM = EOM_FREQUENCY / (2 * np.pi * 1e6)          # 17.5 MHz
REGIONS = {"R1": (11.5, 208.5), "R2": (11.5, 191.0)}


def compute():
    import sympy as smp
    from centrex_tlf import states
    prepared, rho0, ground_idx, info, obe = build(return_obe=True)
    QN, H = obe.QN, obe.H_symbolic
    deltas = [s for s in H.free_symbols if str(s).startswith("δ")]
    diag = {i: smp.expand(H[i, i]) for i in range(len(QN))}
    rows, pairs = [], []
    for k, c in enumerate(obe.couplings):
        mat = sum(np.abs(f.field) ** 2 for f in c.fields)
        excited = sorted({e for g in ground_idx for e in np.nonzero(mat[:, g] > 1e-4 * mat.max())[0]
                          if QN[e].largest.electronic_state == states.ElectronicState.B})
        for e in excited:
            pos, w, glab = [], [], []
            for g in ground_idx:
                if mat[e, g] <= 0:
                    continue
                gl = QN[g].largest
                glab.append((gl.F1, gl.F, gl.mF))
                dd = diag[e] - diag[g]
                coeff = next(float(dd.coeff(d)) for d in deltas if float(dd.coeff(d)) != 0)
                const = complex(smp.N(dd.subs({d: 0 for d in deltas}))).real
                pos.append(-const / coeff / (2 * np.pi * 1e6)); w.append(mat[e, g])
            # frame of beam 0: beam-1 lines are offset by Delta_main
            pos = np.asarray(pos) - (info["delta_main_MHz"] if k == 1 else 0.0)
            w = np.asarray(w)
            lab = QN[e].largest
            rows.append((k, lab.mF, lab.P, np.average(pos, weights=w), w.sum(), pos.max() - pos.min()))
            pairs += [(len(rows) - 1, p_, w_, *gl_) for p_, w_, gl_ in zip(pos, w, glab)]
    rows, pairs = np.asarray(rows, dtype=float), np.asarray(pairs, dtype=float)
    for k in (0, 1):
        sel = rows[:, 0] == k
        scale = rows[sel, 4].max()
        pairs[np.isin(pairs[:, 0], np.flatnonzero(sel)), 2] /= scale
        rows[sel, 4] /= scale
    np.savez(CACHE, rows=rows, pairs=pairs, delta_main=info["delta_main_MHz"])
    for k in (0, 1):
        level = "F'=3" if k == 0 else "F'=4"
        print(f"beam {k} ({level}): mF', P', position (MHz), strength, spread")
        for r in rows[rows[:, 0] == k]:
            print(f"   {r[1]:+.0f}  {'+' if r[2] > 0 else '-'}  {r[3]:7.2f}  {r[4]:.3f}  {r[5]:.2f}")


def excitation_rates(rows, pairs, k, delta, gsel=None, gamma=1.6076):
    """Weak-field excitation-rate proxy per excited sublevel (row index -> rate).

    sum over ground sublevels and comb lines of |M_eg|^2 * J_n(beta)^2 * Lorentzian(FWHM = Gamma/2pi).
    """
    from scipy.special import jv
    n = np.arange(-3, 4)
    power = jv(n, BETA) ** 2
    lines = delta + n * F_EOM
    rate = np.zeros(len(rows))
    for e, pos, w, F1, F, mF in pairs:
        e = int(e)
        if rows[e, 0] != k or (gsel is not None and not gsel(F1, F, mF)):
            continue
        rate[e] += w * np.sum(power / (1 + (2 * (pos - lines) / gamma) ** 2))
    return rate


def plot_slides():
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.special import jv

    # style of the short-IR thesis figures: sans serif, bold panel titles, light grid
    INK, GRAY, LGRAY, GRID = "#111111", "#555555", "#c8c8c8", "#e6e6e6"
    BLUE, ORANGE, GREEN, YELLOW = "#2a78d6", "#e8673a", "#2e8b57", "#f6e3a1"
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "mathtext.fontset": "dejavusans",
        "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "axes.linewidth": 0.7, "axes.edgecolor": GRAY, "axes.labelcolor": INK,
        "xtick.color": GRAY, "ytick.color": GRAY, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
        "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": "white",
    })
    z = np.load(CACHE)
    rows, pairs, dmain = z["rows"], z["pairs"], float(z["delta_main"])
    orders = np.arange(-3, 4)
    power = jv(orders, BETA) ** 2
    d1 = {r: REGIONS[r][0] + REGIONS[r][1] - dmain for r in REGIONS}
    ADDRESSED = 0.05                                      # relative rate above which a level is "driven"

    panels = [
        (0, REGIONS["R1"][0], BLUE, "(a) F′ = 3, beam 0", ["R1 and R2", "δ₀/2π = +11.5 MHz"], (-45, 95)),
        (1, d1["R1"], ORANGE, "(b) F′ = 4, beam 1, R1",
         ["$f_\\mathrm{AOM}$ = 208.5 MHz", f"δ₁/2π = {d1['R1']:+.2f} MHz".replace("-", "−")], (-65, 65)),
        (1, d1["R2"], ORANGE, "(c) F′ = 4, beam 1, R2",
         ["$f_\\mathrm{AOM}$ = 191.0 MHz", f"δ₁/2π = {d1['R2']:+.2f} MHz".replace("-", "−")], (-65, 65)),
    ]
    fig = plt.figure(figsize=(6.1, 6.05))
    gs = fig.add_gridspec(2, 3, height_ratios=[4.4, 1.95], width_ratios=[7, 9, 9],
                          left=0.09, right=0.955, top=0.885, bottom=0.01, wspace=0.30, hspace=0.27)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    axes[2].sharey(axes[1])
    for ax, (k, delta, color, title, subtitle, ylim) in zip(axes, panels):
        rate = excitation_rates(rows, pairs, k, delta)
        sel = np.flatnonzero(rows[:, 0] == k)
        rmax = rate[sel].max()
        ax.set_ylim(*ylim)
        ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
        if k == 1:                                         # lines shared by R1 and R2
            for y in (d1["R1"], d1["R2"]):
                ax.axhspan(y - 1.6, y + 1.6, color=YELLOW, zorder=0.5, lw=0)
        for n, pw in zip(orders, power):
            y = delta + n * F_EOM
            if ylim[0] < y < ylim[1]:
                ax.axhline(y, color=color, lw=0.35 + 6.5 * pw, alpha=0.9, zorder=1)
                ax.text(1.005, y, "0" if n == 0 else f"{n:+d}".replace("-", "−"), transform=ax.get_yaxis_transform(),
                        fontsize=6.5, color=color, va="center", ha="left", clip_on=False)
        mmax = int(np.max(np.abs(rows[sel, 1])))
        for i in sel:
            mF, P, pos = rows[i, 1], rows[i, 2], rows[i, 3]
            r = rate[i] / rmax
            driven = r >= ADDRESSED
            shade = INK if driven else ("#a9a9a9" if r >= 5e-3 else "#d6d6d6")
            ax.plot([mF - 0.36, mF + 0.36], [pos, pos], color=shade, lw=3.2 if driven else 1.8,
                    solid_capstyle="butt", zorder=3)
            if driven and mF >= 0:
                ax.text(mF + 0.42, pos + 1.2, f"{100*r:.0f} %", fontsize=6.3, color=INK, va="bottom",
                        ha="left", zorder=4)
        # branch labels: same-mF levels of the two Lambda-doublet components are stacked vertically
        for P in (-1.0, 1.0):
            end = sel[(rows[sel, 1] == -mmax) & (rows[sel, 2] == P)][0]
            ax.text(-mmax + 0.45, rows[end, 3], "P′ = −" if P < 0 else "P′ = +", fontsize=6.5, color=GRAY,
                    ha="left", va="center")
        ax.set_xlim(-mmax - 0.7, mmax + 0.7)
        ax.set_xticks(np.arange(-mmax, mmax + 1))
        ax.set_xlabel("m_F′".replace("m_F′", "$m_F'$"))
        ax.text(0.0, 1.115, title, transform=ax.transAxes, fontsize=8.5, fontweight="bold", color=INK, ha="left")
        for j, s in enumerate(subtitle):
            ax.text(0.0, 1.065 - 0.045 * j, s, transform=ax.transAxes, fontsize=7, color=GRAY, ha="left")
    axes[0].set_ylabel("detuning from the beam's main resonance (MHz)")
    axes[2].tick_params(labelleft=False)

    # stretched-state pumping rates through beam 1, to quote in panel (d)
    stretched = lambda F1, F, mF: F1 == 2.5 and F == 3 and abs(mF) == 3
    rs = {r: excitation_rates(rows, pairs, 1, d1[r], stretched).sum() for r in ("R1", "R2")}

    axd = fig.add_subplot(gs[1, :]); axd.axis("off")
    axd.text(0.0, 1.0, "(d) why both AOM frequencies work", fontsize=8.5, fontweight="bold", color=INK,
             transform=axd.transAxes, va="top")
    lines = [
        (INK, "$f_\\mathrm{AOM}$(R1) − $f_\\mathrm{AOM}$(R2) = 208.5 − 191.0 MHz = 17.5 MHz = $f_\\mathrm{EOM}$"),
        (GREEN, "→ the beam-1 lines lie at the same frequencies in R1 and R2; only their order n shifts by one"),
        (ORANGE, f"→ carrier $J_0^2$(1.5) = {power[3]:.2f} ≈ first order $J_1^2$(1.5) = {power[4]:.2f}, "
                 "so $m_F'$ = 0 and ±2 (yellow) get almost the same power"),
        (BLUE, "→ the bottleneck |F = 3, $m_F$ = ±3⟩ reaches F′ = 4 only via $m_F'$ = ±3, ±4; "
               f"R2's −1 order at {d1['R2'] - F_EOM:.1f} MHz".replace("-2", "−2")),
        (BLUE, f"    pumps it {100*(rs['R2']/rs['R1']-1):.0f} % faster than R1, hence 99.85 % (R2) against "
               "99.82 % (R1)"),
        (GRAY, "bars: excited sublevels, same $m_F'$ stacked; black = driven; line width ∝ sideband power"),
        (GRAY, "%: excitation rate relative to the most strongly driven sublevel of the same beam; black if ≥ 5 %"),
        (GRAY, "    rate = weak-field estimate, Σ over J = 2 sublevels and sidebands of |M|² × sideband power × "
               "Lorentzian (width Γ)"),
    ]
    for j, (c, s) in enumerate(lines):
        axd.text(0.0, 0.84 - 0.121 * j, s, fontsize=6.9 if c != GRAY else 6.3, color=c,
                 transform=axd.transAxes, va="top")
    (OUT / "slides").mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / "slides" / f"CU2_200V_level_scheme_R1_R2_slides.{ext}", dpi=220,
                    bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print("stretched-state pumping rate R2/R1 =", rs["R2"] / rs["R1"])
    print("written:", *sorted(p.name for p in (OUT / "slides").glob("CU2_200V_level_scheme*")))


def plot_paper():
    """Report figure in the style of the level diagrams of the CeNTREX overview paper
    (Grasdijk et al. 2021, Figs. 2-4): serif fonts, monochrome levels with one dash per sublevel,
    driven states in bold, a broken frequency axis between the two excited hyperfine levels.
    Only the laser lines are coloured, one colour per beam.
    """
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.special import jv

    INK, LIGHT = "#111111", "#a8a8a8"
    BLUE, ORANGE = "#1f5aa6", "#d4652f"             # short-IR thesis palette (figures_src/style.py)
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True, "font.size": 9, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "xtick.direction": "out", "ytick.direction": "out", "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "pdf.fonttype": 42, "savefig.facecolor": "white",
    })
    z = np.load(CACHE)
    rows, pairs, dmain = z["rows"], z["pairs"], float(z["delta_main"])
    orders = np.arange(-3, 4)
    power = jv(orders, BETA) ** 2
    d1 = {r: REGIONS[r][0] + REGIONS[r][1] - dmain for r in REGIONS}
    YLIM = {0: (-45, 95), 1: (160, 275)}             # beam-0 frame: MHz above the F'=3 main resonance

    fig = plt.figure(figsize=(6.1, 5.3))
    gs = fig.add_gridspec(2, 2, height_ratios=[YLIM[1][1] - YLIM[1][0], YLIM[0][1] - YLIM[0][0]],
                          left=0.105, right=0.935, top=0.90, bottom=0.085, hspace=0.05, wspace=0.07)
    for col, reg in enumerate(("R1", "R2")):
        ax4 = fig.add_subplot(gs[0, col])
        ax3 = fig.add_subplot(gs[1, col], sharex=ax4)
        for ax, k, delta, color in ((ax3, 0, REGIONS[reg][0], BLUE), (ax4, 1, d1[reg], ORANGE)):
            offset = 0.0 if k == 0 else dmain
            rate = excitation_rates(rows, pairs, k, delta)
            sel = np.flatnonzero(rows[:, 0] == k)
            rmax = rate[sel].max()
            ylim = YLIM[k]
            for n, pw in zip(orders, power):            # carrier and sidebands of the beam
                y = offset + delta + n * F_EOM
                if ylim[0] + 2 < y < ylim[1] - 2:
                    ax.axhline(y, color=color, lw=0.3 + 3.2 * pw, zorder=1)
                    ax.text(4.78, y + 0.8, ("$n=0$" if n == 0 else f"${n:+d}$"), color=color, fontsize=6.5,
                            ha="right", va="bottom", zorder=2)
            mmax = int(np.max(np.abs(rows[sel, 1])))
            for i in sel:                                # one dash per excited sublevel
                mF, P, pos = rows[i, 1], rows[i, 2], rows[i, 3] + offset
                driven = rate[i] / rmax >= 0.05
                ax.plot([mF - 0.3, mF + 0.3], [pos, pos], color=INK if driven else LIGHT,
                        lw=2.0 if driven else 0.9, solid_capstyle="butt", zorder=3)
            for P in (-1.0, 1.0):                        # the two Lambda-doublet branches
                end = sel[(rows[sel, 1] == -mmax) & (rows[sel, 2] == P)][0]
                x, ha = (-mmax - 0.4, "right") if mmax < 4 else (-mmax + 0.38, "left")
                ax.text(x, rows[end, 3] + offset, "$P'=-$" if P < 0 else "$P'=+$", fontsize=7, color=INK,
                        ha=ha, va="center")
            ax.set_ylim(*ylim)
            ax.set_xlim(-4.8, 4.8)
            if col == 1:
                ax.text(1.02, 0.5, f"$F'={3 if k == 0 else 4}$", transform=ax.transAxes, rotation=90,
                        ha="left", va="center", fontsize=9)
                ax.tick_params(labelleft=False)
        # broken axis between the two excited hyperfine levels
        ax4.spines["bottom"].set_visible(False)
        ax4.tick_params(bottom=False, labelbottom=False)
        kw = dict(color=INK, lw=0.6, clip_on=False)
        ax4.plot([-0.012, 0.012], [-0.012, 0.012], transform=ax4.transAxes, **kw)
        ax3.plot([-0.012, 0.012], [0.988, 1.012], transform=ax3.transAxes, **kw)
        ax3.set_xticks(np.arange(-4, 5))
        ax3.set_xlabel("$m_F'$")
        title = (f"({'ab'[col]}) {reg}: $f_\\mathrm{{AOM}}={REGIONS[reg][1]:.1f}$ MHz\n"
                 f"$\\delta_0/2\\pi=+{REGIONS[reg][0]:.1f}$ MHz, $\\delta_1/2\\pi={d1[reg]:+.2f}$ MHz")
        ax4.set_title(title, loc="left", fontsize=8.5)
    fig.text(0.015, 0.5, "frequency above the $F'=3$ main resonance (MHz)", rotation=90, ha="left",
             va="center", fontsize=9)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"CU2_200V_level_diagram_R1_R2.{ext}", dpi=300 if ext == "png" else None)
    plt.close(fig)
    print("written:", *sorted(p.name for p in OUT.glob("CU2_200V_level_diagram*")))


if __name__ == "__main__":
    {"compute": compute, "paper": plot_paper, "slides": plot_slides}[sys.argv[1]]()
