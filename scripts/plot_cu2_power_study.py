"""Results of the CU2 power study (60 mW, 5 % loss per pass) in one figure.

    python scripts/plot_cu2_power_study.py   # docs/state_cleanup/figures/CU2_power_study.{pdf,png}

Every number is read from the checkpoints of scripts/cu2_power_check.py in tmp/cu2_power. Panels:
(a) power split, (b) vx dependence on axis, (c) laser detuning, (d) beam height, (e) loss per pass, (f) E field.
(a), (b), (d) are full acceptance grids; (c), (e), (f) use the 6-point proxy (reads 0.1-0.7 % high).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys  # noqa: E402
sys.path.insert(0, str(ROOT))
from scripts import vz_paths  # noqa: E402
sys.path.insert(0, str(ROOT))
from scripts import cu2_power_check as pc  # noqa: E402
from scripts.plot_cleanup_maps_current import SY, acceptance_mean  # noqa: E402

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator  # noqa: E402

OUT = ROOT / "docs" / "state_cleanup" / "figures"
INK, MUTED, GRID, SURFACE = "#262624", "#6b6b68", "#e4e3dc", "#ffffff"
COLOR = {"spb": "#2a78d6", "uniform": "#eb6834"}
LABEL = {"spb": "SPB start ($m_J=0$)", "uniform": "uniform start (all $J=2$)"}
INITS = ("spb", "uniform")
RES_TICKS = [0.03, 0.1, 0.3, 1, 3, 10]              # remaining J=2 (%)

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "axes.labelcolor": INK, "text.color": INK, "savefig.facecolor": SURFACE,
})


def rem(stem):
    parts = vz_paths.glob(pc.OUT, f"{stem}_*.npz") or [vz_paths.find(pc.OUT / f"{stem}.npz")]
    return np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in parts])


def acc_grid(name, p0, p1, init, tag=""):
    return rem(f"acc_loss5_{name}_{p0}_{p1}{tag}_{init}").reshape(len(pc.VX), len(pc.Y0))


def acc_pct(*args, **kw):
    """Remaining J=2 (%) averaged over the detection acceptance."""
    return 100 * acceptance_mean(acc_grid(*args, **kw), pc.VX, pc.Y0, SY["CU2"])


def proxy_pct(r):
    """Remaining J=2 (%) in the 6-point proxy, one value per setting."""
    return 100 * pc.proxy(r)


def grid_proxy(stem, grid):
    d0s, aoms = pc.GRIDS[grid]
    return proxy_pct(rem(stem)).reshape(len(d0s), len(aoms)), d0s, aoms


def style(ax, xlabel, title, note, top=False):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(FixedLocator(RES_TICKS)); ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticklabels([f"{t:g}" for t in RES_TICKS])
    ax.set_xlabel(xlabel); ax.set_title(title, loc="left")
    ax.text(0.5 if top else 0.03, 0.96 if top else 0.04, note, transform=ax.transAxes, fontsize=6.8, color=MUTED,
            va="top" if top else "bottom", ha="center" if top else "left",
            bbox=dict(fc=SURFACE, ec="none", pad=0.8), zorder=2)
    ax.grid(True, color=GRID, lw=0.5, zorder=0); ax.set_axisbelow(True)


def right_ticks(ax):
    """Efficiency labels on the right: the same quantity as the left axis, 100 % - remaining."""
    lo, hi = ax.get_ylim()
    ticks = [t for t in RES_TICKS if lo <= t <= hi]
    ax2 = ax.twinx(); ax2.set_yscale("log"); ax2.set_ylim(lo, hi)
    ax2.yaxis.set_major_locator(FixedLocator(ticks)); ax2.yaxis.set_minor_locator(NullLocator())
    ax2.set_yticklabels([f"{100 - t:g}" for t in ticks])
    for s in ax2.spines.values():
        s.set_visible(False)
    return ax2


def line(ax, x, y, init, ls="-", filled=True, label=None, ms=3.6):
    ax.plot(x, y, color=COLOR[init], lw=2 if ls == "-" else 1.2, ls=ls, marker="o", ms=ms,
            mfc=COLOR[init] if filled else SURFACE, mec=COLOR[init], mew=1.0, label=label, zorder=3)


def main():
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.7), sharey=True)
    ylim = (0.02, 12)

    # (a) power split, full acceptance mean, R1 and the new set point
    ax = axes[0, 0]
    r1 = {"spb": [(20, 40), (25, 35), (30, 30), (35, 25), (40, 20), (45, 15)],
          "uniform": [(20, 40), (25, 35), (30, 30), (35, 25)]}
    new = [(25, 35), (30, 30), (35, 25), (40, 20)]
    for init in INITS:
        line(ax, [p1 for _, p1 in r1[init]], [acc_pct("R1", a, b, init) for a, b in r1[init]], init,
             ls="--", filled=False)
        line(ax, [p1 for _, p1 in new], [acc_pct("d12p5_f207p5", a, b, init) for a, b in new], init)
    ax.axvline(25, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(25.6, 9.0, "present\n35 + 25 mW", fontsize=6.8, color=MUTED, va="top")
    ax.set_xticks([15, 20, 25, 30, 35, 40])
    style(ax, "$P$(F$'$=4) (mW), total 60 mW", "(a) power split", "acceptance mean")

    # (b) vx dependence on the beam axis (y0 = 0), 35 + 25 mW
    ax = axes[0, 1]
    for init in INITS:
        line(ax, pc.VX, 100 * acc_grid("R1", 35, 25, init)[:, 0], init, ls="--", filled=False, ms=3)
        line(ax, pc.VX, 100 * acc_grid("d12p5_f207p5", 35, 25, init)[:, 0], init, ms=3)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    style(ax, "transverse velocity $v_x$ (m/s)", "(b) transverse velocity", "$y_0$ = 0, 35 + 25 mW", top=True)

    # (c) laser detuning at f_AOM = 207.5 MHz (expected in-situ laser scan), proxy
    ax = axes[0, 2]
    for init in INITS:
        d, e = [], []
        for grid, col in (("fine", 1), ("fine2", 3)):        # f_AOM = 207.5 MHz column of each grid
            tab, d0s, aoms = grid_proxy(f"{grid}_loss5_35_25_{init}", grid)
            assert aoms[col] == 207.5
            d += list(d0s); e += list(tab[:, col])
        order = np.argsort(d); d, e = np.asarray(d)[order], np.asarray(e)[order]
        keep = np.r_[True, np.diff(d) > 0]
        line(ax, d[keep], e[keep], init)
    ax.axvline(12.5, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.axvline(11.5, color=MUTED, lw=0.6, ls=":", zorder=1)
    ax.text(12.6, 9.0, "new", fontsize=6.8, color=MUTED, va="top")
    ax.text(11.4, 9.0, "R1", fontsize=6.8, color=MUTED, va="top", ha="right")
    style(ax, r"F$'$=3 laser detuning $\delta_0/2\pi$ (MHz)", "(c) laser detuning",
          r"proxy, $f_\mathrm{AOM}$ = 207.5 MHz")

    # (d) beam height at the new set point, same power, full acceptance mean
    ax = axes[1, 0]
    heights = [(1.6, "_sy4"), (2.0, "_sy5"), (2.4, "_sy6"), (3.0, ""), (4.0, "_sy10")]
    for init in INITS:
        line(ax, [h for h, _ in heights],
             [acc_pct("d12p5_f207p5", 35, 25, init, tag=t) for _, t in heights], init)
    ax.axvline(3.0, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(3.06, 9.0, "present", fontsize=6.8, color=MUTED, va="top")
    style(ax, r"beam height $4\sigma_y$ (cm)", "(d) beam height", "acceptance mean, same power")

    # (e) loss per pass, proxy
    ax = axes[1, 1]
    for init in INITS:
        tab, d0s, aoms = grid_proxy(f"fine2_loss5_35_25_{init}", "fine2")
        r5 = tab[d0s.index(12.5), aoms.index(207.5)]
        loss = [3, 5, 8, 10, 15]
        line(ax, loss, [proxy_pct(rem(f"loss_{L}_{init}"))[0] if L != 5 else r5 for L in loss], init)
    ax.axvline(5, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(5.3, 9.0, "assumed", fontsize=6.8, color=MUTED, va="top")
    style(ax, "window loss per pass (%)", "(e) loss per pass", "proxy")

    # (f) E field with the laser frequencies fixed, proxy
    ax = axes[1, 2]
    fields = [190, 195, 197.5, 200, 202.5, 205, 210]
    n = len(pc.PROXY_VX) * len(pc.PROXY_Y0)
    for init in INITS:
        tab, d0s, aoms = grid_proxy(f"fine2_loss5_35_25_{init}", "fine2")
        r0 = tab[d0s.index(12.5), aoms.index(207.5)]
        line(ax, fields, [r0 if E == 200 else proxy_pct(rem(f"efield_loss5_{E:g}_{init}")[:n])[0]
                          for E in fields], init)
    ax.axvline(200, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(200.4, 9.0, "set point\nfound here", fontsize=6.8, color=MUTED, va="top")
    style(ax, "CU2 field (V/cm), laser fixed", "(f) field error", "proxy")

    for ax in axes.flat:
        ax.set_ylim(*ylim)
    for ax in axes[:, 0]:
        ax.set_ylabel("remaining in $J=2$ (%)")
    for ax in axes[:, 2]:
        right_ticks(ax).set_ylabel("cleanup efficiency (%)")
    fig.tight_layout(rect=(0, 0.115, 1, 1))

    # Legend and the explanation of every quantity, in the figure.
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=COLOR[i], lw=2, marker="o", ms=3.6, label=LABEL[i]) for i in INITS]
    handles += [Line2D([], [], color=INK, lw=2, marker="o", ms=3.6,
                       label=r"new set point $(\delta_0, f_\mathrm{AOM}) = (12.5, 207.5)$ MHz"),
                Line2D([], [], color=INK, lw=1.2, ls="--", marker="o", ms=3.6, mfc=SURFACE,
                       label="R1 = (11.5, 208.5) MHz, (a) and (b) only")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, 0.06), handlelength=2.4, columnspacing=1.6)
    fig.text(0.5, 0.008,
             "Left axis: $J=2$ population remaining after CU2; right axis: cleanup efficiency = 100 % $-$ remaining.\n"
             "Unless varied: 35 + 25 mW, new set point, 13 passes, 5 % loss/pass, 3 cm beam, 200 V/cm, "
             "$v_z$ = 184 m/s.\n"
             r"Acceptance mean: $|v_x| \leq 1$ m/s ($\sigma_{v_x}$ = 0.65 m/s), $|y_0| \leq 8.5$ mm. "
             r"Proxy: 6 points, $v_x$ = 0, ±0.6 m/s × $y_0$ = 0, 6 mm (reads 0.1–0.7 % high).",
             ha="center", va="bottom", fontsize=6.6, color=MUTED, linespacing=1.4)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"CU2_power_study.{ext}", dpi=220)
    print("wrote", OUT / "CU2_power_study.{pdf,png}")


if __name__ == "__main__":
    main()
