"""Plot the CU2 2D cleanup-efficiency maps (laser frequency x AOM frequency, 200 V/cm).

Reads the runs written by scripts/cu2_laser_aom_scan.py from results/SPB_J2_cleanup_scans.h5 and writes
    docs/state_cleanup/figures/CU2_200V_laser_AOM_map.{pdf,png}   coarse grid + fine grids overlaid
    docs/state_cleanup/figures/CU2_200V_fine_regions.{pdf,png}    the three 0.5 MHz grids
Colour encodes the remaining J=2 population on a log scale (darker = better cleanup); the colour
bar is labelled in cleanup efficiency.
"""
from pathlib import Path

import h5py
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
HDF_FILE = ROOT / "results" / "SPB_J2_cleanup_scans.h5"
OUT = ROOT / "docs" / "state_cleanup" / "figures"
PREFIX = "laser_AOM_2D_200V"
EARLIER_AOM = 186.83  # MHz, AOM frequency of the first laser-frequency scan

# Sequential single-hue ramp (light -> dark); darkest = lowest remaining population.
BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
        "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seq_blue_dark_is_good", BLUE[::-1])
NORM = LogNorm(vmin=1e-3, vmax=0.15)               # 99.9 % ... 85 % cleanup
TICKS_REMAINING = [0.1, 0.03, 0.01, 0.005, 0.003, 0.0015]
CONTOURS_REMAINING = [0.01, 0.005, 0.003]          # 99 %, 99.5 %, 99.7 %
INK, MUTED, SURFACE = "#262624", "#6b6b68", "#ffffff"

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "axes.labelcolor": INK, "text.color": INK, "savefig.facecolor": SURFACE,
})


def load(name):
    with h5py.File(HDF_FILE, "r") as h5:
        g = h5["runs"][name]
        return (g["axes"]["Fprime3_detuning_MHz"][...], g["axes"]["AOM_frequency_MHz"][...],
                g["results"][...])


def edges(c):
    c = np.asarray(c, float)
    step = np.median(np.diff(c))
    return np.concatenate((c - step / 2, [c[-1] + step / 2]))


def eff_label(r):
    return f"{100 * (1 - r):.4g} %"


def colorbar(fig, mappable, ax, shrink=1.0):
    cb = fig.colorbar(mappable, ax=ax, pad=0.02, fraction=0.035, aspect=30 * shrink, shrink=shrink)
    cb.set_ticks(TICKS_REMAINING)
    cb.set_ticklabels([eff_label(r) for r in TICKS_REMAINING])
    cb.minorticks_off()
    cb.ax.invert_yaxis()                               # best cleanup at the top
    cb.set_label("cleanup efficiency $\\eta_\\mathrm{c}$")
    cb.outline.set_edgecolor(MUTED); cb.outline.set_linewidth(0.6)
    for r in CONTOURS_REMAINING:
        cb.ax.axhline(r, color=SURFACE, lw=0.8)
    return cb


def cells(ax, d0, aom, rem, lw):
    # x = AOM frequency, y = F'=3 detuning; thin surface-coloured gaps between cells
    return ax.pcolormesh(edges(aom), edges(d0), rem, cmap=CMAP, norm=NORM, shading="flat",
                         edgecolors=SURFACE, linewidth=lw, rasterized=False)


def best(d0, aom, rem):
    i, j = np.unravel_index(np.nanargmin(rem), rem.shape)
    return d0[i], aom[j], rem[i, j]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cd0, caom, crem = load(f"{PREFIX}_coarse_2MHz")
    regions = [load(f"{PREFIX}_fine_0p5MHz_region{k}") for k in (1, 2, 3)]

    # ---- Figure 1: full map --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.1, 3.7), constrained_layout=True)
    m = cells(ax, cd0, caom, crem, lw=0.6)
    for k, (d0, aom, rem) in enumerate(regions, start=1):
        cells(ax, d0, aom, rem, lw=0.25)
        xe, ye = edges(aom), edges(d0)
        ax.add_patch(Rectangle((xe[0], ye[0]), xe[-1] - xe[0], ye[-1] - ye[0], fill=False,
                               edgecolor=INK, lw=0.9, zorder=5))
        bd0, baom, brem = best(d0, aom, rem)
        ax.plot(baom, bd0, "o", ms=5, mfc=SURFACE, mec=INK, mew=0.9, zorder=6)
        label = f"R{k}: {100*(1-brem):.2f} %\n({bd0:.1f}, {baom:.1f}) MHz"
        # labels sit just outside the region box, on the surface side
        if k == 3:
            xy, ha, va = (xe[0] - 0.4, ye[-1] - 0.2), "right", "top"
        else:
            xy, ha, va = (xe[0] + (xe[-1] - xe[0]) / 2, ye[0] - 0.35), "center", "top"
        ax.annotate(label, xy=xy, ha=ha, va=va, fontsize=7.5, color=INK, zorder=7,
                    bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none", alpha=0.9))
    ax.axvline(EARLIER_AOM, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=4)
    ax.text(EARLIER_AOM - 0.25, 4.4, "first laser-frequency scan, $f_\\mathrm{AOM}=186.83$ MHz",
            rotation=90, fontsize=6.8, color=MUTED, va="bottom", ha="right", zorder=7,
            bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec="none", alpha=0.85))
    ax.set_xlim(184, 216); ax.set_ylim(4, 25.25)
    ax.set_xlabel("AOM frequency $f_\\mathrm{AOM}$ (MHz)")
    ax.set_ylabel("$F'=3$ detuning $\\delta_0/2\\pi$ (MHz)")
    ax.set_xticks(np.arange(185, 216, 5)); ax.set_yticks(np.arange(5, 26, 5))
    ax.tick_params(length=2.5, width=0.6)
    ax.set_facecolor("#f3f3f1")                     # not computed (above the coarse grid)
    colorbar(fig, m, ax)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"CU2_200V_laser_AOM_map.{ext}", dpi=220)
    plt.close(fig)

    # ---- Figure 2: the three fine grids ------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(6.1, 2.55), constrained_layout=True)
    for k, (ax, (d0, aom, rem)) in enumerate(zip(axes, regions), start=1):
        m = cells(ax, d0, aom, rem, lw=0.4)
        ax.contour(aom, d0, rem, levels=sorted(CONTOURS_REMAINING), colors=SURFACE, linewidths=0.8)
        bd0, baom, brem = best(d0, aom, rem)
        ax.plot(baom, bd0, "o", ms=5, mfc=SURFACE, mec=INK, mew=0.9, zorder=6)
        ax.set_title(f"Region {k}: best {100*(1-brem):.2f} %\nat ({bd0:.1f}, {baom:.1f}) MHz",
                     fontsize=7.5, color=INK)
        ax.set_xlabel("$f_\\mathrm{AOM}$ (MHz)")
        if k == 1:
            ax.set_ylabel("$\\delta_0/2\\pi$ (MHz)")
        ax.set_xticks(aom[::2]); ax.set_yticks(d0[::2])
        ax.tick_params(length=2.5, width=0.6, labelsize=7)
        ax.set_aspect("equal")
    colorbar(fig, m, axes, shrink=0.78)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"CU2_200V_fine_regions.{ext}", dpi=220)
    plt.close(fig)
    print("written:", *sorted(p.name for p in OUT.glob("CU2_200V_*")))


if __name__ == "__main__":
    main()
