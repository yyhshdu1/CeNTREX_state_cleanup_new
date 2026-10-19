"""Plot the CU1 cleanup versus transverse velocity vx and transverse position y0.

Reads the runs doppler_cu1_spb and doppler_cu1_uniform written by scripts/cleanup_spb_doppler.py from
results/SPB_J1_cleanup_scans.h5 and writes docs/state_cleanup/figures/CU1_doppler_vx_y_map.{pdf,png}:
    (a), (b)  remaining J=1 population versus (vx, y0) for the SPB-related and the uniform initial state,
              log colour scale (darker = better cleanup), mirrored to y0 < 0 (the problem is symmetric)
    (c)       cuts versus vx, with the per-pass Doppler shift k*vx on the top axis
Also prints averages over a few transverse distributions.
"""
from pathlib import Path

import h5py
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

ROOT = Path(__file__).resolve().parents[1]
HDF_FILE = ROOT / "results" / "SPB_J1_cleanup_scans.h5"
OUT = ROOT / "docs" / "state_cleanup" / "figures"
DOPPLER_MHZ_PER_MPS = 3.68          # k*vx / 2pi at 271.7 nm

BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
        "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seq_blue_dark_is_good", BLUE[::-1])
NORM = LogNorm(vmin=2e-4, vmax=0.5)
TICKS_REMAINING = [0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 3e-4]
CONTOURS_REMAINING = [0.1, 0.01, 0.001]           # 90 %, 99 %, 99.9 %
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]        # categorical slots 1-3
INK, MUTED, GRID, SURFACE = "#262624", "#6b6b68", "#e4e3dc", "#ffffff"

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
        return g["axes"]["vx_m_per_s"][...], g["axes"]["y0_mm"][...], g["results"][...]


def mirror(y, r):
    """Extend y0 >= 0 to the full symmetric range."""
    return np.concatenate((-y[:0:-1], y)), np.concatenate((r[:, :0:-1], r), axis=1)


def edges(c):
    c = np.asarray(c, float)
    step = np.median(np.diff(c))
    return np.concatenate((c - step / 2, [c[-1] + step / 2]))


def eff_label(r):
    return f"{100 * (1 - r):.4g} %"


def averages(name, vx, y, r):
    """Mean remaining population for a few transverse distributions (y uniform within |y0| <= 15 mm)."""
    yf, rf = mirror(y, r)
    wy = np.ones_like(yf)
    out = {"uniform vx in [-1, 1] m/s": np.average(np.average(rf, axis=1, weights=wy))}
    for s in (0.25, 0.5):
        wv = np.exp(-0.5 * (vx / s) ** 2)
        out[f"Gaussian vx, sigma = {s} m/s (cut at 1 m/s)"] = np.average(np.average(rf, axis=1, weights=wy),
                                                                            weights=wv)
    print(f"{name}: on axis, vx = 0: {eff_label(r[np.argmin(abs(vx)), 0])}")
    for k, v in out.items():
        print(f"   {k}: remaining {v:.3e}  ->  {eff_label(v)}")


def main():
    runs = {"spb": load("doppler_cu1_spb"), "uniform": load("doppler_cu1_uniform")}
    for k, v in runs.items():
        averages(k, *v)

    fig = plt.figure(figsize=(6.3, 6.0))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.05], height_ratios=[1, 0.9], hspace=0.75, wspace=0.12)
    titles = {"spb": "(a) SPB-related levels", "uniform": "(b) all 12 $J=1$ sublevels"}
    for col, key in enumerate(("spb", "uniform")):
        vx, y, r = runs[key]
        yf, rf = mirror(y, r)
        ax = fig.add_subplot(gs[0, col])
        mesh = ax.pcolormesh(edges(vx), edges(yf), rf.T, cmap=CMAP, norm=NORM, shading="flat",
                             edgecolors=SURFACE, linewidth=0.3)
        levels = np.log10(sorted(CONTOURS_REMAINING))          # contour the log, as the colour scale
        cs = ax.contour(vx, yf, np.log10(rf.T), levels=levels, colors=SURFACE, linewidths=0.8, linestyles="solid")
        ax.clabel(cs, fmt={lv: eff_label(10 ** lv) for lv in levels}, fontsize=6.5, inline_spacing=2)
        ax.set_title(titles[key], loc="left", fontsize=8.5)
        ax.set_xlabel("transverse velocity $v_x$ (m/s)")
        if col == 0:
            ax.set_ylabel("transverse position $y_0$ (mm)")
        else:
            ax.tick_params(labelleft=False)
        ax.set_xticks([-1, -0.5, 0, 0.5, 1])
        ax.set_yticks([-15, -10, -5, 0, 5, 10, 15])
    cax = fig.add_subplot(gs[0, 2])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_ticks(TICKS_REMAINING)
    cb.set_ticklabels([eff_label(t) for t in TICKS_REMAINING])
    cb.minorticks_off()
    cb.ax.invert_yaxis()
    cb.set_label("cleanup efficiency $\\eta_\\mathrm{c}$")
    cb.outline.set_edgecolor(MUTED); cb.outline.set_linewidth(0.6)
    for c in CONTOURS_REMAINING:
        cb.ax.axhline(c, color=SURFACE, lw=0.8)

    ax = fig.add_subplot(gs[1, :2])
    cuts = [("spb", 0.0, "SPB-related, $y_0=0$"), ("spb", 12.0, "SPB-related, $y_0=12$ mm"),
            ("uniform", 0.0, "all 12 sublevels, $y_0=0$")]
    for (key, y0, label), color in zip(cuts, SERIES):
        vx, y, r = runs[key]
        rr = r[:, np.argmin(abs(y - y0))]
        ax.plot(vx, rr, color=color, lw=2, marker="o", ms=3.2, label=label, zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(1.5e-4, 0.6)
    ax.set_xlim(-1.05, 1.05)
    ticks = [0.3, 0.1, 0.01, 0.001, 3e-4]
    ax.set_yticks(ticks)
    ax.set_yticklabels([eff_label(t) for t in ticks])
    ax.minorticks_off()
    ax.set_ylabel("cleanup efficiency $\\eta_\\mathrm{c}$")
    ax.invert_yaxis()                                   # better cleanup up
    ax.set_xlabel("transverse velocity $v_x$ (m/s)")
    ax.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("(c) cuts at fixed $y_0$; top axis: Doppler shift per pass, alternating in sign",
                 loc="left", fontsize=8.5, pad=18)
    top = ax.secondary_xaxis("top", functions=(lambda v: v * DOPPLER_MHZ_PER_MPS,
                                               lambda f: f / DOPPLER_MHZ_PER_MPS))
    top.set_xlabel("$\\pm k v_x/2\\pi$ (MHz)", labelpad=2)
    top.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.legend(loc="lower center", frameon=False, fontsize=7.5, ncol=1, handlelength=1.6)
    fig.savefig(OUT / "CU1_doppler_vx_y_map.pdf", bbox_inches="tight")
    fig.savefig(OUT / "CU1_doppler_vx_y_map.png", dpi=200, bbox_inches="tight")
    print("wrote CU1_doppler_vx_y_map.{pdf,png}")


if __name__ == "__main__":
    main()
