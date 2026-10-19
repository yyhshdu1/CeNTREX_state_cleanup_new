"""Plot CU1 with the 1.6 MHz broadening EOM against the transverse Doppler shift.

Writes docs/state_cleanup/figures/CU1_eom_doppler.{pdf,png}:
    (a) vx-averaged cleanup at y0 = 0 versus EOM modulation depth beta and laser detuning
        (vx Gaussian, FWHM 1.8 m/s, truncated at |vx| <= 1 m/s; checkpoints of cu1_eom_doppler_scan.py coarse)
    (b) cleanup versus vx at y0 = 0 for three settings
    (c), (d) vx x y0 maps without and with the EOM (runs doppler_cu1_spb and the EOM map in
        results/SPB_J1_cleanup_scans.h5), same colour scale as CU1_doppler_vx_y_map
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.plot_cu1_doppler_map import (CMAP, CONTOURS_REMAINING, GRID, INK, MUTED, NORM, OUT, SERIES,  # noqa: E402
                                          SURFACE, TICKS_REMAINING, edges, eff_label, load, mirror)
from scripts.cu1_eom_doppler_scan import VX_COARSE, VX_FWHM, key, vx_weights, y_weights  # noqa: E402
from scripts.cu1_eom_doppler_scan import load as load_ckpt  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

EOM_RUN = "doppler_cu1_spb_eom1p6MHz_beta1p5_det20MHz"
BEST = (1.5, 20.0)


def averaged(vx, y, r):
    """(vx-averaged at y0 = 0, vx- and y-averaged) with the weights of cu1_eom_doppler_scan."""
    wv, wy = vx_weights(vx), y_weights(y)
    return wv @ r[:, 0], wv @ r @ wy


def main():
    done = load_ckpt("coarse")
    betas = sorted({k[1] for k in done}); dets = sorted({k[0] for k in done})
    wv = vx_weights(VX_COARSE)
    A = np.full((len(betas), len(dets)), np.nan)
    for i, b in enumerate(betas):
        for j, d in enumerate(dets):
            try:
                A[i, j] = wv @ np.array([done[key((d, b, v, 0.0))] for v in VX_COARSE])
            except KeyError:
                pass

    vx0, y0, r0 = load("doppler_cu1_spb")
    vx1, y1, r1 = load(EOM_RUN)
    for name, (vx, y, r) in (("no EOM, 20.2 MHz", (vx0, y0, r0)), ("EOM beta 1.5, 20.0 MHz", (vx1, y1, r1))):
        a0, a1 = averaged(vx, y, r)
        print(f"{name}: on axis {eff_label(r[np.argmin(abs(vx)), 0])}; vx-avg y0=0 {eff_label(a0)}; "
              f"vx+y avg {eff_label(a1)}")

    fig = plt.figure(figsize=(6.3, 6.6))
    top = fig.add_gridspec(1, 4, width_ratios=[1, 0.05, 0.8, 1.25], wspace=0.06, top=0.97, bottom=0.62)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.05], wspace=0.12, top=0.47, bottom=0.09)

    # (a) beta x detuning
    ax = fig.add_subplot(top[0, 0])
    norm_a = LogNorm(vmin=0.015, vmax=0.25)
    # irregular steps: draw each cell between midpoints
    dm = np.concatenate(([dets[0] - 0.25], 0.5 * (np.array(dets[1:]) + np.array(dets[:-1])), [dets[-1] + 0.25]))
    bm = np.concatenate(([betas[0] - 0.25], 0.5 * (np.array(betas[1:]) + np.array(betas[:-1])), [betas[-1] + 0.25]))
    mesh_a = ax.pcolormesh(dm, bm, np.ma.masked_invalid(A), cmap=CMAP, norm=norm_a, shading="flat",
                           edgecolors=SURFACE, linewidth=0.4)
    ax.set_facecolor("#f1f0ea")
    i, j = np.unravel_index(np.nanargmin(A), A.shape)
    ax.plot(dets[j], betas[i], marker="o", ms=5, mfc="none", mec=SURFACE, mew=1.2)
    ax.plot(BEST[1], BEST[0], marker="s", ms=5, mfc="none", mec="#eda100", mew=1.4)
    ax.plot(20.2, 0.0, marker="^", ms=5, mfc="none", mec="#eb6834", mew=1.4)
    ax.set_xlabel("laser detuning $\\delta/2\\pi$ (MHz)")
    ax.set_ylabel("EOM modulation depth $\\beta$")
    ax.set_title("(a) $v_x$-averaged, $y_0=0$", loc="left", fontsize=8.5)
    cb = fig.colorbar(mesh_a, cax=fig.add_subplot(top[0, 1]))
    ticks = [0.2, 0.1, 0.05, 0.03, 0.02]
    cb.set_ticks(ticks); cb.set_ticklabels([eff_label(t) for t in ticks]); cb.minorticks_off()
    cb.ax.invert_yaxis(); cb.ax.tick_params(labelsize=7); cb.ax.set_title(r"$\eta_\mathrm{c}$", fontsize=8)
    cb.outline.set_edgecolor(MUTED); cb.outline.set_linewidth(0.6)

    # (b) cuts versus vx
    ax = fig.add_subplot(top[0, 3])
    cuts = [(vx0, r0[:, 0], "no EOM, $\\delta/2\\pi=20.2$ MHz", SERIES[1], "^"),
            (VX_COARSE, np.array([done[key((19.0, 0.0, v, 0.0))] for v in VX_COARSE]),
             "no EOM, $\\delta/2\\pi=19.0$ MHz", SERIES[2], "D"),
            (vx1, r1[:, 0], "EOM $\\beta=1.5$, $\\delta/2\\pi=20.0$ MHz", SERIES[0], "s")]
    for vx, rr, label, color, mk in cuts:
        ax.plot(vx, rr, color=color, lw=2, marker=mk, ms=3.2, label=label, zorder=3)
    ax.set_yscale("log"); ax.set_ylim(1.5e-4, 0.4); ax.invert_yaxis()
    yt = [0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 3e-4]
    ax.set_yticks(yt); ax.set_yticklabels([eff_label(t) for t in yt]); ax.minorticks_off()
    ax.set_xlabel("transverse velocity $v_x$ (m/s)"); ax.set_ylabel("cleanup efficiency $\\eta_\\mathrm{c}$")
    ax.set_title("(b) $y_0=0$", loc="left", fontsize=8.5)
    ax.grid(True, color=GRID, lw=0.5, zorder=0); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(loc="lower center", frameon=False, fontsize=6.5, handlelength=1.2, borderaxespad=0.2)

    # (c), (d) maps
    for col, (vx, y, r, title) in enumerate(((vx0, y0, r0, "(c) no EOM, $\\delta/2\\pi=20.2$ MHz"),
                                             (vx1, y1, r1, "(d) EOM $\\beta=1.5$, $\\delta/2\\pi=20.0$ MHz"))):
        yf, rf = mirror(y, r)
        ax = fig.add_subplot(gs[0, col])
        mesh = ax.pcolormesh(edges(vx), edges(yf), rf.T, cmap=CMAP, norm=NORM, shading="flat",
                             edgecolors=SURFACE, linewidth=0.3)
        levels = np.log10(sorted(CONTOURS_REMAINING))
        cs = ax.contour(vx, yf, np.log10(rf.T), levels=levels, colors=SURFACE, linewidths=0.8, linestyles="solid")
        ax.clabel(cs, fmt={lv: eff_label(10 ** lv) for lv in levels}, fontsize=6.5, inline_spacing=2)
        a0, a1 = averaged(vx, y, r)
        ax.set_title(f"{title}\nmean over $v_x$ and $y_0$: {eff_label(a1)}", loc="left", fontsize=8.5)
        ax.set_xlabel("transverse velocity $v_x$ (m/s)")
        if col == 0:
            ax.set_ylabel("transverse position $y_0$ (mm)")
        else:
            ax.tick_params(labelleft=False)
        ax.set_xticks([-1, -0.5, 0, 0.5, 1]); ax.set_yticks([-15, -10, -5, 0, 5, 10, 15])
    cax = fig.add_subplot(gs[0, 2])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_ticks(TICKS_REMAINING); cb.set_ticklabels([eff_label(t) for t in TICKS_REMAINING]); cb.minorticks_off()
    cb.ax.invert_yaxis(); cb.set_label("cleanup efficiency $\\eta_\\mathrm{c}$")
    cb.outline.set_edgecolor(MUTED); cb.outline.set_linewidth(0.6)
    for c in CONTOURS_REMAINING:
        cb.ax.axhline(c, color=SURFACE, lw=0.8)

    note = [f"Averages: $v_x$ Gaussian, FWHM {VX_FWHM} m/s, truncated at $|v_x|\\leq1$ m/s; "
            "$y_0$ uniform in $|y_0|\\leq15$ mm.",
            "Markers in (a): circle best grid point, square setting of (d), triangle present setting (c); "
            "grey cells not computed.",
            "80 mW, 4.8 cm beam, 1.6 MHz EOM, 8 SPB-related initial levels; single trajectories."]
    fig.text(0.02, 0.005, "\n".join(note), ha="left", va="top", fontsize=7, color=INK)
    fig.savefig(OUT / "CU1_eom_doppler.pdf", bbox_inches="tight")
    fig.savefig(OUT / "CU1_eom_doppler.png", dpi=200, bbox_inches="tight")
    print("wrote CU1_eom_doppler.{pdf,png}")


if __name__ == "__main__":
    main()
