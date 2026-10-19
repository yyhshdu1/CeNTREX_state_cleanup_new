"""CU1 and CU2 cleanup efficiency versus transverse velocity vx and position y0 at the present set points.

    python scripts/plot_cleanup_maps_current.py   # docs/state_cleanup/figures/cleanup_maps_vx_y_current.{pdf,png}

Data (vx -1...1 m/s, y0 >= 0, mirrored to y0 < 0):
    CU1  no EOM, 3 cm beam, SPB-related start, 0.1 m/s x 1.5 mm   runs doppler_cu1_spb_3cm_det*MHz (cu1_3cm_scan.py)
    CU2  setting R1, 3 cm beam, uniform and SPB starts, 0.2 m/s x 3 mm   runs doppler_cu2_R1_*_3cm_coarse
Writes cleanup_map_vx_y_CU1_current and cleanup_map_vx_y_CU2_current (pdf, png) in docs/state_cleanup/figures.
Means over the detection acceptance |vx| <= 1 m/s, |y0| <= 8.5 mm (detection laser 1.7 cm high), weighted with
the 25 kV beam-simulation distributions: vx Gaussian sigma 0.65 m/s; y0 Gaussian sigma 10.95 mm at CU1
(lens + 0.5 m) and 10.35 mm at CU2 (lens + 1.0 m). log(remaining) is interpolated linearly in y0 between grid
points before the y0 integral.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import cu1_eom_doppler_scan as e1  # noqa: E402
from scripts import cu2_doppler_scan as c2  # noqa: E402
from scripts.plot_cu1_doppler_map import (CMAP, MUTED, NORM, OUT, SURFACE, TICKS_REMAINING, edges,  # noqa: E402
                                          eff_label, mirror)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

SVX, VCUT, YCUT = 0.65, 1.0, 8.5e-3
SY = {"CU1": 10.95e-3, "CU2": 10.35e-3}
ACCENT = "#eb6834"
CONTOURS = [0.1, 0.03, 0.01, 0.003, 0.001]          # 90, 97, 99, 99.7, 99.9 %


def acceptance_mean(g, vx, y, sy):
    wv = np.exp(-0.5 * (vx / SVX) ** 2) * (np.abs(vx) <= VCUT + 1e-12)
    t = np.ones_like(vx); t[[0, -1]] = 0.5; wv = wv * t; wv /= wv.sum()
    ym = np.linspace(0.0, YCUT, 171)
    rows = np.array([np.exp(np.interp(ym, y, np.log(r))) for r in g])
    wy = np.exp(-0.5 * (ym / sy) ** 2); t = np.ones_like(ym); t[[0, -1]] = 0.5; wy = wy * t; wy /= wy.sum()
    return wv @ rows @ wy


CU1_DETUNINGS = (19.0,)                              # MHz; maps doppler_cu1_spb_3cm_det*MHz (cu1_3cm_scan.py)


def cu1_run(d):
    return f"doppler_cu1_spb_3cm_det{d:g}MHz".replace(".", "p")


def grids():
    """key -> (vx, y0, remaining[vx, y0]) for every panel."""
    import h5py
    out = {}
    with h5py.File(ROOT / "results" / "SPB_J1_cleanup_scans.h5", "r") as h5:
        for d in CU1_DETUNINGS:
            r = h5["runs"][cu1_run(d)]
            out[f"cu1_{d:g}"] = (r["axes"]["vx_m_per_s"][...], r["axes"]["y0_mm"][...] * 1e-3, r["results"][...])
    vx, y = e1.VX_COARSE, e1.Y_COARSE
    for init in ("uniform", "spb"):
        dd = c2.load(c2.tag("R1", init))
        out[f"cu2_{init}"] = (vx, y, np.array([[dd[c2.key(v, yy)] for yy in y] for v in vx]))
    return out


def figure(panels, g, fname, figsize, note, text_ratio):
    """One figure: maps side by side with a shared colour bar, the conditions under each map, a note below."""
    import textwrap
    n = len(panels)
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, n + 1, width_ratios=[1] * n + [0.06 * n], height_ratios=[1, text_ratio],
                          wspace=0.13, hspace=0.30, left=0.10, right=0.88, top=0.94, bottom=0.20)
    for col, (key, stage, title, cond) in enumerate(panels):
        vx, y, r = g[key]
        i0 = int(np.argmin(abs(vx)))
        yf, rf = mirror(y, r)
        ax = fig.add_subplot(gs[0, col])
        mesh = ax.pcolormesh(edges(vx), edges(yf * 1e3), rf.T, cmap=CMAP, norm=NORM, shading="flat",
                             edgecolors=SURFACE, linewidth=0.3)
        levels = [np.log10(c) for c in sorted(CONTOURS) if rf.min() < c < rf.max()]
        if levels:
            cs = ax.contour(vx, yf * 1e3, np.log10(rf.T), levels=levels, colors=SURFACE, linewidths=0.7,
                            linestyles="solid")
            ax.clabel(cs, fmt={lv: eff_label(10 ** lv) for lv in levels}, fontsize=5.5, inline_spacing=1)
        for sgn in (1, -1):
            ax.axhline(sgn * YCUT * 1e3, color=ACCENT, lw=1.1, ls=(0, (4, 2)))
        ax.add_patch(Ellipse((0, 0), 2 * SVX, 2 * SY[stage] * 1e3, fill=False, ec=ACCENT, lw=0.9, ls=":"))
        ax.set_title(title, loc="left", fontsize=8)
        ax.set_xlabel("$v_x$ (m/s)")
        ax.set_xticks([-1, -0.5, 0, 0.5, 1]); ax.set_yticks([-15, -10, -5, 0, 5, 10, 15])
        if col == 0:
            ax.set_ylabel("$y_0$ (mm)")
        else:
            ax.tick_params(labelleft=False)
        mean = acceptance_mean(r, vx, y, SY[stage])
        tx = fig.add_subplot(gs[1, col]); tx.axis("off")
        lines = cond + [f"on axis: {eff_label(r[i0, 0])}",
                        rf"$\mathbf{{mean\ in\ acceptance:\ {100 * (1 - mean):.1f}\ \%}}$"]
        tx.text(0.0, 0.98, "\n".join(lines), ha="left", va="top", fontsize=6.3, linespacing=1.35,
                transform=tx.transAxes)
        print(f"{fname} {key}: on axis {eff_label(r[i0, 0])}, mean in acceptance {100 * (1 - mean):.2f} %")
    cax = fig.add_subplot(gs[0, n])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_ticks(TICKS_REMAINING); cb.set_ticklabels([eff_label(t) for t in TICKS_REMAINING]); cb.minorticks_off()
    cb.ax.invert_yaxis(); cb.ax.tick_params(labelsize=6.5)
    cb.set_label(r"cleanup efficiency $\eta_\mathrm{c}$", fontsize=7.5)
    cb.outline.set_edgecolor(MUTED); cb.outline.set_linewidth(0.6)
    width = int(figsize[0] * 26)                       # characters per line at 5.8 pt

    def wrap(par):                                     # never break inside $...$
        parts = par.split("$")
        keep = "$".join(p.replace(" ", "\x00") if k % 2 else p for k, p in enumerate(parts))
        return textwrap.fill(keep, width).replace("\x00", " ")
    wrapped = "\n".join(wrap(par) for par in note)
    fig.text(0.10, 0.19, wrapped, ha="left", va="top", fontsize=5.8, linespacing=1.4)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 220})):
        fig.savefig(OUT / f"{fname}.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)


def main():
    g = grids()
    d0, faom = c2.REGIONS["R1"]
    d1 = d0 + faom - 211.829

    def common(grid):
        return (r"13 passes, 10 % loss per pass, X/Z polarization switching at 0.5 MHz, Doppler shift $\pm k v_x$ "
                rf"alternating pass by pass; single trajectories at $v_z=184$ m/s; grid {grid}, computed "
                r"for $y_0\geq0$ and mirrored.")

    def note(stage, plane, grid):
        return [rf"Dashed orange: detection acceptance $|y_0|\leq8.5$ mm (detection laser 1.7 cm high); the grid "
                rf"spans the accepted $|v_x|\leq1$ m/s. Dotted orange: $1\sigma$ of the beam distribution at "
                rf"{stage} ({plane}, lens at 25 kV): $\sigma_{{v_x}}=0.65$ m/s, $\sigma_y={SY[stage] * 1e3:.2f}$ mm; "
                rf"it weights the mean in acceptance. White lines: contours of $\eta_\mathrm{{c}}$. On axis: "
                rf"$v_x=0$, $y_0=0$.", common(grid)]

    labels = "abcd"
    cu1_panels = [(f"cu1_{d:g}", "CU1", (f"({labels[k]}) " if len(CU1_DETUNINGS) > 1 else "") + f"CU1, +{d:g} MHz",
                   [r"R(1) $F_1'=5/2$, $F'=2$; 150 V/cm",
                    rf"$\delta/2\pi=+{d:.1f}$ MHz, no EOM, 80 mW",
                    r"beam 3 cm high ($\sigma_y=7.5$ mm)",
                    r"start: 8 SPB-related $J=1$ levels"]) for k, d in enumerate(CU1_DETUNINGS)]
    figure(cu1_panels, g, "cleanup_map_vx_y_CU1_current", (3.6 if len(cu1_panels) == 1 else 5.4, 4.9),
           note("CU1", "lens + 0.5 m", r"0.1 m/s $\times$ 1.5 mm"), 0.36)

    common_cu2 = [r"R(2) $F_1'=7/2$, $F'=3$ + $F'=4$; 200 V/cm",
                  rf"setting R1: $f_\mathrm{{AOM}}={faom:.1f}$ MHz,",
                  rf"$\delta_0/2\pi=+{d0:.1f}$, $\delta_1/2\pi={d1:+.2f}$ MHz",
                  r"80 + 24 mW; EOM 17.5 MHz, $\beta=1.5$",
                  r"beam 3 cm high ($\sigma_y=7.5$ mm)"]
    figure([("cu2_uniform", "CU2", "(a) CU2, all $J=2$ sublevels", common_cu2 + ["start: all 20 $J=2$ sublevels"]),
            ("cu2_spb", "CU2", "(b) CU2, SPB source levels",
             common_cu2 + ["start: 4 $J=2$, $m_J=0$ levels", "(lens weights 2/3, 1/9, 1/9, 1/9)"])],
           g, "cleanup_map_vx_y_CU2_current", (5.4, 5.4), note("CU2", "lens + 1.0 m", r"0.2 m/s $\times$ 3 mm"), 0.66)


if __name__ == "__main__":
    main()
