"""CU1 with the laser tilted by an angle theta from perpendicular to the molecular beam (in the plane of v_z).

    python scripts/cu1_tilt_scan.py compute [--power P --height H]   # extend the coarse cube to |vx_eff| <= 2 m/s
    python scripts/cu1_tilt_scan.py archive [--power P --height H]   # store the extended cube in the HDF5 file
    python scripts/cu1_tilt_scan.py analyze [--power P --height H] [--parallel]   # table versus tilt
    python scripts/cu1_tilt_scan.py plot    [--power P --height H] [--parallel]   # figure versus tilt

Defaults: 80 mW, 3 cm (the present beam); the base cube is detuning_vx_y_cu1_spb_3cm_coarse, other conditions
use detuning_vx_y_cu1_spb_<H>cm_<P>mW_coarse from scripts/cu1_conditions_scan.py. The extension to 1 < |vx| <= 2 m/s
costs 520 trajectories. Figure: docs/state_cleanup/figures/CU1_laser_tilt[_<H>cm_<P>mW].{pdf,png}.

A tilt adds k v_z sin(theta) to the Doppler shift of every pass: v_z/lambda = 0.68 MHz per mrad at 184 m/s.
Two cases:
  antiparallel  retro-reflecting two-prism multipass, as modelled: the shift flips sign with the transverse one,
                k.v = +-k (v_x + v_z tan theta), i.e. every molecule has an effective transverse velocity
                v_x + v_z theta in CU1, while detection still selects on the true v_x.
  parallel      (--parallel) the z component of k has the same sign in every pass (e.g. a plane-mirror zig-zag):
                a detuning offset v_z theta / lambda, removable by re-tuning the laser.
Model: build_cu1, no EOM, 8 SPB-related levels, same grids and acceptance weights as cu1_3cm_scan.py. The
forward-velocity spread (spread of v_z theta) and the change of the pass positions are neglected.
"""
import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cu1_3cm_scan import DETUNINGS, SY_BEAM, VX_C, Y_C, key  # noqa: E402
from scripts.cu1_conditions_scan import cond_tag, run, setup  # noqa: E402
from scripts.plot_cleanup_maps_current import acceptance_mean  # noqa: E402

HDF = ROOT / "results" / "SPB_J1_cleanup_scans.h5"
TMP = ROOT / "tmp" / "cu1_tilt"
VZ, LAMBDA = 184.0, 271.7e-9
VX_EXT = np.round(np.concatenate([np.arange(-2.0, -1.01, 0.2), np.arange(1.2, 2.01, 0.2)]), 6)
VX_ALL = np.round(np.arange(-2.0, 2.001, 0.2), 6)
MHZ_PER_MRAD = VZ / LAMBDA * 1e-3 / 1e6
I19 = int(np.argmin(abs(DETUNINGS - 19.0)))


def names(power_mw, height_cm):
    """Base run, extension checkpoint folder, archived run and figure name for one condition."""
    if (power_mw, height_cm) == (80.0, 3.0):            # the present beam keeps its original names
        base, ext, fig = "detuning_vx_y_cu1_spb_3cm_coarse", TMP, "CU1_laser_tilt"
    else:
        ct = cond_tag(power_mw, height_cm)
        base, ext, fig = f"detuning_vx_y_cu1_spb_{ct}_coarse", TMP / ct, f"CU1_laser_tilt_{ct}"
    return dict(base=base, ext=ext, archived=base + "_vx2", fig=fig)


def compute(power_mw, height_cm):
    prepared, gidx, rho0 = setup(height_cm)
    pts = [(d, v, y) for d in DETUNINGS for v in VX_EXT for y in Y_C]
    run(names(power_mw, height_cm)["ext"], "ext", pts, prepared, gidx, rho0, power_mw)


def cube(power_mw, height_cm):
    """Remaining J=1 population [detuning, vx_eff in VX_ALL, y0 in Y_C]."""
    import h5py
    from scripts.cu1_conditions_scan import load
    n = names(power_mw, height_cm)
    with h5py.File(HDF, "r") as h5:
        g = h5["runs"][n["base"]]
        base, vx_b = g["results"][...], np.round(g["axes"]["vx_m_per_s"][...], 6)
    ext = load(n["ext"], "ext")
    c = np.empty((len(DETUNINGS), len(VX_ALL), len(Y_C)))
    for j, v in enumerate(VX_ALL):
        if v in vx_b:
            c[:, j] = base[:, int(np.argmin(abs(vx_b - v)))]
        else:
            c[:, j] = [[ext[key((d, v, y))] for y in Y_C] for d in DETUNINGS]
    return c


def shifted(c_d, dv):
    """Grid on the true vx (VX_C) x y0 for an effective-velocity shift dv; log-linear interpolation in vx."""
    return np.array([[np.exp(np.interp(v + dv, VX_ALL, np.log(c_d[:, k]))) for k in range(len(Y_C))] for v in VX_C])


def curves(thetas, c, parallel=False):
    """Mean in acceptance (remaining population) versus tilt; edge flag when the best detuning is at the grid edge."""
    def mean(g):
        return acceptance_mean(g, VX_C, Y_C, SY_BEAM)
    out = {k: [] for k in ("alt19", "altbest", "detbest", "edge", "common")}
    det_curve = None
    if parallel:
        cols = [list(VX_ALL).index(v) for v in VX_C]
        det_curve = np.array([mean(c[i][cols]) for i in range(len(DETUNINGS))])
    for th in thetas:
        dv = th * VZ * 1e-3
        alt = np.array([mean(shifted(c[i], dv)) for i in range(len(DETUNINGS))])
        ib = int(np.argmin(alt))
        out["alt19"].append(alt[I19]); out["altbest"].append(alt[ib]); out["detbest"].append(DETUNINGS[ib])
        out["edge"].append(ib in (0, len(DETUNINGS) - 1))
        if parallel:
            f = 19.0 + th * MHZ_PER_MRAD
            inside = DETUNINGS[0] <= f <= DETUNINGS[-1]
            out["common"].append(np.exp(np.interp(f, DETUNINGS, np.log(det_curve))) if inside else np.nan)
    return {k: np.asarray(v) for k, v in out.items() if len(v)}


def analyze(power_mw, height_cm, parallel=False):
    thetas = np.array([-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5], dtype=float)
    r = curves(thetas, cube(power_mw, height_cm), parallel)
    print(f"CU1 {power_mw:g} mW, {height_cm:g} cm; {MHZ_PER_MRAD:.3f} MHz and {VZ * 1e-3:.3f} m/s per mrad")
    print("theta (mrad) | shift (MHz) | dvx (m/s) | antiparallel, +19 MHz | antiparallel, best detuning"
          + (" | parallel, laser at +19 MHz" if parallel else ""))
    for k, th in enumerate(thetas):
        line = (f"   {th:+5.1f}  {th * MHZ_PER_MRAD:+6.2f}  {th * VZ * 1e-3:+5.2f}   {100 * (1 - r['alt19'][k]):7.2f} %   "
                f"{'>=' if r['edge'][k] else '  '}{100 * (1 - r['altbest'][k]):7.2f} % at {r['detbest'][k]:+.1f} MHz"
                f"{' (grid edge)' if r['edge'][k] else ''}")
        if parallel:
            line += ("   outside scanned detunings" if np.isnan(r["common"][k]) else
                     f"   {100 * (1 - r['common'][k]):7.2f} %")
        print(line)


def plot(power_mw, height_cm, parallel=False):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from scripts.plot_cu1_doppler_map import GRID, INK, MUTED, OUT, SERIES   # also sets the report rcParams
    n = names(power_mw, height_cm)
    th_line = np.round(np.arange(-5.0, 5.001, 0.25), 3)
    th_mark = np.arange(-5.0, 5.001, 1.0)
    c = cube(power_mw, height_cm)
    rl, rm = curves(th_line, c, parallel), curves(th_mark, c, parallel)
    eff = {k: 100 * (1 - rl[k]) for k in ("alt19", "altbest") + (("common",) if parallel else ())}
    effm = {k: 100 * (1 - rm[k]) for k in eff}
    e0 = eff["alt19"][int(np.argmin(abs(th_line)))]
    ref = None
    if (power_mw, height_cm) != (80.0, 3.0):            # the present beam as a reference curve
        ref = 100 * (1 - curves(th_line, cube(80.0, 3.0))["alt19"])

    fig, ax = plt.subplots(figsize=(6.3, 4.1))
    lo = np.nanmin([np.nanmin(v) for v in eff.values()] + ([ref.min()] if ref is not None else []))
    hi = max(np.nanmax([np.nanmax(v) for v in eff.values()]), ref.max() if ref is not None else 0)
    y0, y1 = np.floor(lo) - 0.6, hi + 1.0
    ax.axvspan(-1, 1, color="#ecebe6", lw=0, zorder=0)
    ax.text(0, y1 - 0.1, r"$|\theta|\leq1$ mrad", ha="center", va="top", fontsize=7.5, color=MUTED)
    ax.axhline(e0, color=MUTED, lw=1.0, ls=(0, (1.5, 2)), zorder=1)
    specs = [("alt19", SERIES[0], "o", "antiparallel passes, laser at +19.0 MHz", "+19.0 MHz"),
             ("altbest", SERIES[1], "s", "antiparallel passes, detuning re-optimized", "re-tuned")]
    if parallel:
        specs.append(("common", SERIES[2], "^", "parallel passes, laser kept at +19.0 MHz", "parallel, +19.0 MHz"))
    end_labels = []                                     # (y, text, colour) at the right end, spread out below
    if ref is not None:
        ax.plot(th_line, ref, color="#9c9b96", lw=1.3, ls="-", zorder=2)
        end_labels.append((ref[-1], "80 mW, 3 cm, +19.0 MHz", MUTED))
    for k, col, mk, lab, short in specs:
        ax.plot(th_line, eff[k], color=col, lw=1.8, zorder=3)
        edge = rm["edge"] if k == "altbest" else np.zeros(len(th_mark), bool)
        ax.plot(th_mark[~edge], effm[k][~edge], ls="none", marker=mk, ms=5.5, color=col, mec="white", mew=0.8,
                zorder=4)
        ax.plot(th_mark[edge], effm[k][edge], ls="none", marker=mk, ms=5.5, mfc="white", mec=col, mew=1.3, zorder=4)
        end_labels.append((eff[k][-1], short, INK))
    gap = 0.05 * (y1 - y0)
    placed = []
    for y, txt, col in sorted(end_labels):
        y = max(y, placed[-1] + gap) if placed else y
        placed.append(y)
        ax.text(5.18, y, txt, ha="left", va="center", fontsize=7.2, color=col)
    ax.text(5.18, e0 + 0.15, "no tilt", ha="left", va="bottom", fontsize=7.2, color=MUTED)
    ax.set_xlim(-5.3, 5.3)
    ax.set_ylim(y0, y1)
    ax.set_xlabel(r"laser tilt from perpendicular, $\theta$ (mrad)")
    ax.set_ylabel(r"mean cleanup efficiency in acceptance (%)")
    ax.set_xticks(np.arange(-5, 6, 1))
    ax.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    top = ax.secondary_xaxis("top", functions=(lambda t: t * MHZ_PER_MRAD, lambda f: f / MHZ_PER_MRAD))
    top.set_xlabel(r"Doppler shift per pass, $v_z\theta/\lambda$ (MHz)", labelpad=4)
    top.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    handles = [Line2D([], [], color=col, lw=1.8, marker=mk, ms=5.5, mec="white", mew=0.8,
                      label=lab + (f", {power_mw:g} mW, {height_cm:g} cm" if ref is not None else ""))
               for _, col, mk, lab, _ in specs]
    if ref is not None:
        handles.append(Line2D([], [], color="#9c9b96", lw=1.3, label="same, 80 mW, 3 cm beam (present), +19.0 MHz"))
    handles += [Line2D([], [], color=MUTED, lw=1.0, ls=(0, (1.5, 2)), label=f"no tilt ({e0:.2f} %)"),
                Line2D([], [], ls="none", marker="s", ms=5.5, mfc="white", mec=SERIES[1], mew=1.3,
                       label="optimum at the edge of the +17 to +23 MHz scan (lower bound)")]
    leg = ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize=6.8,
                    frameon=False, columnspacing=1.5, handlelength=2.2)
    note = (rf"CU1, {power_mw:g} mW, {height_cm:g} cm beam ($\sigma_y={height_cm * 10 / 4:g}$ mm), no EOM, start in the 8 "
            r"SPB-related $J=1$ levels; mean over $v_x$ (Gaussian, $\sigma=0.65$ m/s, $|v_x|\leq1$ m/s) and $y_0$ "
            rf"($\sigma=10.95$ mm, $|y_0|\leq8.5$ mm). A tilt adds $k v_z\sin\theta$ to every pass, "
            rf"{MHZ_PER_MRAD:.2f} MHz/mrad at $v_z={VZ:.0f}$ m/s. Antiparallel passes (retro-reflecting prisms, as in "
            r"the model): the shift alternates in sign like the transverse one, i.e. an effective $v_x$ shift of "
            rf"$v_z\theta$ = {VZ * 1e-3:.3f} m/s per mrad that re-tuning cannot remove. Forward-velocity spread "
            r"neglected.")
    if parallel:
        note += (" Parallel passes: a common detuning offset; with the laser re-tuned there is no loss; missing points "
                 "lie outside the scanned detunings.")
    fig.tight_layout()
    fig.canvas.draw()                                   # place the note just under the legend
    y_leg = leg.get_window_extent().transformed(fig.transFigure.inverted()).y0
    fig.text(0.02, y_leg - 0.03, textwrap.fill(note, 150), ha="left", va="top", fontsize=6.6, color=INK,
             linespacing=1.4)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 220})):
        fig.savefig(OUT / f"{n['fig']}.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)
    print(f"wrote {n['fig']}.{{pdf,png}}")


def archive(power_mw, height_cm):
    """Store the extended cube (vx_eff to +-2 m/s) in results/SPB_J1_cleanup_scans.h5."""
    from scripts.cleanup_spb_doppler import CU1_SOLVER
    from scripts.scan_hdf import save_scan_hdf
    n = names(power_mw, height_cm)
    save_scan_hdf(HDF, n["archived"], cube(power_mw, height_cm),
                  axes={"laser_detuning_MHz": (DETUNINGS, "MHz"), "vx_m_per_s": (VX_ALL, "m/s"),
                        "y0_mm": (Y_C * 1e3, "mm")},
                  quantity="final_population_remaining_in_X_J1",
                  fixed_parameters={"laser_power_mW": power_mw, "beam_height_4sigma_cm": height_cm, "eom": "none",
                                    "initial_state": "spb (8 levels, equal weights)",
                                    "doppler": "k*vx, sign alternating pass by pass",
                                    "vx_range": f"|vx| <= 1 m/s from {n['base']}, 1 < |vx| <= 2 m/s added by "
                                                "cu1_tilt_scan.py (centrex-tlf 0.3.0)"},
                  model={"script": "scripts/cu1_tilt_scan.py", "build": "cleanup_spb_doppler.build_cu1"},
                  solver=dict(CU1_SOLVER), notes="CU1 coarse cube extended to |vx| <= 2 m/s, for the laser-tilt study "
                  "(alternating tilt Doppler = effective vx shift v_z theta).", overwrite=False)
    print(f"saved {n['archived']}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["compute", "archive", "analyze", "plot"])
    ap.add_argument("--power", type=float, default=80.0, help="laser power, mW")
    ap.add_argument("--height", type=float, default=3.0, help="full 4-sigma beam height, cm")
    ap.add_argument("--parallel", action="store_true", help="also the parallel-pass (common shift) case")
    a = ap.parse_args()
    if a.stage == "compute":
        compute(a.power, a.height)
    elif a.stage == "archive":
        archive(a.power, a.height)
    elif a.stage == "analyze":
        analyze(a.power, a.height, a.parallel)
    else:
        plot(a.power, a.height, a.parallel)
