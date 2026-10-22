"""Alignment test: CU1 with the first prism removed (two antiparallel passes), scan the laser and look for a
double peak in the F=2 depletion seen by the R(1) F'=3 probe.

    python scripts/cu1_two_pass_test.py compute     # grid detuning x vx_eff x power at y0 = 0 (~7500 trajectories)
    python scripts/cu1_two_pass_test.py analyze     # spectra versus tilt and vx spread; dip visibility; figure

The two passes see the laser at delta -+ k (v_x + v_z theta): a tilt theta splits the depletion line into two
peaks 2 v_z theta / lambda = 1.35 MHz/mrad apart, each broadened by the transverse Doppler spread k sigma_vx.
Model: build_cu1(n_pass=2) (passes 1 and 2, the second 10 % weaker), 150 V/cm, 2 cm beam (sigma_y = 5 mm),
X/Z switching, no EOM, molecule on the beam axis (y0 = 0). Start: the AP1 triplet levels F1=3/2, F=2, mF=0,+-1
(RC on P(2) F'=2, so no singlet). Observable: F=2 population after CU1 (what R(1) F'=3 counts), relative to
CU1 off. The spectrum for tilt theta and a Gaussian v_x distribution (sigma, |v_x| <= 1 m/s) is the weighted
mean of the grid at v_x + v_z theta.
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cleanup_spb_doppler import CU1_SOLVER, build_cu1, describe  # noqa: E402
from scripts.cu2_laser_aom_scan import FIRST_PASS, PASS_SPACING, SIGMA_Z, VZ, log  # noqa: E402

OUT = ROOT / "tmp" / "cu1_two_pass"
DETS = np.round(np.arange(11.0, 26.001, 0.25), 3)            # MHz
VXE = np.round(np.arange(-2.0, 2.001, 0.1), 6)              # m/s, effective v_x + v_z theta
POWERS = (0.3, 1.0, 3.0)                                    # mW, laser power into the two-pass setup
SIGMA_Y = 5e-3                                              # m, 2 cm beam
T_END = (FIRST_PASS + PASS_SPACING + 8 * SIGMA_Z) / VZ      # light is negligible after the second pass
LAMBDA = 271.7e-9
MHZ_PER_MRAD = VZ / LAMBDA * 1e-9                           # per-pass shift, MHz per mrad
SLOTS = ["laser_detuning", "vx", "y0", "laser_power"]
CHUNK, THREADS = 400, 20


def setup():
    prepared, gidx, obe = build_cu1(return_obe=True, doppler=True, sigma_y_m=SIGMA_Y, n_pass=2)
    lab = [(describe(obe.QN[i])["F1"], describe(obe.QN[i])["F"], round(describe(obe.QN[i])["mF"])) for i in gidx]
    trip = [gidx[lab.index((1.5, 2, m))] for m in (-1, 0, 1)]
    f2 = [gidx[k] for k, l in enumerate(lab) if l[:2] == (1.5, 2)]
    rho0 = np.zeros((len(obe.QN),) * 2, dtype=np.complex128)
    rho0[trip, trip] = 1.0 / 3
    return prepared, rho0, f2


def compute():
    from centrex_tlf import lindblad
    OUT.mkdir(parents=True, exist_ok=True)
    prepared, rho0, f2 = setup()
    for p in POWERS:
        f = OUT / f"grid_{p:g}mW.npy"
        if f.exists():
            continue
        pts = np.array([(d, v) for d in DETS for v in VXE])
        vals = []
        for s in range(0, len(pts), CHUNK):
            c = pts[s:s + CHUNK]
            rows = np.column_stack([2 * np.pi * 1e6 * c[:, 0], c[:, 1], np.zeros(len(c)), np.full(len(c), p * 1e-3)])
            t0 = time.perf_counter()
            res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=SLOTS,
                                          parameter_batch=rows.astype(np.complex128), threads=THREADS, **CU1_SOLVER)
            vals.append(np.asarray(res.values).real[:, f2].sum(axis=1))
            log(f"{p:g} mW: {s + len(c)}/{len(pts)} ({time.perf_counter() - t0:.0f} s)")
        np.save(f, np.concatenate(vals).reshape(len(DETS), len(VXE)))


def spectrum(grid, theta_mrad, sigma_vx):
    """Mean F=2 depletion versus detuning for a tilt and a Gaussian v_x distribution cut at |v_x| <= 1 m/s."""
    vx = np.round(np.arange(-1.0, 1.0001, 0.05), 6)
    w = np.exp(-0.5 * (vx / sigma_vx) ** 2)
    w /= w.sum()
    shift = theta_mrad * VZ * 1e-3
    remain = np.array([np.interp(vx + shift, VXE, grid[i]) for i in range(len(DETS))])
    return 1 - remain @ w


def dip(spec):
    """(two maxima?, dip depth relative to the smaller peak, splitting MHz) on a 0.25 MHz grid."""
    imax = [i for i in range(1, len(spec) - 1) if spec[i] >= spec[i - 1] and spec[i] > spec[i + 1]]
    imax = sorted(imax, key=lambda i: -spec[i])[:2]
    if len(imax) < 2:
        return False, 0.0, 0.0
    a, b = sorted(imax)
    lo = spec[a:b + 1].min()
    return True, 1 - lo / min(spec[a], spec[b]), DETS[b] - DETS[a]


def fwhm(spec):
    half = spec.max() / 2
    above = np.where(spec >= half)[0]
    return DETS[above[-1]] - DETS[above[0]]


def analyze():
    grids = {p: np.load(OUT / f"grid_{p:g}mW.npy") for p in POWERS}
    thetas = (0, 1, 2, 3, 4, 5)
    for sig in (0.65, 0.35, 0.2):
        print(f"\nsigma_vx = {sig} m/s (Doppler FWHM per pass {2.3548 * sig / LAMBDA * 1e-6:.2f} MHz)")
        print("  power | theta (mrad): peak depletion %, FWHM MHz, dip (relative depth) ...")
        for p, g in grids.items():
            cells = []
            for th in thetas:
                s = spectrum(g, th, sig)
                two, depth, split = dip(s)
                cells.append(f"{th}: {100 * s.max():4.1f}%/{fwhm(s):4.1f}/" + (f"dip {100 * depth:4.1f}%" if two else "  one "))
            print(f"  {p:4.1f} mW | " + " | ".join(cells))
    plot(grids)


def plot(grids):
    import matplotlib.pyplot as plt
    from scripts.plot_cu1_doppler_map import GRID, INK, MUTED, OUT as FIG, SERIES
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.0), sharey=False)
    p = 1.0
    for ax, sig in zip(axes, (0.65, 0.2)):
        for th, col in zip((0, 2, 4, 5), ("#9c9b96", SERIES[0], SERIES[1], SERIES[2])):
            s = spectrum(grids[p], th, sig)
            ax.plot(DETS, 100 * s, color=col, lw=1.6, label=rf"$\theta={th}$ mrad")
        ax.set_title(rf"(${'a' if sig > 0.5 else 'b'}$) $\sigma_{{v_x}}={sig}$ m/s", loc="left", fontsize=8.5)
        ax.set_xlabel(r"CU1 laser detuning $\delta/2\pi$ (MHz)")
        ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel(r"$F=2$ depletion seen by R(1) $F'=3$ (%)")
    axes[1].legend(fontsize=7, frameon=False, loc="upper right")
    note = (rf"Two antiparallel passes (first prism removed), {p:g} mW, 2 cm beam, 150 V/cm, X/Z switching, molecule on "
            r"the beam axis; start in the AP1 triplet levels. A tilt $\theta$ splits each line into two peaks "
            rf"$2v_z\theta/\lambda$ = {2 * MHZ_PER_MRAD:.2f} MHz/mrad apart; each is broadened by the transverse Doppler "
            r"spread. (a) the 25 kV lens beam ($\sigma_{v_x}=0.65$ m/s), (b) a collimated beam ($\sigma_{v_x}=0.2$ m/s).")
    import textwrap
    fig.text(0.02, -0.04, textwrap.fill(note, 135), ha="left", va="top", fontsize=6.6, color=INK, linespacing=1.4)
    fig.tight_layout()
    for ext, kw in (("pdf", {}), ("png", {"dpi": 220})):
        fig.savefig(FIG / f"CU1_two_pass_tilt_test.{ext}", bbox_inches="tight", **kw)
    plt.close(fig)
    print("wrote CU1_two_pass_tilt_test.{pdf,png}")


if __name__ == "__main__":
    {"compute": compute, "analyze": analyze}[sys.argv[1]]()
