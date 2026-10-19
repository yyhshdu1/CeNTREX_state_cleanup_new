"""CU1 with a 1.6 MHz phase EOM to broaden the laser against the transverse Doppler shift.

    python scripts/cu1_eom_doppler_scan.py test              # beta = 0 reproduces the unmodulated model; timing
    python scripts/cu1_eom_doppler_scan.py coarse            # beta x detuning x vx at y0 = 0
    python scripts/cu1_eom_doppler_scan.py yavg B D [B D...] # vx x y0 at a few (beta, detuning) points
    python scripts/cu1_eom_doppler_scan.py map B D           # full 21 x 11 vx x y0 map at one point
    python scripts/cu1_eom_doppler_scan.py summary           # tables from the checkpoints

Model: scripts/cleanup_spb_doppler.build_cu1 (R(1) F1'=5/2 F'=2, 150 V/cm, 80 mW, 4.8 cm beam, 13 passes,
10 % loss per pass, X/Z switching at 0.5 MHz), drive multiplied by exp(i beta sin(2 pi 1.6 MHz t)).
Initial state: the 8 SPB-related J=1 levels (|m_J|=1 content > 0.5), equal weights.
Transverse distribution: vx Gaussian with FWHM 1.8 m/s, truncated at |vx| <= 1 m/s (molecules outside
do not reach the detection region); y0 uniform in |y0| <= 15 mm. Both are weights applied afterwards, so
the stored grids can be re-averaged with any other distribution.
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cleanup_spb_doppler import CU1_DETUNING, CU1_SOLVER, build_cu1, describe  # noqa: E402
from scripts.cu2_laser_aom_scan import T_END, log  # noqa: E402
from scripts.scan_hdf import save_scan_hdf  # noqa: E402

OUT = ROOT / "tmp" / "cu1_eom_doppler"
F_EOM_MHZ = 1.6
VX_FWHM = 1.8                                           # m/s
VX_SIGMA = VX_FWHM / (2 * np.sqrt(2 * np.log(2)))
VX_CUT = 1.0                                            # m/s
BETAS = np.round(np.arange(0.0, 4.01, 0.5), 3)
DETUNINGS = np.array([19.0, 19.5, 20.0, 20.2, 20.5, 21.0, 21.5])     # MHz
VX_COARSE = np.round(np.linspace(-1.0, 1.0, 11), 6)                  # m/s
VX_FINE = np.round(np.linspace(-1.0, 1.0, 21), 6)
Y_COARSE = np.round(np.linspace(0.0, 15e-3, 6), 9)                   # m
Y_FINE = np.round(np.linspace(0.0, 15e-3, 11), 9)
CHUNK, THREADS = 40, 8
SLOTS = ["laser_detuning", "eom_beta", "vx", "y0"]


def vx_weights(vx):
    """Trapezoid weights times the truncated Gaussian, normalised."""
    w = np.exp(-0.5 * (vx / VX_SIGMA) ** 2) * (np.abs(vx) <= VX_CUT + 1e-9)
    trap = np.ones_like(vx); trap[[0, -1]] = 0.5
    w = w * trap
    return w / w.sum()


def y_weights(y):
    """Uniform in |y0| <= 15 mm, from the y0 >= 0 half (trapezoid)."""
    w = np.ones_like(y); w[[0, -1]] = 0.5
    return w / w.sum()


def setup():
    from centrex_tlf import lindblad  # noqa: F401
    prepared, gidx, obe = build_cu1(return_obe=True, doppler=True)
    spb = [i for i in gidx if sum(v for m, v in describe(obe.QN[i])["p_mj"].items() if abs(m) == 1) > 0.5]
    rho0 = np.zeros((len(obe.QN),) * 2, dtype=np.complex128)
    rho0[spb, spb] = 1.0 / len(spb)
    return prepared, gidx, rho0


def key(p):
    return tuple(round(float(v), 9) for v in p)


def load(tag):
    done = {}
    for f in sorted(OUT.glob(f"{tag}_*.npz")):
        z = np.load(f)
        done.update({key(p): r for p, r in zip(z["points"], z["remaining"])})
    return done


def run(tag, points, prepared, gidx, rho0):
    """points: rows (detuning MHz, beta, vx m/s, y0 m). Skips points already checkpointed."""
    from centrex_tlf import lindblad
    OUT.mkdir(parents=True, exist_ok=True)
    done = load(tag)
    todo = np.asarray([p for p in points if key(p) not in done]).reshape(-1, 4)
    log(f"{tag}: {len(points)} points, {len(todo)} to compute")
    for c, s in enumerate(range(0, len(todo), CHUNK)):
        chunk = todo[s:s + CHUNK]
        batch = np.column_stack([2 * np.pi * 1e6 * chunk[:, 0], chunk[:, 1], chunk[:, 2], chunk[:, 3]])
        t0 = time.perf_counter()
        res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=SLOTS,
                                      parameter_batch=batch.astype(np.complex128), threads=THREADS, **CU1_SOLVER)
        rem = np.asarray(res.values).real[:, gidx].sum(axis=1)
        np.savez(OUT / f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}_{c:03d}.npz", points=chunk, remaining=rem)
        done.update({key(p): r for p, r in zip(chunk, rem)})
        log(f"{tag}: chunk {c} ({len(chunk)} pts) {time.perf_counter() - t0:.0f} s; "
            f"{sum(key(p) in done for p in points)}/{len(points)} done")
    return done


def grid(done, d, b, vxs, ys):
    return np.array([[done[key((d, b, v, y))] for y in ys] for v in vxs])


def eff(r):
    return f"{100 * (1 - r):.3f} %"


def stage_test():
    t0 = time.perf_counter()
    prepared, gidx, rho0 = setup()
    log(f"model ready ({time.perf_counter() - t0:.0f} s)")
    pts = [(CU1_DETUNING, b, 0.0, 0.0) for b in (0.0, 1.0, 2.0, 3.0)] + \
          [(CU1_DETUNING, b, 0.5, 0.0) for b in (0.0, 1.0, 2.0, 3.0)]
    done = run("test", pts, prepared, gidx, rho0)
    for p in pts:
        log(f"   detuning {p[0]} MHz, beta {p[1]}, vx {p[2]} m/s: cleanup {eff(done[key(p)])}")
    log("reference (unmodulated map, vx = 0 / 0.5 m/s): 99.965 % / 97.65 %")


def stage_coarse():
    prepared, gidx, rho0 = setup()
    pts = [(d, b, v, 0.0) for b in BETAS for d in DETUNINGS for v in VX_COARSE]
    run("coarse", pts, prepared, gidx, rho0)
    summary()


def stage_yavg(pairs):
    prepared, gidx, rho0 = setup()
    pts = [(d, b, v, y) for b, d in pairs for v in VX_COARSE for y in Y_COARSE]
    done = run("yavg", pts, prepared, gidx, rho0)
    wv, wy = vx_weights(VX_COARSE), y_weights(Y_COARSE)
    for b, d in pairs:
        g = grid(done, d, b, VX_COARSE, Y_COARSE)
        log(f"beta {b}, detuning {d} MHz: vx-averaged at y0=0 {eff(wv @ g[:, 0])}; "
            f"vx- and y-averaged {eff(wv @ g @ wy)}; on axis vx=0 {eff(g[5, 0])}")


def stage_map(b, d):
    prepared, gidx, rho0 = setup()
    pts = [(d, b, v, y) for v in VX_FINE for y in Y_FINE]
    done = run("map", pts, prepared, gidx, rho0)
    g = grid(done, d, b, VX_FINE, Y_FINE)
    wv, wy = vx_weights(VX_FINE), y_weights(Y_FINE)
    name = f"doppler_cu1_spb_eom1p6MHz_beta{b:g}_det{d:g}MHz".replace(".", "p")
    save_scan_hdf(ROOT / "results" / "SPB_J1_cleanup_scans.h5", name, g,
                  axes={"vx_m_per_s": (VX_FINE, "m/s"), "y0_mm": (Y_FINE * 1e3, "mm")},
                  quantity="final_population_remaining_in_X_J1",
                  fixed_parameters={"laser_detuning_MHz": d, "eom_frequency_MHz": F_EOM_MHZ, "eom_beta": b,
                                    "laser_power_mW": 80.0, "initial_state": "spb (8 levels, equal weights)",
                                    "doppler": "k*vx, sign alternating pass by pass (two-prism multipass)",
                                    "y_symmetry": "results for y0 >= 0; symmetric in y0"},
                  model={"script": "scripts/cu1_eom_doppler_scan.py", "build": "cleanup_spb_doppler.build_cu1"},
                  solver=dict(CU1_SOLVER), notes="CU1 with a 1.6 MHz broadening EOM.", overwrite=False)
    log(f"saved {name}; vx-averaged (FWHM {VX_FWHM} m/s, |vx|<={VX_CUT}) at y0=0 {eff(wv @ g[:, 0])}, "
        f"vx- and y-averaged {eff(wv @ g @ wy)}, on axis {eff(g[10, 0])}")


def stage_extend():
    """Lower detunings for small beta (the vx-averaged optimum moved to the edge) and beta 1.25, 1.75."""
    prepared, gidx, rho0 = setup()
    pts = [(d, b, v, 0.0) for b in (0.0, 0.5, 1.0, 1.5, 2.0) for d in (16.0, 16.5, 17.0, 17.5, 18.0, 18.5)
           for v in VX_COARSE]
    pts += [(d, b, v, 0.0) for b in (1.25, 1.75) for d in (18.0, 18.5, 19.0, 19.5, 20.0, 20.5) for v in VX_COARSE]
    run("coarse", pts, prepared, gidx, rho0)
    summary()


def summary():
    done = load("coarse")
    wv = vx_weights(VX_COARSE)
    BETAS = sorted({k[1] for k in done}); DETUNINGS = sorted({k[0] for k in done})
    print(f"vx-averaged cleanup at y0 = 0 (Gaussian FWHM {VX_FWHM} m/s, |vx| <= {VX_CUT} m/s), SPB init")
    print("beta \\ detuning " + "".join(f"{d:>9.1f}" for d in DETUNINGS))
    best = None
    for b in BETAS:
        row = []
        for d in DETUNINGS:
            try:
                g = np.array([done[key((d, b, v, 0.0))] for v in VX_COARSE])
            except KeyError:
                row.append("      --"); continue
            a = wv @ g
            row.append(f"{100 * (1 - a):9.3f}")
            if best is None or a < best[0]:
                best = (a, b, d, g)
        print(f"{b:>15.1f} " + "".join(row))
    if best:
        a, b, d, g = best
        print(f"best: beta {b}, detuning {d} MHz: {eff(a)}; vs vx: " +
              ", ".join(f"{v:+.1f}:{100 * (1 - r):.2f}" for v, r in zip(VX_COARSE, g)))


if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "test":
        stage_test()
    elif stage == "coarse":
        stage_coarse()
    elif stage == "extend":
        stage_extend()
    elif stage == "yavg":
        vals = [float(v) for v in sys.argv[2:]]
        stage_yavg(list(zip(vals[::2], vals[1::2])))
    elif stage == "map":
        stage_map(float(sys.argv[2]), float(sys.argv[3]))
    else:
        summary()
