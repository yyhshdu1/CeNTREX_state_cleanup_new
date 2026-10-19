"""CU1 under the present conditions: 3 cm beam, 25 kV beam distribution, detection acceptance |y0| <= 8.5 mm.

    python scripts/cu1_3cm_scan.py detuning        # detuning x vx x y0 (coarse), averaged over the acceptance
    python scripts/cu1_3cm_scan.py map D [D ...]   # full 21 x 11 vx x y0 map at detuning(s) D (MHz), stored

Model: cleanup_spb_doppler.build_cu1(sigma_y_m=7.5e-3): R(1) F1'=5/2 F'=2, 150 V/cm, 80 mW, 3 cm full height at
constant power, 13 passes, 10 % loss, X/Z switching at 0.5 MHz, no EOM, Doppler +-k vx alternating pass by pass.
Start: the 8 SPB-related J=1 levels, equal weights. Means: vx Gaussian sigma 0.65 m/s, |vx| <= 1 m/s; y0
Gaussian sigma 10.95 mm (lens + 0.5 m, 25 kV), |y0| <= 8.5 mm (detection laser 1.7 cm high).
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cleanup_spb_doppler import CU1_SOLVER, build_cu1, describe  # noqa: E402
from scripts.cu2_laser_aom_scan import T_END, log  # noqa: E402
from scripts.plot_cleanup_maps_current import acceptance_mean  # noqa: E402
from scripts.scan_hdf import save_scan_hdf  # noqa: E402

OUT = ROOT / "tmp" / "cu1_3cm"
SIGMA_Y = 7.5e-3                                            # m, 3 cm full 4-sigma height
SY_BEAM = 10.95e-3                                          # m, molecular beam at CU1
DETUNINGS = np.round(np.arange(17.0, 23.01, 0.5), 3)        # MHz
VX_C = np.round(np.linspace(-1.0, 1.0, 11), 6)
Y_C = np.round(np.array([0.0, 3e-3, 6e-3, 9e-3]), 9)        # covers the 8.5 mm acceptance
VX_F = np.round(np.linspace(-1.0, 1.0, 21), 6)
Y_F = np.round(np.linspace(0.0, 15e-3, 11), 9)
CHUNK, THREADS = 40, 8


def key(p):
    return tuple(round(float(v), 9) for v in p)


def load(tag):
    done = {}
    for f in sorted(OUT.glob(f"{tag}_*.npz")):
        z = np.load(f)
        done.update({key(p): r for p, r in zip(z["points"], z["remaining"])})
    return done


def setup():
    t0 = time.perf_counter()
    prepared, gidx, obe = build_cu1(return_obe=True, doppler=True, sigma_y_m=SIGMA_Y)
    spb = [i for i in gidx if sum(v for m, v in describe(obe.QN[i])["p_mj"].items() if abs(m) == 1) > 0.5]
    rho0 = np.zeros((len(obe.QN),) * 2, dtype=np.complex128)
    rho0[spb, spb] = 1.0 / len(spb)
    log(f"CU1 model ready ({time.perf_counter() - t0:.0f} s), sigma_y {SIGMA_Y * 1e3} mm, {len(spb)} SPB levels")
    return prepared, gidx, rho0


def run(tag, pts, prepared, gidx, rho0):
    """pts: rows (detuning MHz, vx m/s, y0 m)."""
    from centrex_tlf import lindblad
    OUT.mkdir(parents=True, exist_ok=True)
    done = load(tag)
    todo = np.asarray([p for p in pts if key(p) not in done]).reshape(-1, 3)
    log(f"{tag}: {len(todo)} of {len(pts)} points to compute")
    for c, s in enumerate(range(0, len(todo), CHUNK)):
        chunk = todo[s:s + CHUNK]
        batch = np.column_stack([2 * np.pi * 1e6 * chunk[:, 0], chunk[:, 1], chunk[:, 2]])
        t0 = time.perf_counter()
        res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=["laser_detuning", "vx", "y0"],
                                      parameter_batch=batch.astype(np.complex128), threads=THREADS, **CU1_SOLVER)
        rem = np.asarray(res.values).real[:, gidx].sum(axis=1)
        np.savez(OUT / f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}_{c:03d}.npz", points=chunk, remaining=rem)
        done.update({key(p): r for p, r in zip(chunk, rem)})
        log(f"{tag}: chunk {c} ({len(chunk)} pts) {time.perf_counter() - t0:.0f} s")
    return done


def stage_detuning():
    prepared, gidx, rho0 = setup()
    pts = [(d, v, y) for d in DETUNINGS for v in VX_C for y in Y_C]
    done = run("detuning", pts, prepared, gidx, rho0)
    print("detuning (MHz) | on axis | vx-averaged, y0=0 | mean in acceptance (|vx|<=1 m/s, |y0|<=8.5 mm)")
    best = None
    for d in DETUNINGS:
        g = np.array([[done[key((d, v, y))] for y in Y_C] for v in VX_C])
        wv = np.exp(-0.5 * (VX_C / 0.65) ** 2); wv[[0, -1]] *= 0.5; wv /= wv.sum()
        m = acceptance_mean(g, VX_C, Y_C, SY_BEAM)
        print(f"   {d:5.1f}  {100 * (1 - g[5, 0]):8.3f} %  {100 * (1 - wv @ g[:, 0]):8.3f} %  {100 * (1 - m):8.3f} %")
        if best is None or m < best[0]:
            best = (m, d)
    print(f"best detuning on the 0.5 MHz grid: {best[1]} MHz, {100 * (1 - best[0]):.3f} %")


def stage_map(dets):
    prepared, gidx, rho0 = setup()
    for d in dets:
        tag = f"map_det{d:g}".replace(".", "p")
        done = run(tag, [(d, v, y) for v in VX_F for y in Y_F], prepared, gidx, rho0)
        g = np.array([[done[key((d, v, y))] for y in Y_F] for v in VX_F])
        name = f"doppler_cu1_spb_3cm_det{d:g}MHz".replace(".", "p")
        try:
            save_scan_hdf(ROOT / "results" / "SPB_J1_cleanup_scans.h5", name, g,
                          axes={"vx_m_per_s": (VX_F, "m/s"), "y0_mm": (Y_F * 1e3, "mm")},
                          quantity="final_population_remaining_in_X_J1",
                          fixed_parameters={"laser_detuning_MHz": d, "laser_power_mW": 80.0,
                                            "beam_height_4sigma_cm": 4 * SIGMA_Y * 100, "eom": "none",
                                            "initial_state": "spb (8 levels, equal weights)",
                                            "doppler": "k*vx, sign alternating pass by pass",
                                            "y_symmetry": "results for y0 >= 0; symmetric in y0"},
                          model={"script": "scripts/cu1_3cm_scan.py", "build": "cleanup_spb_doppler.build_cu1"},
                          solver=dict(CU1_SOLVER), notes="CU1 at the present 3 cm beam height.", overwrite=False)
            log(f"saved {name}")
        except Exception as exc:
            log(f"could not save {name}: {exc}")
        m = acceptance_mean(g, VX_F, Y_F, SY_BEAM)
        log(f"{name}: on axis {100 * (1 - g[10, 0]):.3f} %; mean in acceptance {100 * (1 - m):.3f} %")


if __name__ == "__main__":
    if sys.argv[1] == "detuning":
        stage_detuning()
    else:
        stage_map([float(v) for v in sys.argv[2:]])
