"""CU2 at the two 200 V/cm optima versus transverse velocity vx and transverse position y0 (coarse grid).

    python scripts/cu2_doppler_scan.py run        # R2 and R1, uniform and SPB initial states; checkpointed
    python scripts/cu2_doppler_scan.py summary    # on-axis value and averages from the checkpoints

Beam height 3 cm (sigma_y = 7.5 mm), as in the experiment for both CU1 and CU2, at constant power (80 + 24 mW):
cu2_laser_aom_scan.build(sigma_y_m=...) scales the intensity normalisation with the height. The optima
(delta0, f_AOM) are those found for the 2.4 cm beam. Doppler shift k*vx with the sign alternating pass by
pass, common to both beams (same multipass). Grid: vx -1...1 m/s in 0.2 m/s, y0 0...15 mm in 3 mm (11 x 6).
Initial states: all 20 J=2 sublevels (uniform) and the four J=2, m_J=0 levels with lens weights
(singlet 2/3, triplets 1/9 each). Averages as for CU1: vx Gaussian, FWHM 1.8 m/s, |vx| <= 1 m/s;
y0 uniform in |y0| <= 15 mm.
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cleanup_spb_doppler import describe  # noqa: E402
from scripts.cu1_eom_doppler_scan import VX_COARSE, Y_COARSE, vx_weights, y_weights  # noqa: E402
from scripts.cu2_laser_aom_scan import SCAN_SOLVER, T_END, build, log  # noqa: E402
from scripts.scan_hdf import save_scan_hdf  # noqa: E402

OUT = ROOT / "tmp" / "cu2_doppler"
REGIONS = {"R2": (11.5, 191.0), "R1": (11.5, 208.5)}     # (delta0/2pi, f_AOM) in MHz
SIGMA_Y_CU2 = 7.5e-3                                      # m, 3 cm full 4-sigma height (experiment, CU1 and CU2)
INITS = ("uniform", "spb")
CHUNK, THREADS = 40, 8
SLOTS = ["detuning0", "detuning1", "vx", "y0"]


def key(v, y):
    return (round(float(v), 6), round(float(y), 9))


def tag(reg, init):
    return f"doppler_cu2_{reg}_{init}_3cm_coarse"


def load(t):
    done = {}
    for f in sorted(OUT.glob(f"{t}_*.npz")):
        z = np.load(f)
        done.update({key(*p): r for p, r in zip(z["points"], z["remaining"])})
    return done


def setup():
    t0 = time.perf_counter()
    prepared, _, gidx, info, obe = build(return_obe=True, doppler=True, sigma_y_m=SIGMA_Y_CU2)
    dim = len(obe.QN)
    rho_u = np.zeros((dim, dim), dtype=np.complex128)
    rho_u[gidx, gidx] = 1.0 / len(gidx)
    rho_s = np.zeros((dim, dim), dtype=np.complex128)
    spb = [i for i in gidx if describe(obe.QN[i])["p_mj"].get(0, 0) > 0.5]
    for i in spb:
        rho_s[i, i] = 2 / 3 if describe(obe.QN[i])["spin"]["singlet"] > 0.5 else 1 / 9
    assert len(spb) == 4 and abs(np.trace(rho_s).real - 1) < 1e-12
    log(f"CU2 model ready ({time.perf_counter() - t0:.0f} s, {dim} levels, sigma_y {SIGMA_Y_CU2 * 1e3:.0f} mm); "
        f"Delta_main {info['delta_main_MHz']:.3f} MHz")
    return prepared, gidx, info, {"uniform": rho_u, "spb": rho_s}


def run_one(reg, init, prepared, gidx, info, rho0):
    from centrex_tlf import lindblad
    OUT.mkdir(parents=True, exist_ok=True)
    d0, faom = REGIONS[reg]
    d1 = d0 + faom - info["delta_main_MHz"]
    t = tag(reg, init)
    # small |vx| first, so that the first chunks already give the central part of the map
    pts = sorted(((v, y) for v in VX_COARSE for y in Y_COARSE), key=lambda p: (abs(p[0]), p[0], p[1]))
    done = load(t)
    todo = np.asarray([p for p in pts if key(*p) not in done]).reshape(-1, 2)
    log(f"{t}: delta0 {d0} MHz, delta1 {d1:+.2f} MHz; {len(todo)} of {len(pts)} points to compute")
    for c, s in enumerate(range(0, len(todo), CHUNK)):
        chunk = todo[s:s + CHUNK]
        batch = np.column_stack([np.full(len(chunk), 2 * np.pi * 1e6 * d0), np.full(len(chunk), 2 * np.pi * 1e6 * d1),
                                 chunk[:, 0], chunk[:, 1]])
        t0 = time.perf_counter()
        res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=SLOTS,
                                      parameter_batch=batch.astype(np.complex128), threads=THREADS, **SCAN_SOLVER)
        rem = np.asarray(res.values).real[:, gidx].sum(axis=1)
        np.savez(OUT / f"{t}_{time.strftime('%Y%m%d_%H%M%S')}_{c:03d}.npz", points=chunk, remaining=rem)
        done.update({key(*p): r for p, r in zip(chunk, rem)})
        log(f"{t}: chunk {c} ({len(chunk)} pts) {time.perf_counter() - t0:.0f} s; "
            f"{sum(key(*p) in done for p in pts)}/{len(pts)} done")
    grid = np.array([[done[key(v, y)] for y in Y_COARSE] for v in VX_COARSE])
    try:
        save_scan_hdf(ROOT / "results" / "SPB_J2_cleanup_scans.h5", t, grid,
                      axes={"vx_m_per_s": (VX_COARSE, "m/s"), "y0_mm": (Y_COARSE * 1e3, "mm")},
                      quantity="final_population_remaining_in_X_J2",
                      fixed_parameters={"region": reg, "delta0_MHz": d0, "f_AOM_MHz": faom, "delta1_MHz": d1,
                                        "E_V_per_cm": 200.0, "powers_mW": [80.0, 24.0],
                                        "beam_height_4sigma_cm": 4 * SIGMA_Y_CU2 * 100,
                                        "initial_state": init,
                                        "doppler": "k*vx common to both beams, sign alternating pass by pass",
                                        "y_symmetry": "results for y0 >= 0; symmetric in y0"},
                      model={"script": "scripts/cu2_doppler_scan.py", "levels": 56},
                      solver=dict(SCAN_SOLVER), notes="CU2 transverse velocity/position, coarse grid, 3 cm beam.",
                      overwrite=False)
        log(f"saved run {t!r}")
    except Exception as exc:                                  # keep the checkpoints in any case
        log(f"could not save {t!r}: {exc}")
    report(t, grid)


def report(t, g):
    wv, wy = vx_weights(VX_COARSE), y_weights(Y_COARSE)
    i0 = int(np.argmin(abs(VX_COARSE)))
    log(f"{t}: on axis {100 * (1 - g[i0, 0]):.3f} %; vx-averaged at y0=0 {100 * (1 - wv @ g[:, 0]):.2f} %; "
        f"vx and y0 averaged {100 * (1 - wv @ g @ wy):.2f} %")
    log(f"   y0=0 vs vx: " + ", ".join(f"{v:+.1f}:{100 * (1 - r):.2f}" for v, r in zip(VX_COARSE, g[:, 0])))
    log(f"   vx=0 vs y0: " + ", ".join(f"{y * 1e3:.0f}mm:{100 * (1 - r):.2f}" for y, r in zip(Y_COARSE, g[i0])))


def main():
    prepared, gidx, info, rhos = setup()
    for init in INITS:
        for reg in REGIONS:
            run_one(reg, init, prepared, gidx, info, rhos[init])


def summary():
    for init in INITS:
        for reg in REGIONS:
            t = tag(reg, init)
            done = load(t)
            try:
                g = np.array([[done[key(v, y)] for y in Y_COARSE] for v in VX_COARSE])
            except KeyError:
                log(f"{t}: {len(done)} of {len(VX_COARSE) * len(Y_COARSE)} points"); continue
            report(t, g)


if __name__ == "__main__":
    {"run": main, "summary": summary}[sys.argv[1]]()
