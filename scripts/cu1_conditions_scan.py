"""CU1 at a chosen laser power and beam height: detuning scan over the acceptance, then (vx, y0) maps.

    python scripts/cu1_conditions_scan.py detuning --power 50 --height 2         # coarse, 572 trajectories
    python scripts/cu1_conditions_scan.py map 19.0 [D ...] --power 50 --height 2  # 21 x 11 maps, stored
    python scripts/cu1_conditions_scan.py plot 19.0 --power 50 --height 2         # figure vs 80 mW, 3 cm

Same model and grids as scripts/cu1_3cm_scan.py (which is the 80 mW, 3 cm case): build_cu1 with the beam
height at constant power (intensity normalisation scales with sigma_y), laser power set through the runtime
parameter laser_power; 13 passes, 10 % loss, X/Z switching at 0.5 MHz, no EOM, Doppler +-k vx alternating pass
by pass; start in the 8 SPB-related J=1 levels with equal weights. Averages: vx Gaussian sigma 0.65 m/s,
|vx| <= 1 m/s; y0 Gaussian sigma 10.95 mm (lens + 0.5 m, 25 kV), |y0| <= 8.5 mm (detection laser 1.7 cm high).
Height is the full 4-sigma beam height. Results go to results/SPB_J1_cleanup_scans.h5.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cleanup_spb_doppler import CU1_SOLVER, build_cu1, describe  # noqa: E402
from scripts.cu1_3cm_scan import DETUNINGS, SY_BEAM, VX_C, VX_F, Y_C, Y_F, key  # noqa: E402
from scripts.cu2_laser_aom_scan import T_END, log  # noqa: E402
from scripts.plot_cleanup_maps_current import acceptance_mean  # noqa: E402
from scripts.scan_hdf import save_scan_hdf  # noqa: E402

HDF = ROOT / "results" / "SPB_J1_cleanup_scans.h5"
CHUNK, THREADS = 60, 20
SLOTS = ["laser_detuning", "vx", "y0", "laser_power"]


def cond_tag(power_mw, height_cm):
    return f"{height_cm:g}cm_{power_mw:g}mW".replace(".", "p")


def setup(height_cm):
    t0 = time.perf_counter()
    prepared, gidx, obe = build_cu1(return_obe=True, doppler=True, sigma_y_m=height_cm / 100 / 4)
    spb = [i for i in gidx if sum(v for m, v in describe(obe.QN[i])["p_mj"].items() if abs(m) == 1) > 0.5]
    rho0 = np.zeros((len(obe.QN),) * 2, dtype=np.complex128)
    rho0[spb, spb] = 1.0 / len(spb)
    log(f"CU1 model ready ({time.perf_counter() - t0:.0f} s), height {height_cm} cm, {len(spb)} SPB levels")
    return prepared, gidx, rho0


def load(out, tag):
    done = {}
    for f in sorted(out.glob(f"{tag}_*.npz")):
        z = np.load(f)
        done.update({key(p): r for p, r in zip(z["points"], z["remaining"])})
    return done


def run(out, tag, pts, prepared, gidx, rho0, power_mw):
    """pts: rows (detuning MHz, vx m/s, y0 m); checkpointed in out/."""
    from centrex_tlf import lindblad
    out.mkdir(parents=True, exist_ok=True)
    done = load(out, tag)
    todo = np.asarray([p for p in pts if key(p) not in done]).reshape(-1, 3)
    log(f"{tag}: {len(todo)} of {len(pts)} points to compute")
    for c, s in enumerate(range(0, len(todo), CHUNK)):
        chunk = todo[s:s + CHUNK]
        batch = np.column_stack([2 * np.pi * 1e6 * chunk[:, 0], chunk[:, 1], chunk[:, 2],
                                 np.full(len(chunk), power_mw * 1e-3)])
        t0 = time.perf_counter()
        res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=SLOTS,
                                      parameter_batch=batch.astype(np.complex128), threads=THREADS, **CU1_SOLVER)
        rem = np.asarray(res.values).real[:, gidx].sum(axis=1)
        np.savez(out / f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}_{c:03d}.npz", points=chunk, remaining=rem)
        done.update({key(p): r for p, r in zip(chunk, rem)})
        log(f"{tag}: chunk {c} ({len(chunk)} pts) {time.perf_counter() - t0:.0f} s")
    return done


def fixed(power_mw, height_cm, **extra):
    return {"laser_power_mW": power_mw, "beam_height_4sigma_cm": height_cm, "eom": "none",
            "initial_state": "spb (8 levels, equal weights)", "doppler": "k*vx, sign alternating pass by pass",
            "y_symmetry": "results for y0 >= 0; symmetric in y0", **extra}


def stage_detuning(power_mw, height_cm):
    ct = cond_tag(power_mw, height_cm)
    out = ROOT / "tmp" / f"cu1_{ct}"
    prepared, gidx, rho0 = setup(height_cm)
    pts = [(d, v, y) for d in DETUNINGS for v in VX_C for y in Y_C]
    done = run(out, "detuning", pts, prepared, gidx, rho0, power_mw)
    cube = np.array([[[done[key((d, v, y))] for y in Y_C] for v in VX_C] for d in DETUNINGS])
    print(f"CU1 {power_mw} mW, {height_cm} cm: detuning (MHz) | on axis | vx-averaged, y0=0 | mean in acceptance")
    wv = np.exp(-0.5 * (VX_C / 0.65) ** 2); wv[[0, -1]] *= 0.5; wv /= wv.sum()
    means = []
    for d, g in zip(DETUNINGS, cube):
        m = acceptance_mean(g, VX_C, Y_C, SY_BEAM)
        means.append(m)
        print(f"   {d:5.1f}  {100 * (1 - g[len(VX_C) // 2, 0]):8.3f} %  {100 * (1 - wv @ g[:, 0]):8.3f} %  "
              f"{100 * (1 - m):8.3f} %")
    best = DETUNINGS[int(np.argmin(means))]
    print(f"best detuning on the 0.5 MHz grid: {best} MHz, {100 * (1 - min(means)):.3f} %")
    name = f"detuning_vx_y_cu1_spb_{ct}_coarse"
    try:
        save_scan_hdf(HDF, name, cube,
                      axes={"laser_detuning_MHz": (DETUNINGS, "MHz"), "vx_m_per_s": (VX_C, "m/s"),
                            "y0_mm": (Y_C * 1e3, "mm")},
                      quantity="final_population_remaining_in_X_J1", fixed_parameters=fixed(power_mw, height_cm),
                      model={"script": "scripts/cu1_conditions_scan.py", "build": "cleanup_spb_doppler.build_cu1"},
                      solver=dict(CU1_SOLVER), notes=f"CU1 detuning scan, coarse vx x y0, {power_mw} mW, "
                      f"{height_cm} cm.", overwrite=False)
        log(f"saved {name}")
    except Exception as exc:
        log(f"could not save {name}: {exc}")


def stage_map(dets, power_mw, height_cm):
    ct = cond_tag(power_mw, height_cm)
    out = ROOT / "tmp" / f"cu1_{ct}"
    prepared, gidx, rho0 = setup(height_cm)
    for d in dets:
        tag = f"map_det{d:g}".replace(".", "p")
        done = run(out, tag, [(d, v, y) for v in VX_F for y in Y_F], prepared, gidx, rho0, power_mw)
        g = np.array([[done[key((d, v, y))] for y in Y_F] for v in VX_F])
        name = f"doppler_cu1_spb_{ct}_det{d:g}MHz".replace(".", "p")
        try:
            save_scan_hdf(HDF, name, g, axes={"vx_m_per_s": (VX_F, "m/s"), "y0_mm": (Y_F * 1e3, "mm")},
                          quantity="final_population_remaining_in_X_J1",
                          fixed_parameters=fixed(power_mw, height_cm, laser_detuning_MHz=d),
                          model={"script": "scripts/cu1_conditions_scan.py", "build": "cleanup_spb_doppler.build_cu1"},
                          solver=dict(CU1_SOLVER), notes=f"CU1 (vx, y0) map, {power_mw} mW, {height_cm} cm.",
                          overwrite=False)
            log(f"saved {name}")
        except Exception as exc:
            log(f"could not save {name}: {exc}")
        m = acceptance_mean(g, VX_F, Y_F, SY_BEAM)
        log(f"{name}: on axis {100 * (1 - g[len(VX_F) // 2, 0]):.3f} %; mean in acceptance {100 * (1 - m):.3f} %")


def stage_plot(d, power_mw, height_cm):
    """Side-by-side (vx, y0) maps: the present 80 mW, 3 cm beam and the chosen conditions, at detuning d."""
    import h5py
    import scripts.plot_cleanup_maps_current as pcm
    ct = cond_tag(power_mw, height_cm)
    runs = {"base": f"doppler_cu1_spb_3cm_det{d:g}MHz".replace(".", "p"),
            "new": f"doppler_cu1_spb_{ct}_det{d:g}MHz".replace(".", "p")}
    g = {}
    with h5py.File(HDF, "r") as h5:
        for k, name in runs.items():
            r = h5["runs"][name]
            g[k] = (r["axes"]["vx_m_per_s"][...], r["axes"]["y0_mm"][...] * 1e-3, r["results"][...])
    common = [r"R(1) $F_1'=5/2$, $F'=2$; 150 V/cm", rf"$\delta/2\pi=+{d:.1f}$ MHz, no EOM"]
    panels = [("base", "CU1", "(a) 80 mW, 3 cm beam", common + [r"80 mW, 3 cm high ($\sigma_y=7.5$ mm)",
                                                               r"start: 8 SPB-related $J=1$ levels"]),
              ("new", "CU1", f"(b) {power_mw:g} mW, {height_cm:g} cm beam",
               common + [rf"{power_mw:g} mW, {height_cm:g} cm high ($\sigma_y={height_cm * 10 / 4:g}$ mm)",
                         r"start: 8 SPB-related $J=1$ levels"])]
    note = [r"Dashed orange: detection acceptance $|y_0|\leq8.5$ mm (detection laser 1.7 cm high); the grid spans "
            r"the accepted $|v_x|\leq1$ m/s. Dotted orange: $1\sigma$ of the beam distribution at CU1 (lens + 0.5 m, "
            r"lens at 25 kV): $\sigma_{v_x}=0.65$ m/s, $\sigma_y=10.95$ mm; it weights the mean in acceptance. "
            r"White lines: contours of $\eta_\mathrm{c}$. On axis: $v_x=0$, $y_0=0$. Beam heights are full $4\sigma$ "
            r"widths; the power normalisation scales with the height.",
            r"13 passes, 10 % loss per pass, X/Z polarization switching at 0.5 MHz, Doppler shift $\pm k v_x$ "
            r"alternating pass by pass; single trajectories at $v_z=184$ m/s; grid 0.1 m/s $\times$ 1.5 mm, "
            r"computed for $y_0\geq0$ and mirrored. (a) centrex-tlf 0.2.5, (b) 0.3.0 (CU1 agrees to $10^{-6}$)."]
    pcm.figure(panels, g, f"cleanup_map_vx_y_CU1_{ct}_vs_3cm_80mW", (5.4, 4.9), note, 0.36)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["detuning", "map", "plot"])
    ap.add_argument("detunings", nargs="*", type=float)
    ap.add_argument("--power", type=float, required=True, help="laser power, mW")
    ap.add_argument("--height", type=float, required=True, help="full 4-sigma beam height, cm")
    a = ap.parse_args()
    if a.stage == "detuning":
        stage_detuning(a.power, a.height)
    elif a.stage == "map":
        stage_map(a.detunings, a.power, a.height)
    else:
        stage_plot(a.detunings[0], a.power, a.height)
