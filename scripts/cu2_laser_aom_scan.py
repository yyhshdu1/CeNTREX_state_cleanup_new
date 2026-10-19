"""CU2 two-dimensional scan at 200 V/cm in the experimental coordinates (laser frequency, AOM frequency).

The two CU2 beams (R(2) F'=3 and F'=4) come from one laser through a single-pass AOM, so
    detuning1 = detuning0 + 2*pi*(f_AOM - Delta_main),
with Delta_main the separation of the two main resonances (computed from the OBE system).

Stages (run one per call; every chunk is checkpointed and skipped on re-run):
    test    first three chunks of the coarse grid with different thread counts
    coarse  remaining chunks of the 2 MHz coarse grid
    fine    0.5 MHz grids around the best coarse local maxima
    confirm best fine points again with tighter tolerances

Usage (from the repository root):
    python scripts/cu2_laser_aom_scan.py test
    python scripts/cu2_laser_aom_scan.py coarse --threads 12 --chunk 40
    python scripts/cu2_laser_aom_scan.py fine --threads 12 --chunk 40
    python scripts/cu2_laser_aom_scan.py confirm
"""
from __future__ import annotations

import argparse
import json
import logging
import site
import sys
import time
import warnings
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
# Same convention as the notebooks: the installed centrex-tlf wheel shadows the source tree.
sys.path.insert(0, site.getsitepackages()[-1])

import numpy as np
import sympy as smp

from centrex_tlf import couplings, lindblad, states, transitions, utils
from scripts.scan_hdf import save_scan_hdf

logging.disable(logging.CRITICAL)

CHECKPOINT_DIR = ROOT / "tmp" / "cu2_laser_aom_scan"
HDF_FILE = ROOT / "results" / "SPB_J2_cleanup_scans.h5"
RUN_PREFIX = "laser_AOM_2D_200V"

# ---- grids (MHz) -----------------------------------------------------------
AOM_MIN, AOM_MAX = 185.0, 215.0
COARSE_D0 = np.arange(5.0, 23.0 + 1e-9, 2.0)            # 10 points
COARSE_AOM = np.arange(AOM_MIN, AOM_MAX + 1e-9, 2.0)    # 16 points
FINE_HALF_WIDTH, FINE_STEP, FINE_REGIONS = 2.0, 0.5, 3
TEST_THREADS = (20, 12, 8)
TEST_CHUNK = 40

# ---- experimental settings (baseline of docs/state_cleanup/main.tex) --------
E_FIELD = np.array([0.0, 0.0, 200.0])   # V/cm
B_FIELD = np.array([0.0, 0.0, 1e-3])    # G
VZ, VX = 184.0, 0.0                     # m/s; transverse velocity not considered
LASER_FREQUENCY = 1.103e15              # Hz
P0, P1 = 80e-3, 24e-3                   # W, F'=3 and F'=4 beams
EOM_FREQUENCY, BETA = 17.5 * 2 * np.pi * 1e6, 1.5
POL_FREQUENCY = 0.5 * 2 * np.pi * 1e6
SIGMA_Z, SIGMA_Y = 0.75e-3, 6e-3        # m; beam height 4*sigma_y = 2.4 cm
N_PASS, LOSS, FIRST_PASS, PASS_SPACING = 13, 0.10, 3.5e-3, 3.5e-3
DIPOLE = 2.6675506e-30
T_END = 0.05 / VZ

SCAN_SOLVER = dict(solver="dopri5", execution_mode="expanded_sparse", dt=2e-9, reltol=1e-6,
                   abstol=1e-8, maxiters=600_000, output="populations", output_when="final",
                   parallel=True)
CONFIRM_SOLVER = {**SCAN_SOLVER, "reltol": 1e-7, "abstol": 1e-9, "maxiters": 1_000_000}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def doppler_sign(t):
    """+1 / -1 for the alternating propagation direction of the passes along the molecular beam.

    In the two-prism multipass consecutive optical passes run in opposite directions, and the spatial
    order (1, 12, 3, 10, ...) alternates odd/even passes, so the direction flips from one pass to the
    next along the beam. The sign changes half-way between neighbouring pass centres.
    """
    return lindblad.alternating_sign(VZ * t, FIRST_PASS - PASS_SPACING / 2, PASS_SPACING)


def build(return_obe=False, doppler=False, sigma_y_m=None):
    """Build and prepare the 56-level CU2 model exactly as in the CU2 notebook.

    sigma_y_m sets the Gaussian beam height (sigma, m) *and* the power normalisation, so a different
    height is compared at constant power (default SIGMA_Y = 6 mm, i.e. a 2.4 cm full 4-sigma height).

    doppler=True adds a runtime parameter "vx" (m/s): a transverse velocity along the laser beams,
    which shifts both detunings by +-k*vx with the sign of the pass the molecule is in.
    """
    warnings.simplefilter("ignore")
    trans = [transitions.R2_F1_7o2_F3, transitions.R2_F1_7o2_F4]
    pols = [[couplings.polarization_X, couplings.polarization_Z]] * 2
    X, B = states.ElectronicState.X, states.ElectronicState.B
    ground_mains = [
        1 * states.CoupledBasisState(J=2, F=3, F1=5/2, mF=2, I1=1/2, I2=1/2, Ω=0, P=1, electronic_state=X),
        1 * states.CoupledBasisState(J=2, F=3, F1=5/2, mF=0, I1=1/2, I2=1/2, Ω=0, P=1, electronic_state=X),
    ]
    excited_mains = [
        1 * states.CoupledBasisState(J=3, F=3, F1=7/2, mF=3, I1=1/2, I2=1/2, Ω=1, P=-1, electronic_state=B),
        1 * states.CoupledBasisState(J=3, F=4, F1=7/2, mF=1, I1=1/2, I2=1/2, Ω=1, P=-1, electronic_state=B),
    ]
    sel = couplings.generate_transition_selectors(trans, pols, ground_mains=ground_mains,
                                                  excited_mains=excited_mains)
    obe = lindblad.generate_OBE_system_transitions(trans, sel, qn_compact=True, E=E_FIELD, B=B_FIELD,
                                                   retain_opposite_parity_levels=True, Jmax_X=5, Jmax_B=5)

    order = utils.multipass.multipass_prism_order(N_PASS)
    xlocs0 = [FIRST_PASS + i * PASS_SPACING for i in range(N_PASS)]
    pass_power = [(1.0 - LOSS) ** (p - 1) for p in order]
    sy = SIGMA_Y if sigma_y_m is None else float(sigma_y_m)
    intensities_per_watt = tuple(np.asarray(pass_power) / (2 * np.pi * SIGMA_Z * sy))

    params = lindblad.LindbladParameters()
    t = params.time()
    power0 = params.real("power0", P0); power1 = params.real("power1", P1)
    detuning0 = params.real("detuning0", 18.5 * 2 * np.pi * 1e6)
    detuning1 = params.real("detuning1", -6.5 * 2 * np.pi * 1e6)
    eom0 = params.real("eom0", EOM_FREQUENCY); eom1 = params.real("eom1", EOM_FREQUENCY)
    beta0 = params.real("beta0", BETA); beta1 = params.real("beta1", BETA)
    pol_frequency = params.real("pol_frequency", POL_FREQUENCY)
    y0 = params.real("y0", 0.0)
    sigma_z = params.real("sigma_z", SIGMA_Z); sigma_y = params.real("sigma_y", sy)
    profile = params.real("intensities_per_watt", intensities_per_watt)
    xlocs = params.real("xlocs", tuple(xlocs0)); ylocs = params.real("ylocs", tuple([0.0] * N_PASS))
    drives = []
    for k, (power, beta, eom) in enumerate(((power0, beta0, eom0), (power1, beta1, eom1))):
        mc = params.real(f"main_coupling{k}", float(abs(obe.couplings[k].main_coupling)))
        spatial = lindblad.multipass_2d_rabi(VZ * t, y0, profile, xlocs, ylocs, sigma_z, sigma_y, mc, DIPOLE)
        drives.append(power**0.5 * lindblad.phase_modulation(t, beta, eom) * spatial)
    doppler_per_mps = float(utils.detuning.velocity_to_detuning(1.0, LASER_FREQUENCY))
    if doppler:
        vx = params.real("vx", VX)
        shift = vx * doppler_per_mps * doppler_sign(t)     # same for both beams (same multipass)
    else:
        shift = VX * doppler_per_mps
    for s, drive, det in zip(sel, drives, (detuning0, detuning1)):
        params.bind(s.Ω, drive, finalize=False)
        params.bind(s.δ, det + shift, finalize=False)
    # Orthogonal polarizations through one polarization EOM: antiphase X/Z gates.
    z_gate = lindblad.square_wave(t, pol_frequency, 0.0)
    px0, pz0 = sel[0].polarization_symbols; px1, pz1 = sel[1].polarization_symbols
    params.bind(px0, 1.0 - z_gate, finalize=False); params.bind(pz0, z_gate, finalize=False)
    params.bind(px1, z_gate, finalize=False); params.bind(pz1, 1.0 - z_gate, finalize=False)
    params._finalize()
    prepared = lindblad.prepare_lindblad_problem(obe, params, backend="rust",
                                                 hamiltonian_representation="decomposed")

    ground_idx = np.asarray(states.QuantumSelector(J=2, electronic=X).get_indices(obe.QN), dtype=int).ravel()
    rho0 = np.zeros(obe.H_symbolic.shape, dtype=np.complex128)
    rho0[ground_idx, ground_idx] = 1.0 / ground_idx.size

    # Delta_main: F'=4 main resonance relative to the F'=3 main resonance (frame of beam 0).
    def find(st):
        lab = st.largest
        return next(i for i, q in enumerate(obe.QN) if q.largest == lab)
    g1, e1 = find(obe.couplings[1].ground_main), find(obe.couplings[1].excited_main)
    dd = smp.expand(obe.H_symbolic[e1, e1] - obe.H_symbolic[g1, g1])
    delta_main = complex(smp.N(dd.subs({s.δ: 0 for s in sel}))).real / (2 * np.pi * 1e6)

    info = dict(levels=len(obe.QN), main_couplings=[float(abs(c.main_coupling)) for c in obe.couplings],
                delta_main_MHz=delta_main, pass_power=pass_power)
    if return_obe:
        return prepared, rho0, ground_idx, info, obe
    return prepared, rho0, ground_idx, info


def coarse_points():
    d0, aom = np.meshgrid(COARSE_D0, COARSE_AOM, indexing="ij")
    pts = np.column_stack((d0.ravel(), aom.ravel()))
    order = np.random.default_rng(20261005).permutation(len(pts))  # chunks are random samples
    return pts[order]


def key(a, b):
    return (round(float(a), 6), round(float(b), 6))


def run_chunks(tag, points, prepared, rho0, ground_idx, info, threads, chunk, solver=SCAN_SOLVER,
               n_chunks=None, threads_per_chunk=None):
    """Run the points of `points` not yet in any checkpoint of `tag`, `chunk` at a time."""
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    done = {key(a, b) for a, b in load(tag)[0]}
    todo = np.asarray([p for p in points if key(*p) not in done]).reshape(-1, 2)
    if len(todo) < len(points):
        log(f"{tag}: {len(points) - len(todo)} of {len(points)} points already done")
    starts = list(range(0, len(todo), chunk))
    if n_chunks is not None:
        starts = starts[:n_chunks]
    for ic, start in enumerate(starts):
        path = CHECKPOINT_DIR / f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}_{ic:03d}.npz"
        nthreads = threads_per_chunk[ic] if threads_per_chunk else threads
        pts = todo[start:start + chunk]
        d1 = pts[:, 0] + pts[:, 1] - info["delta_main_MHz"]
        batch = 2 * np.pi * 1e6 * np.column_stack((pts[:, 0], d1))
        t0 = time.perf_counter()
        res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=["detuning0", "detuning1"],
                                      parameter_batch=batch, threads=nthreads, **solver)
        wall = time.perf_counter() - t0
        pops = np.asarray(res.values).real
        remaining = pops[:, ground_idx].sum(axis=-1)
        np.savez(path, points=pts, detuning1=d1, populations=pops, remaining=remaining, wall=wall,
                 threads=-1 if nthreads is None else nthreads, chunk=len(pts), reltol=solver["reltol"])
        best = int(np.argmin(remaining))
        log(f"{tag}: chunk {ic} ({len(pts)} pts, threads={nthreads}) {wall:.0f} s = "
            f"{wall/len(pts):.1f} s/trajectory; best here {100*(1-remaining[best]):.3f}% at "
            f"d0={pts[best,0]:.2f}, fAOM={pts[best,1]:.2f} MHz")


def load(tag):
    pts, rem, pops, walls = [], [], [], []
    for path in sorted(CHECKPOINT_DIR.glob(f"{tag}_*.npz")):
        z = np.load(path, allow_pickle=True)  # early checkpoints stored threads=None
        pts.append(z["points"]); rem.append(z["remaining"]); pops.append(z["populations"])
        th = z["threads"].item()
        walls.append((-1 if th is None else int(th), int(z["chunk"]), float(z["wall"])))
    if not pts:
        return np.empty((0, 2)), np.empty(0), np.empty((0, 0)), walls
    return np.concatenate(pts), np.concatenate(rem), np.concatenate(pops), walls


def metadata(info, solver, notes_extra=""):
    fixed = {
        "electric_field_V_per_cm": E_FIELD.tolist(), "magnetic_field_G": B_FIELD.tolist(),
        "interaction_length_m": 0.05, "molecular_vz_m_per_s": VZ, "molecular_vx_m_per_s": VX,
        "transverse_velocity": "not considered", "laser_frequency_Hz": LASER_FREQUENCY,
        "Fprime3_power_mW": P0 * 1e3, "Fprime4_power_mW": P1 * 1e3,
        "eom_frequency_MHz": EOM_FREQUENCY / (2 * np.pi * 1e6), "modulation_depth_rad": BETA,
        "eom": "single EOM before the AOM, identical modulation on both beams",
        "polarization_switch_frequency_MHz": POL_FREQUENCY / (2 * np.pi * 1e6),
        "polarization_scheme": "orthogonal polarizations through one polarization EOM (antiphase X/Z)",
        "aom": "single pass", "aom_range_MHz": [AOM_MIN, AOM_MAX],
        "Delta_main_MHz": info["delta_main_MHz"],
        "detuning_relationship": "detuning1 = detuning0 + 2*pi*(f_AOM - Delta_main)",
        "sigma_z_mm": SIGMA_Z * 1e3, "sigma_y_mm": SIGMA_Y * 1e3, "beam_height_4sigma_cm": 4 * SIGMA_Y * 1e2,
        "multipass_count": N_PASS, "loss_per_pass": LOSS, "multipass_relative_power": info["pass_power"],
        "main_couplings": info["main_couplings"],
    }
    model = {
        "script": "scripts/cu2_laser_aom_scan.py", "centrex_tlf_version": version("centrex-tlf"),
        "transitions": ["OpticalTransition(R(2) F1'=7/2 F'=3)", "OpticalTransition(R(2) F1'=7/2 F'=4)"],
        "retain_opposite_parity_levels": True, "number_of_levels": info["levels"],
        "initial_manifold": "X, J=2, uniform population", "initial_state_count": 20,
    }
    return fixed, model, {k: v for k, v in solver.items()}


def save_grid(run_name, d0_axis, aom_axis, lookup, info, solver, notes):
    grid = np.full((len(d0_axis), len(aom_axis)), np.nan)
    for i, a in enumerate(d0_axis):
        for j, b in enumerate(aom_axis):
            grid[i, j] = lookup.get((round(a, 6), round(b, 6)), np.nan)
    fixed, model, solv = metadata(info, solver)
    try:
        save_scan_hdf(HDF_FILE, run_name, grid,
                      axes={"Fprime3_detuning_MHz": (d0_axis, "MHz"), "AOM_frequency_MHz": (aom_axis, "MHz")},
                      quantity="final_population_remaining_in_X_J2", fixed_parameters=fixed, model=model,
                      solver=solv, notes=notes, overwrite=False)
        log(f"saved run {run_name!r} to {HDF_FILE.name}")
    except Exception as exc:  # keep going: the checkpoints hold the data
        log(f"could not save {run_name!r} ({exc}); data remain in {CHECKPOINT_DIR}")
    return grid


def local_maxima(grid, d0_axis, aom_axis, n):
    eff = 1.0 - grid
    cands = []
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v = eff[i, j]
            if np.isnan(v):
                continue
            nb = eff[max(i-1, 0):i+2, max(j-1, 0):j+2]
            if v >= np.nanmax(nb):
                cands.append((v, d0_axis[i], aom_axis[j]))
    cands.sort(reverse=True)
    chosen = []
    for v, a, b in cands:
        if all(abs(a - c[1]) > FINE_HALF_WIDTH or abs(b - c[2]) > FINE_HALF_WIDTH for c in chosen):
            chosen.append((v, a, b))
        if len(chosen) == n:
            break
    return chosen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["test", "coarse", "fine", "confirm"])
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--chunk", type=int, default=TEST_CHUNK)
    args = ap.parse_args()

    log(f"centrex-tlf {version('centrex-tlf')}; stage {args.stage}")
    t0 = time.perf_counter()
    prepared, rho0, ground_idx, info = build()
    log(f"model: {info['levels']} levels, Delta_main = {info['delta_main_MHz']:.3f} MHz, "
        f"build+prepare {time.perf_counter() - t0:.0f} s")
    (CHECKPOINT_DIR).mkdir(parents=True, exist_ok=True)
    (CHECKPOINT_DIR / "info.json").write_text(json.dumps(info, indent=1))
    pts = coarse_points()

    if args.stage == "test":
        run_chunks("coarse", pts, prepared, rho0, ground_idx, info, None, TEST_CHUNK,
                   n_chunks=len(TEST_THREADS), threads_per_chunk=list(TEST_THREADS))
        for th, n, wall in load("coarse")[3]:
            log(f"timing: threads={th:2d} chunk={n}: {wall/n:.1f} s/trajectory")
        return

    if args.stage == "coarse":
        run_chunks("coarse", pts, prepared, rho0, ground_idx, info, args.threads, args.chunk)
        cpts, crem, _, _ = load("coarse")
        lookup = {(round(a, 6), round(b, 6)): r for (a, b), r in zip(cpts, crem)}
        grid = save_grid(f"{RUN_PREFIX}_coarse_2MHz", COARSE_D0, COARSE_AOM, lookup, info, SCAN_SOLVER,
                         "Coarse 2 MHz grid in (laser detuning of F'=3 beam, AOM frequency) at 200 V/cm; "
                         "detuning1 follows from the single-pass AOM.")
        for v, a, b in local_maxima(grid, COARSE_D0, COARSE_AOM, FINE_REGIONS):
            log(f"coarse local maximum: {100*v:.3f}% at d0={a:.1f} MHz, fAOM={b:.1f} MHz")
        return

    cpts, crem, _, _ = load("coarse")
    lookup = {(round(a, 6), round(b, 6)): r for (a, b), r in zip(cpts, crem)}
    grid = np.full((len(COARSE_D0), len(COARSE_AOM)), np.nan)
    for i, a in enumerate(COARSE_D0):
        for j, b in enumerate(COARSE_AOM):
            grid[i, j] = lookup.get((round(a, 6), round(b, 6)), np.nan)
    if np.isnan(grid).any():
        raise SystemExit("coarse grid incomplete; run the coarse stage first")
    regions = local_maxima(grid, COARSE_D0, COARSE_AOM, FINE_REGIONS)

    if args.stage == "fine":
        offs = np.arange(-FINE_HALF_WIDTH, FINE_HALF_WIDTH + 1e-9, FINE_STEP)
        axes, new = [], []
        for v, a, b in regions:
            d0_axis = a + offs
            aom_axis = np.clip(b + offs, AOM_MIN, AOM_MAX); aom_axis = np.unique(aom_axis)
            axes.append((d0_axis, aom_axis, v, a, b))
            for x in d0_axis:
                for y in aom_axis:
                    kxy = key(x, y)
                    if kxy not in lookup and kxy not in {key(p, q) for p, q in new}:
                        new.append((x, y))
        new = np.asarray(new)
        log(f"fine: {len(regions)} regions, {len(new)} new points")
        run_chunks("fine", new, prepared, rho0, ground_idx, info, args.threads, args.chunk)
        fpts, frem, _, _ = load("fine")
        lookup.update({(round(a, 6), round(b, 6)): r for (a, b), r in zip(fpts, frem)})
        for k, (d0_axis, aom_axis, v, a, b) in enumerate(axes, start=1):
            g = save_grid(f"{RUN_PREFIX}_fine_0p5MHz_region{k}", d0_axis, aom_axis, lookup, info, SCAN_SOLVER,
                          f"0.5 MHz grid around coarse maximum d0={a:.1f} MHz, fAOM={b:.1f} MHz "
                          f"({100*v:.3f}%); coarse points reused where they coincide.")
            i, j = np.unravel_index(np.nanargmin(g), g.shape)
            log(f"region {k}: best {100*(1-g[i,j]):.3f}% at d0={d0_axis[i]:.2f}, fAOM={aom_axis[j]:.2f} MHz")
        return

    if args.stage == "confirm":
        fpts, frem, _, _ = load("fine")
        allpts = np.concatenate((cpts, fpts)); allrem = np.concatenate((crem, frem))
        best = allpts[np.argsort(allrem)[:3]]
        log("confirm: " + "; ".join(f"d0={a:.2f}, fAOM={b:.2f}" for a, b in best))
        run_chunks("confirm", best, prepared, rho0, ground_idx, info, None, len(best), solver=CONFIRM_SOLVER)
        qpts, qrem, _, _ = load("confirm")
        scan_rem = {(round(a, 6), round(b, 6)): r for (a, b), r in zip(allpts, allrem)}
        fixed, model, solv = metadata(info, CONFIRM_SOLVER)
        fixed["AOM_frequency_MHz_per_point"] = qpts[:, 1].tolist()
        fixed["remaining_at_scan_tolerance"] = [float(scan_rem[(round(a, 6), round(b, 6))]) for a, b in qpts]
        try:
            save_scan_hdf(HDF_FILE, f"{RUN_PREFIX}_confirm_rtol1e-7", qrem,
                          axes={"Fprime3_detuning_MHz": (qpts[:, 0], "MHz")},
                          quantity="final_population_remaining_in_X_J2", fixed_parameters=fixed, model=model,
                          solver=solv, notes="Best points of the 2D scan re-run with rtol 1e-7, atol 1e-9; "
                          "AOM frequency of each point in fixed_parameters.", overwrite=False)
        except Exception as exc:
            log(f"could not save confirm run ({exc}); data remain in {CHECKPOINT_DIR}")
        for (a, b), r in zip(qpts, qrem):
            log(f"confirmed: d0={a:.2f} MHz, fAOM={b:.2f} MHz: {100*(1-r):.4f}% "
                f"(scan tolerance: {100*(1-scan_rem[(round(a, 6), round(b, 6))]):.4f}%)")


if __name__ == "__main__":
    main()
