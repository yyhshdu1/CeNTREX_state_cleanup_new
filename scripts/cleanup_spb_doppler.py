"""Cleanup efficiency of the SPB-related states, and transverse velocity / position maps.

    python scripts/cleanup_spb_doppler.py states                    # task 1: per-state cleanup, CU1 and CU2
    python scripts/cleanup_spb_doppler.py doppler cu1               # task 9: vx x y map at the CU1 optimum
    python scripts/cleanup_spb_doppler.py doppler cu2 R1|R2         # task 9: vx x y map at a CU2 optimum

SPB-related states
    CU2 (after SPB): the four X, J=2, m_J=0 states that SPB acts on (nuclear-spin singlet and triplet);
                     whatever SPB does not transfer is what CU2 has to remove.
    CU1 (before SPB): the eight X, J=1, |m_J|=1 states, the states SPB transfers into; population left
                     there by CU1 mixes with the transferred molecules.
The states are identified from their m_J content in the uncoupled basis, not from labels. Because the
optical-Bloch equations are linear in rho0, the per-state results give any incoherent mixture.

Transverse motion: "vx" is the velocity along the laser beams; the Doppler shift k*vx changes sign from
one multipass pass to the next (scripts.cu2_laser_aom_scan.doppler_sign). "y0" is the transverse offset
along the beam height. Each (vx, y0) is one trajectory; the maps are stored for later averaging with any
beam distribution. y0 >= 0 only: the problem is symmetric in y0.
"""
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cu2_laser_aom_scan import (B_FIELD, CONFIRM_SOLVER, DIPOLE, FIRST_PASS, LASER_FREQUENCY,  # noqa: E402
                                        LOSS, N_PASS, PASS_SPACING, POL_FREQUENCY, SCAN_SOLVER, SIGMA_Z, T_END,
                                        VZ, doppler_sign, log)
from scripts.cu2_laser_aom_scan import build as build_cu2  # noqa: E402
from scripts.scan_hdf import save_scan_hdf  # noqa: E402

OUT_DIR = ROOT / "tmp" / "spb_states_doppler"
CU2_POINTS = {"R1": (11.5, 208.5), "R2": (11.5, 191.0)}    # (delta0/2pi, f_AOM) in MHz
CU1_E, CU1_POWER, CU1_SIGMA_Y, CU1_DETUNING = 150.0, 80e-3, 12e-3, 20.2   # V/cm, W, m, MHz
CU1_EOM_FREQUENCY = 1.6 * 2 * np.pi * 1e6                                    # rad/s, optional broadening EOM
CU1_SOLVER = dict(solver="dopri5", execution_mode="expanded_sparse", dt=2e-9, reltol=1e-7, abstol=1e-9,
                  maxiters=500_000, output="populations", output_when="final", parallel=True)
VX_GRID = np.round(np.linspace(-1.0, 1.0, 21), 6)          # m/s
Y_GRID = np.round(np.linspace(0.0, 15e-3, 11), 9)          # m
CHUNK, THREADS = 40, 8


def build_cu1(return_obe=False, doppler=False, sigma_y_m=None):
    """CU1 model of SPB_cleanup_J1-J2F2.ipynb: R(1) F1'=5/2 F'=2 at 150 V/cm, both parities (26 levels).

    sigma_y_m sets the beam height (sigma, m) and the power normalisation together, so a different height is
    compared at constant power (default CU1_SIGMA_Y = 12 mm, 4.8 cm full 4-sigma height).
    """
    from centrex_tlf import couplings, lindblad, states, transitions, utils
    warnings.simplefilter("ignore")
    X, B = states.ElectronicState.X, states.ElectronicState.B
    trans = [transitions.R1_F1_5o2_F2]
    sel = couplings.generate_transition_selectors(
        trans, [[couplings.polarization_X, couplings.polarization_Z]],
        ground_mains=[1 * states.CoupledBasisState(J=1, F=2, F1=3/2, mF=1, I1=1/2, I2=1/2, Ω=0, P=-1,
                                                   electronic_state=X)],
        excited_mains=[1 * states.CoupledBasisState(J=2, F=2, F1=5/2, mF=2, I1=1/2, I2=1/2, Ω=1, P=1,
                                                    electronic_state=B)])
    obe = lindblad.generate_OBE_system_transitions(trans, sel, qn_compact=True, E=np.array([0.0, 0.0, CU1_E]),
                                                   B=B_FIELD, retain_opposite_parity_levels=True)
    order = utils.multipass.multipass_prism_order(N_PASS)
    pass_power = [(1.0 - LOSS) ** (p - 1) for p in order]
    params = lindblad.LindbladParameters()
    t = params.time()
    power = params.real("laser_power", CU1_POWER)
    detuning = params.real("laser_detuning", CU1_DETUNING * 2 * np.pi * 1e6)
    pol_frequency = params.real("polarization_frequency", POL_FREQUENCY)
    y0 = params.real("y0", 0.0)
    sy = CU1_SIGMA_Y if sigma_y_m is None else float(sigma_y_m)
    sigma_z = params.real("sigma_z", SIGMA_Z); sigma_y = params.real("sigma_y", sy)
    profile = params.real("intensities_per_watt",
                          tuple(np.asarray(pass_power) / (2 * np.pi * SIGMA_Z * sy)))
    xlocs = params.real("xlocs", tuple(FIRST_PASS + i * PASS_SPACING for i in range(N_PASS)))
    ylocs = params.real("ylocs", tuple([0.0] * N_PASS))
    mc = params.real("main_coupling", float(abs(obe.couplings[0].main_coupling)))
    spatial = lindblad.multipass_2d_rabi(VZ * t, y0, profile, xlocs, ylocs, sigma_z, sigma_y, mc, DIPOLE)
    # optional phase EOM (CU1 has none in the baseline: eom_beta = 0 gives exactly the unmodulated drive)
    eom_beta = params.real("eom_beta", 0.0)
    eom_frequency = params.real("eom_frequency", CU1_EOM_FREQUENCY)
    params.bind(sel[0].Ω, power**0.5 * lindblad.phase_modulation(t, eom_beta, eom_frequency) * spatial,
                finalize=False)
    doppler_per_mps = float(utils.detuning.velocity_to_detuning(1.0, LASER_FREQUENCY))
    shift = params.real("vx", 0.0) * doppler_per_mps * doppler_sign(t) if doppler else 0.0
    params.bind(sel[0].δ, detuning + shift, finalize=False)
    gate = lindblad.square_wave(t, pol_frequency, 0.0)
    px, pz = sel[0].polarization_symbols
    params.bind(px, 1.0 - gate, finalize=False); params.bind(pz, gate, finalize=False)
    params._finalize()
    prepared = lindblad.prepare_lindblad_problem(obe, params, backend="rust",
                                                 hamiltonian_representation="decomposed")
    ground_idx = np.asarray(states.QuantumSelector(J=1, electronic=X).get_indices(obe.QN), dtype=int).ravel()
    return (prepared, ground_idx, obe) if return_obe else (prepared, ground_idx)


def uncoupled_content(state):
    """Probabilities of (mJ, m_Tl, m_F) in a field-dressed coupled state (I1 = Tl, I2 = 19F)."""
    amps = {}
    for amp, b in state.transform_to_uncoupled().data:
        k = (b.J, b.mJ, b.m1, b.m2)
        amps[k] = amps.get(k, 0) + amp
    p = {k: abs(a) ** 2 for k, a in amps.items()}
    norm = sum(p.values())
    return {k: v / norm for k, v in p.items()}


def describe(state):
    lab = state.largest
    content = uncoupled_content(state)
    p_mj = {}
    for (J, mJ, m1, m2), v in content.items():
        p_mj[mJ] = p_mj.get(mJ, 0) + v
    # nuclear-spin singlet / triplet character of the mJ = 0 part
    a = {}
    for amp, b in state.transform_to_uncoupled().data:
        if b.mJ == 0:
            a[(b.m1, b.m2)] = a.get((b.m1, b.m2), 0) + amp
    norm = sum(abs(v) ** 2 for v in a.values()) or 1.0
    spin = {"singlet": abs(a.get((0.5, -0.5), 0) - a.get((-0.5, 0.5), 0)) ** 2 / 2 / norm,
            "triplet m=0": abs(a.get((0.5, -0.5), 0) + a.get((-0.5, 0.5), 0)) ** 2 / 2 / norm,
            "triplet m=+1": abs(a.get((0.5, 0.5), 0)) ** 2 / norm,
            "triplet m=-1": abs(a.get((-0.5, -0.5), 0)) ** 2 / norm}
    best = max(content.items(), key=lambda kv: kv[1])
    return dict(F1=lab.F1, F=lab.F, mF=lab.mF, p_mj=p_mj, spin=spin,
                dominant=(best[0][1], best[0][2], best[0][3], best[1]))


def solve_states(prepared, dim, idx, slots, rows, solver, threads):
    """One trajectory per (initial state, parameter row); returns final populations [row, level]."""
    from centrex_tlf import lindblad
    rho_batch, par_batch = [], []
    for i, row in zip(idx, rows):
        rho = np.zeros((dim, dim), dtype=np.complex128)
        rho[i, i] = 1.0
        rho_batch.append(prepared.layout.pack(rho)); par_batch.append(row)
    res = lindblad.solve_lindblad_batch(prepared, np.asarray(rho_batch), (0.0, T_END), parameter_slots=slots,
                                        parameter_batch=np.asarray(par_batch, dtype=np.complex128),
                                        threads=threads, **solver)
    return np.asarray(res.values).real


def stage_states():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}

    # ---- CU2: J=2, m_J=0 ---------------------------------------------------------------------------
    t0 = time.perf_counter()
    prepared, _, gidx2, info, obe2 = build_cu2(return_obe=True)
    dim2 = len(obe2.QN)
    desc2 = {int(i): describe(obe2.QN[i]) for i in gidx2}
    spb2 = [i for i in gidx2 if desc2[int(i)]["p_mj"].get(0, 0) > 0.5]
    log(f"CU2 model ready ({time.perf_counter() - t0:.0f} s); J=2 m_J=0 states: {len(spb2)}")
    for i in spb2:
        d = desc2[int(i)]
        log(f"   F1={d['F1']}, F={d['F']}, mF={d['mF']:+.0f}: P(mJ=0)={d['p_mj'][0]:.3f}; "
            + ", ".join(f"{k} {v:.3f}" for k, v in d["spin"].items() if v > 0.02))
    idx, rows, tags = [], [], []
    for reg, (d0, faom) in CU2_POINTS.items():
        d1 = d0 + faom - info["delta_main_MHz"]
        for i in spb2:
            idx.append(i); rows.append(2 * np.pi * 1e6 * np.array([d0, d1])); tags.append((reg, int(i)))
    t0 = time.perf_counter()
    pops = solve_states(prepared, dim2, idx, ["detuning0", "detuning1"], rows, CONFIRM_SOLVER, threads=len(idx))
    log(f"CU2 per-state trajectories: {time.perf_counter() - t0:.0f} s")
    from centrex_tlf import states
    X = states.ElectronicState.X
    other = {J: np.asarray(states.QuantumSelector(J=J, electronic=X).get_indices(obe2.QN), dtype=int).ravel()
             for J in (1, 3, 4, 5)}
    cu2 = []
    for (reg, i), p in zip(tags, pops):
        d = desc2[i]
        cu2.append(dict(region=reg, F1=d["F1"], F=d["F"], mF=d["mF"], spin=d["spin"],
                        remaining=float(p[gidx2].sum()), **{f"J{J}": float(p[v].sum()) for J, v in other.items()}))
    results["cu2"] = cu2

    # ---- CU1: J=1, |m_J|=1 ----------------------------------------------------------------------------
    t0 = time.perf_counter()
    prepared1, gidx1, obe1 = build_cu1(return_obe=True)
    dim1 = len(obe1.QN)
    desc1 = {int(i): describe(obe1.QN[i]) for i in gidx1}
    spb1 = [i for i in gidx1 if desc1[int(i)]["p_mj"].get(1, 0) + desc1[int(i)]["p_mj"].get(-1, 0) > 0.5]
    log(f"CU1 model ready ({time.perf_counter() - t0:.0f} s, {dim1} levels); J=1 |m_J|=1 states: {len(spb1)}")
    rows1 = [np.array([2 * np.pi * 1e6 * CU1_DETUNING])] * len(spb1)
    t0 = time.perf_counter()
    pops1 = solve_states(prepared1, dim1, spb1, ["laser_detuning"], rows1, CU1_SOLVER, threads=len(spb1))
    log(f"CU1 per-state trajectories: {time.perf_counter() - t0:.0f} s")
    other1 = {J: np.asarray(states.QuantumSelector(J=J, electronic=X).get_indices(obe1.QN), dtype=int).ravel()
              for J in (0, 2, 3, 4)}
    cu1 = []
    for i, p in zip(spb1, pops1):
        d = desc1[int(i)]
        mJ, mTl, mFl, pd = d["dominant"]
        cu1.append(dict(F1=d["F1"], F=d["F"], mF=d["mF"], mJ=mJ, m_Tl=mTl, m_F19=mFl, p_dominant=pd,
                        p_absmJ1=d["p_mj"].get(1, 0) + d["p_mj"].get(-1, 0),
                        remaining=float(p[gidx1].sum()), **{f"J{J}": float(p[v].sum()) for J, v in other1.items()}))
    results["cu1"] = cu1
    (OUT_DIR / "spb_states.json").write_text(json.dumps(results, indent=1, default=float))

    for r in cu2:
        log(f"CU2 {r['region']}: F1={r['F1']}, F={r['F']}, mF={r['mF']:+.0f}: cleanup {100*(1-r['remaining']):.4f} %")
    for r in cu1:
        log(f"CU1 |mJ={r['mJ']:+d}, mTl={r['m_Tl']:+.1f}, mF={r['m_F19']:+.1f}> (F1={r['F1']}, F={r['F']}, "
            f"mF={r['mF']:+.0f}, p={r['p_dominant']:.2f}): cleanup {100*(1-r['remaining']):.4f} %")


def stage_doppler(stage, region=None, init="spb"):
    """vx x y0 map; init='spb' starts from the SPB-related states (equal weights), 'uniform' from all."""
    from centrex_tlf import lindblad, states
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"doppler_{stage}{'_' + region if region else ''}_{init}"
    if stage == "cu2":
        prepared, _, gidx, info, obe = build_cu2(return_obe=True, doppler=True)
        d0, faom = CU2_POINTS[region]
        fixed = {"detuning0": 2 * np.pi * 1e6 * d0, "detuning1": 2 * np.pi * 1e6 * (d0 + faom - info["delta_main_MHz"])}
        spb = [i for i in gidx if describe(obe.QN[i])["p_mj"].get(0, 0) > 0.5]
        solver = SCAN_SOLVER
    else:
        prepared, gidx, obe = build_cu1(return_obe=True, doppler=True)
        fixed = {"laser_detuning": 2 * np.pi * 1e6 * CU1_DETUNING}
        spb = [i for i in gidx if sum(v for m, v in describe(obe.QN[i])["p_mj"].items() if abs(m) == 1) > 0.5]
        solver = CU1_SOLVER
    start = spb if init == "spb" else list(gidx)
    rho0 = np.zeros((len(obe.QN),) * 2, dtype=np.complex128)
    if stage == "cu2" and init == "spb":
        # lens-focused J=2, m_J=0 composition (p2_f1_3o2_f1 fit report): 2/3 singlet, 1/9 each triplet
        for i in start:
            rho0[i, i] = 2 / 3 if describe(obe.QN[i])["spin"]["singlet"] > 0.5 else 1 / 9
    else:
        rho0[start, start] = 1.0 / len(start)
    assert abs(np.trace(rho0).real - 1) < 1e-12
    log(f"{tag}: model ready; initial state over {len(start)} sublevels")
    vx, y = np.meshgrid(VX_GRID, Y_GRID, indexing="ij")
    pts = np.column_stack((vx.ravel(), y.ravel()))
    done = {}
    for f in sorted(OUT_DIR.glob(f"{tag}_*.npz")):
        z = np.load(f)
        done.update({(round(a, 6), round(b, 9)): r for (a, b), r in zip(z["points"], z["remaining"])})
    todo = np.asarray([p for p in pts if (round(p[0], 6), round(p[1], 9)) not in done]).reshape(-1, 2)
    slots = list(fixed) + ["vx", "y0"]
    for c, start_i in enumerate(range(0, len(todo), CHUNK)):
        chunk = todo[start_i:start_i + CHUNK]
        batch = np.column_stack([np.full(len(chunk), v) for v in fixed.values()] + [chunk[:, 0], chunk[:, 1]])
        t0 = time.perf_counter()
        res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=slots,
                                      parameter_batch=batch.astype(np.complex128), threads=THREADS, **solver)
        rem = np.asarray(res.values).real[:, gidx].sum(axis=1)
        np.savez(OUT_DIR / f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}_{c:03d}.npz", points=chunk, remaining=rem)
        done.update({(round(a, 6), round(b, 9)): r for (a, b), r in zip(chunk, rem)})
        log(f"{tag}: chunk {c} ({len(chunk)} pts) {time.perf_counter() - t0:.0f} s; "
            f"{len(done)}/{len(pts)} done")
    grid = np.array([[done[(round(a, 6), round(b, 9))] for b in Y_GRID] for a in VX_GRID])
    hdf = ROOT / "results" / ("SPB_J2_cleanup_scans.h5" if stage == "cu2" else "SPB_J1_cleanup_scans.h5")
    try:
        save_scan_hdf(hdf, tag, grid, axes={"vx_m_per_s": (VX_GRID, "m/s"), "y0_mm": (Y_GRID * 1e3, "mm")},
                      quantity=f"final_population_remaining_in_X_J{2 if stage == 'cu2' else 1}",
                      fixed_parameters={"stage": stage, "region": region, "fixed_MHz": {k: v / (2e6 * np.pi) for k, v in fixed.items()},
                                        "initial_state": init, "initial_sublevels": len(start),
                                        "doppler": "k*vx, sign alternating pass by pass (two-prism multipass)",
                                        "y_symmetry": "results for y0 >= 0; symmetric in y0",
                                        "molecular_vz_m_per_s": VZ, "loss_per_pass": LOSS},
                      model={"script": "scripts/cleanup_spb_doppler.py", "levels": len(obe.QN)},
                      solver=dict(solver), notes="Transverse velocity and position map (task 9).", overwrite=False)
        log(f"saved run {tag!r} to {hdf.name}")
    except Exception as exc:
        log(f"could not save {tag!r}: {exc}")
    log(f"{tag}: box average (uniform in vx in [-1,1] m/s, |y0| <= 15 mm): cleanup {100*(1-grid.mean()):.3f} %; "
        f"on axis, vx=0: {100*(1-grid[10, 0]):.4f} %")


if __name__ == "__main__":
    if sys.argv[1] == "states":
        stage_states()
    else:
        stage_doppler(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] in CU2_POINTS else None,
                      init=sys.argv[-1] if sys.argv[-1] in ("spb", "uniform") else "spb")
