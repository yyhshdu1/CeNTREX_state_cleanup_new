"""CU1 under the conditions of the 2026-10 test on the SPB setup: does the model optimum move toward the measured one?

    python scripts/cu1_test_conditions.py timing  old          # 40 trajectories, time per trajectory
    python scripts/cu1_test_conditions.py scan    old|test     # detuning x vx x y0 at 50 mW
    python scripts/cu1_test_conditions.py power   old|test D [D ...]   # power scan (13, 30, 40, 50 mW) at detunings D
    python scripts/cu1_test_conditions.py analyze              # tables of both configurations
    python scripts/cu1_test_conditions.py levels               # main-pair resonance at 150 and 156 V/cm

Configurations (both: AP1 triplet start |J=1, F1=3/2, F=2, mF=0,+-1> with equal weights, 2 cm beam (4 sigma),
17 passes with 3.5 % loss per pass, X/Z switching at 0.5 MHz, no EOM, Doppler +-k vx alternating pass by pass):
    old   150 V/cm, v_z = 145 m/s   the earlier comparison with the test (report Sec. "Comparison with the model")
    test  156 V/cm, v_z = 156 m/s   the field of the test and the TOF forward velocity
Detected quantity: population left in the five |J=1, F1=3/2, F=2> sublevels (R(1) F1'=5/2 F'=3 sees only these),
relative to the start. vx weights: measured distribution, Gaussian sigma 1.28 m/s, |vx| <= 4 m/s. y0 is reported
on the axis and averaged with the lensed-beam weight of the other CU1 studies (Gaussian sigma 10.95 mm,
|y0| <= 8.5 mm); the y0 distribution of the test molecules (not focused by the lens) is not known.
The path is followed from 0 to the last pass + 3.5 mm (17 passes reach 59.5 mm, beyond the 5 cm of the
13-pass model). Detuning: from the main pair |J=1,F1=3/2,F=2,mF=1> -> |J'=2,F1'=5/2,F'=2,mF'=2,+> at that field.
Results: tmp/cu1_test_conditions/<config>/ (npz).
"""
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cleanup_spb_doppler import CU1_SOLVER, describe  # noqa: E402
from scripts.cu2_laser_aom_scan import (B_FIELD, DIPOLE, FIRST_PASS, LASER_FREQUENCY, PASS_SPACING,  # noqa: E402
                                        POL_FREQUENCY, SIGMA_Z, log)

OUT = ROOT / "tmp" / "cu1_test_conditions"
CONFIGS = {"old": dict(E=150.0, vz=145.0), "test": dict(E=156.0, vz=156.0)}
N_PASS, LOSS, SIGMA_Y = 17, 0.035, 5e-3                     # 2 cm full 4-sigma height
TRIPLET = [(1.5, 2, -1), (1.5, 2, 0), (1.5, 2, 1)]
DETUNINGS = np.round(np.arange(10.0, 24.01, 1.0), 6)          # MHz (about 1.4 s per trajectory)
VX = np.round(np.arange(-4.0, 4.001, 0.25), 6)                # m/s
Y0 = np.array([0.0, 4.5e-3, 8.5e-3])                          # m
POWERS = np.array([13.0, 30.0, 40.0, 50.0])                   # mW
DATA_POWER = {13.0: 0.150, 30.0: 0.063, 40.0: 0.045, 50.0: 0.032}   # measured 10-21, alignment C
SIGMA_VX, SY_BEAM, Y_ACC = 1.28, 10.95e-3, 8.5e-3
CHUNK, THREADS = 80, 20
SLOTS = ["laser_detuning", "vx", "y0", "laser_power"]


def build(E, vz):
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
    obe = lindblad.generate_OBE_system_transitions(trans, sel, qn_compact=True, E=np.array([0.0, 0.0, E]),
                                                   B=B_FIELD, retain_opposite_parity_levels=True)
    order = utils.multipass.multipass_prism_order(N_PASS)
    pass_power = [(1.0 - LOSS) ** (p - 1) for p in order]
    params = lindblad.LindbladParameters()
    t = params.time()
    power = params.real("laser_power", 50e-3)
    detuning = params.real("laser_detuning", 0.0)
    pol_frequency = params.real("polarization_frequency", POL_FREQUENCY)
    y0 = params.real("y0", 0.0)
    sigma_z = params.real("sigma_z", SIGMA_Z); sigma_y = params.real("sigma_y", SIGMA_Y)
    profile = params.real("intensities_per_watt", tuple(np.asarray(pass_power) / (2 * np.pi * SIGMA_Z * SIGMA_Y)))
    xlocs = params.real("xlocs", tuple(FIRST_PASS + i * PASS_SPACING for i in range(N_PASS)))
    ylocs = params.real("ylocs", tuple([0.0] * N_PASS))
    mc = params.real("main_coupling", float(abs(obe.couplings[0].main_coupling)))
    spatial = lindblad.multipass_2d_rabi(vz * t, y0, profile, xlocs, ylocs, sigma_z, sigma_y, mc, DIPOLE)
    params.bind(sel[0].Ω, power**0.5 * spatial, finalize=False)
    doppler_per_mps = float(utils.detuning.velocity_to_detuning(1.0, LASER_FREQUENCY))
    sign = lindblad.alternating_sign(vz * t, FIRST_PASS - PASS_SPACING / 2, PASS_SPACING)
    params.bind(sel[0].δ, detuning + params.real("vx", 0.0) * doppler_per_mps * sign, finalize=False)
    gate = lindblad.square_wave(t, pol_frequency, 0.0)
    px, pz = sel[0].polarization_symbols
    params.bind(px, 1.0 - gate, finalize=False); params.bind(pz, gate, finalize=False)
    params._finalize()
    prepared = lindblad.prepare_lindblad_problem(obe, params, backend="rust", hamiltonian_representation="decomposed")
    gidx = np.asarray(states.QuantumSelector(J=1, electronic=X).get_indices(obe.QN), dtype=int).ravel()
    labs = []
    for i in gidx:
        d = describe(obe.QN[i])
        labs.append((float(d["F1"]), int(d["F"]), int(round(d["mF"]))))
    idx = [gidx[labs.index(k)] for k in TRIPLET]
    rho0 = np.zeros((len(obe.QN),) * 2, dtype=np.complex128)
    rho0[idx, idx] = 1.0 / len(idx)
    t_end = (FIRST_PASS + N_PASS * PASS_SPACING) / vz
    return prepared, rho0, gidx, labs, t_end


def run(prepared, rho0, gidx, t_end, pts):
    """pts: rows (detuning MHz, vx m/s, y0 m, power mW) -> final J=1 sublevel populations."""
    from centrex_tlf import lindblad
    out = []
    for s in range(0, len(pts), CHUNK):
        c = pts[s:s + CHUNK]
        batch = np.column_stack([2 * np.pi * 1e6 * c[:, 0], c[:, 1], c[:, 2], 1e-3 * c[:, 3]]).astype(np.complex128)
        t0 = time.perf_counter()
        res = lindblad.parameter_scan(prepared, rho0, (0.0, t_end), parameter_slots=SLOTS, parameter_batch=batch,
                                      threads=THREADS, **CU1_SOLVER)
        out.append(np.asarray(res.values).real[:, gidx])
        log(f"{s + len(c)}/{len(pts)} points, chunk {time.perf_counter() - t0:.1f} s")
    return np.concatenate(out)


def f2_left(pops, labs):
    k = [i for i, l in enumerate(labs) if l[0] == 1.5 and l[1] == 2]
    return pops[..., k].sum(-1)


def vx_weights():
    w = np.exp(-0.5 * (VX / SIGMA_VX) ** 2)
    return w / w.sum()


def y_weights():
    """trapezoid-like weights on Y0 (>= 0, symmetric) for the Gaussian sigma 10.95 mm cut at 8.5 mm."""
    yy = np.linspace(0, Y_ACC, 401)
    g = np.exp(-0.5 * (yy / SY_BEAM) ** 2)
    w = np.zeros(len(Y0))
    for k, y in enumerate(yy):                  # linear interpolation onto the Y0 nodes
        j = min(np.searchsorted(Y0, y, side="right") - 1, len(Y0) - 2)
        a = (y - Y0[j]) / (Y0[j + 1] - Y0[j])
        w[j] += g[k] * (1 - a); w[j + 1] += g[k] * a
    return w / w.sum()


def timing(cfg):
    prepared, rho0, gidx, labs, t_end = build(**CONFIGS[cfg])
    pts = np.array([(19.0, v, 0.0, 50.0) for v in np.linspace(-4, 4, 40)])
    t0 = time.perf_counter()
    run(prepared, rho0, gidx, t_end, pts)
    log(f"{cfg}: {(time.perf_counter() - t0) / len(pts):.3f} s per trajectory (wall, {THREADS} threads)")


def scan(cfg):
    d = OUT / cfg
    d.mkdir(parents=True, exist_ok=True)
    f = d / "scan_50mW.npz"
    if f.exists():
        log(f"exists: {f}"); return
    prepared, rho0, gidx, labs, t_end = build(**CONFIGS[cfg])
    pts = np.array([(D, v, y, 50.0) for D in DETUNINGS for v in VX for y in Y0])
    pops = run(prepared, rho0, gidx, t_end, pts).reshape(len(DETUNINGS), len(VX), len(Y0), len(gidx))
    np.savez(f, pops=pops, labels=np.array(labs), detunings=DETUNINGS, vx=VX, y0=Y0, power_mW=50.0,
             config=json.dumps(dict(CONFIGS[cfg], n_pass=N_PASS, loss=LOSS, sigma_y=SIGMA_Y, start="AP1 triplet")))
    log(f"saved {f}")


def power(cfg, dets):
    d = OUT / cfg
    d.mkdir(parents=True, exist_ok=True)
    prepared, rho0, gidx, labs, t_end = build(**CONFIGS[cfg])
    for D in dets:
        f = d / f"power_D{D:g}.npz".replace(".", "p", 1).replace("pnpz", ".npz")
        if f.exists():
            log(f"exists: {f}"); continue
        pts = np.array([(D, v, y, P) for P in POWERS for v in VX for y in Y0])
        pops = run(prepared, rho0, gidx, t_end, pts).reshape(len(POWERS), len(VX), len(Y0), len(gidx))
        np.savez(f, pops=pops, labels=np.array(labs), detuning=D, vx=VX, y0=Y0, powers_mW=POWERS)
        log(f"saved {f}")


def means(grid):
    """grid[..., vx, y0] -> (on-axis y0 = 0, y-averaged), both vx-weighted."""
    wv, wy = vx_weights(), y_weights()
    axis = np.tensordot(grid[..., 0], wv, axes=([-1], [0]))
    avg = np.tensordot(np.tensordot(grid, wy, axes=([-1], [0])), wv, axes=([-1], [0]))
    return axis, avg


def analyze():
    res = {}
    for cfg in CONFIGS:
        f = OUT / cfg / "scan_50mW.npz"
        if not f.exists():
            continue
        z = np.load(f)
        labs = [tuple(x) for x in z["labels"].tolist()]
        left = f2_left(z["pops"], labs)                    # (D, vx, y0)
        axis, avg = means(left)
        print(f"\n== {cfg}: {CONFIGS[cfg]}, 50 mW, F=2 left (vx-weighted, sigma {SIGMA_VX} m/s)")
        print("  D/MHz   y0=0     y-avg    vx=0,y0=0")
        for k, D in enumerate(z["detunings"]):
            print(f"  {D:5.1f}  {axis[k]:.4f}  {avg[k]:.4f}  {left[k, int(np.argmin(abs(VX))), 0]:.5f}")
        ka, kb = int(np.argmin(axis)), int(np.argmin(avg))
        print(f"  optimum: y0=0 at {z['detunings'][ka]:g} MHz ({axis[ka]:.4f}); y-avg at {z['detunings'][kb]:g} MHz"
              f" ({avg[kb]:.4f}); at 19 MHz: {axis[list(z['detunings']).index(19.0)]:.4f} /"
              f" {avg[list(z['detunings']).index(19.0)]:.4f}")
        res[cfg] = dict(detunings=z["detunings"].tolist(), axis=axis.tolist(), avg=avg.tolist())
        for g in sorted((OUT / cfg).glob("power_D*.npz")):
            p = np.load(g)
            l = f2_left(p["pops"], labs)
            a, b = means(l)
            row = "  ".join(f"{P:g} mW: {x:.3f}/{y:.3f} ({DATA_POWER[P]:.3f})" for P, x, y in zip(p["powers_mW"], a, b))
            print(f"  power at {float(p['detuning']):g} MHz, y0=0/y-avg (data): {row}")
            res[cfg][f"power_{float(p['detuning']):g}"] = dict(axis=a.tolist(), avg=b.tolist())
    (OUT / "summary.json").write_text(json.dumps(res, indent=1))


def levels():
    """Main-pair transition frequency versus field, from the adiabatically tracked dressed levels."""
    from centrex_tlf import transitions
    from centrex_tlf.utils.plotting import calculate_transition_level_structure
    res = {}
    for E in (150.0, 156.0):
        st = calculate_transition_level_structure(transitions.R1_F1_5o2_F2, E=E, B=1e-3)
        g = [lv for lv in st.ground if (lv.F1, lv.F, lv.mF) == (1.5, 2, 1)]
        e = [lv for lv in st.excited if (lv.F1, lv.F, lv.mF, lv.P) == (2.5, 2, 2, 1)]
        assert len(g) == 1 and len(e) == 1, (len(g), len(e))
        res[E] = (g[0].energy_MHz, e[0].energy_MHz)
        print(f"E = {E:g} V/cm: ground {g[0].energy_MHz:.4f} MHz, excited {e[0].energy_MHz:.4f} MHz, "
              f"difference {e[0].energy_MHz - g[0].energy_MHz:.4f} MHz")
    d = [(res[E][1] - res[E][0]) for E in (150.0, 156.0)]
    print(f"main-pair resonance moves by {d[1] - d[0]:+.3f} MHz from 150 to 156 V/cm "
          f"(ground {res[156.0][0] - res[150.0][0]:+.3f}, excited {res[156.0][1] - res[150.0][1]:+.3f})")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "timing":
        timing(sys.argv[2])
    elif cmd == "scan":
        scan(sys.argv[2])
    elif cmd == "power":
        power(sys.argv[2], [float(x) for x in sys.argv[3:]])
    elif cmd == "analyze":
        analyze()
    elif cmd == "levels":
        levels()
