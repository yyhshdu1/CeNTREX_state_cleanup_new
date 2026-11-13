"""CU2 power split for a total of 60 mW, with 5 % loss per pass.

    python scripts/cu2_power_check.py split         # on axis: R1 and R2, P(F'=3) + P(F'=4) = 60 mW, both starts
    python scripts/cu2_power_check.py split2        # on axis: 15 + 45 and 10 + 50 mW
    python scripts/cu2_power_check.py acceptance R1 35 25   # vx x y0 grid in the detection acceptance at one setting
    python scripts/cu2_power_check.py opt 30 30     # (delta0, f_AOM) grid at one split, 6-point acceptance proxy
    python scripts/cu2_power_check.py fine 35 25    # finer (delta0, f_AOM) grid around R1 (fine2: beyond its edge)
    python scripts/cu2_power_check.py acceptance 13,206 35 25   # a set point given as delta0,f_AOM in MHz
    (append "spb" or "uniform" to acceptance/opt/fine to run one start only)
    python scripts/cu2_power_check.py beta          # EOM modulation depth at two set points, 35 + 25 mW, proxy
    python scripts/cu2_power_check.py splitproxy    # power split at the new set point NEW, proxy
    python scripts/cu2_power_check.py confirm       # a few points again with the tighter CONFIRM_SOLVER
    python scripts/cu2_power_check.py near 35 25 loss10   # grid stages at another loss per pass (%)
    (append "sy10" to acceptance to use a beam of sigma_y = 10 mm instead of 7.5 mm, same power)
    python scripts/cu2_power_check.py yfill sy6     # extra y0 = 7.5 mm row at NEW (edge of a small beam)
    python scripts/cu2_power_check.py near 35 25 sy6   # set point grid near NEW for a smaller beam, proxy
    python scripts/cu2_power_check.py efield 210    # model at another field: laser fixed (drift) and re-tuned, proxy
    python scripts/cu2_power_check.py vz 150        # forward velocity 150 m/s instead of 184 m/s at NEW, proxy
    python scripts/cu2_power_check.py loss 10       # another loss per pass (%) at NEW, 35 + 25 mW, proxy
    python scripts/cu2_power_check.py power         # more total power at NEW (same 7:5 ratio, and 80 + 24), proxy
    python scripts/cu2_power_check.py eomfreq 14.5,209 32 43 e217   # EOM frequency at one set point, proxy
    python scripts/cu2_power_check.py polfreq       # X/Z switching frequency at NEW, 35 + 25 mW, proxy
    python scripts/cu2_power_check.py polphase      # ... and its phase at the molecule's arrival (random in the lab)
    python scripts/cu2_power_check.py opt155 35 25 e155   # grid stages in a model built at another field (V/cm)
    python scripts/cu2_power_check.py opt217 32 43 e217 spb   # coarse grid at the measured CU2 field
    python scripts/cu2_power_check.py acc24 35 25 e155 12.5,208 6,200   # candidates on the reduced acceptance grid
    python scripts/cu2_power_check.py setpoints 44 31 e155 12.5,208 6,200   # the same set points, proxy only
    python scripts/cu2_power_check.py setpoints 32 43 e152 from155 12.5,200  # set points of 155 V/cm, field at 152
    python scripts/cu2_power_check.py perstate 12.5,208 35 25 e155      # each J=2 level alone, reduced grid
    python scripts/cu2_power_check.py summary

Forward velocity: base.VZ = 156 m/s (measured). Append e.g. "vz184" to any stage to run at another vz; those
results go to tmp/cu2_power/vz184/ under the same names (everything before 2026-11-09 is there).
Same model as scripts/cu2_doppler_scan.py (3 cm beam, 200 V/cm, 13 passes, set points R1/R2 found for 80 + 24 mW
with a 2.4 cm beam), except the loss per pass: the windows lose less than assumed, so LOSS = 0.05 instead of 0.10.
The powers are runtime parameters (power0 = F'=3 beam, power1 = F'=4 beam), so one prepared problem serves every
split. 80 + 24 mW at 5 % loss is included as a reference. Checkpoints in tmp/cu2_power (one file per call).
The acceptance grid is vx -1...1 m/s in 0.2 m/s and y0 0...9 mm in 3 mm (acceptance |y0| <= 8.5 mm), averaged
as in scripts/plot_cleanup_maps_current.py.
The set points were optimised on axis at 80 + 24 mW, so stage opt re-scans (delta0, f_AOM) at 60 mW with a cheap
proxy for the acceptance mean: vx = -0.6, 0, 0.6 m/s x y0 = 0, 6 mm, Gaussian weights in vx, plain mean in y0.
On the first three full grids it reads 0.1-0.5 % high but ranks them in the same order.
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import cu2_doppler_scan as c2  # noqa: E402
from scripts import cu2_laser_aom_scan as base  # noqa: E402
from scripts import vz_paths  # noqa: E402
from scripts.cu2_laser_aom_scan import CONFIRM_SOLVER, SCAN_SOLVER, T_END, log  # noqa: E402

LOSS = 0.05
base.LOSS = LOSS                                            # read by base.build() when the model is prepared
OUT = ROOT / "tmp" / "cu2_power"
SLOTS = ["detuning0", "detuning1", "vx", "y0", "power0", "power1", "beta0", "beta1", "pol_frequency", "pol_phase",
         "eom0", "eom1"]
TOTAL = 60
SPLITS = [(TOTAL - p1, p1) for p1 in (10, 15, 20, 25, 30, 35, 40)]   # mW, (F'=3, F'=4)
SPLITS2 = [(15, 45), (10, 50)]                                     # extension: uniform start still rising
REFERENCE = [(80, 24)]
VX = np.round(np.arange(-1.0, 1.0 + 1e-9, 0.2), 6)
Y0 = np.array([0.0, 3e-3, 6e-3, 9e-3])
VX24 = VX[::2]                                                     # reduced grid: <= 0.08 % from the 44-point mean
THREADS = 16
GRIDS = {                                                          # (delta0/2pi, f_AOM) in MHz
    "opt": ([5.0, 8.0, 11.5, 15.0, 18.5], [186.0, 191.0, 197.0, 203.0, 208.5, 214.0]),   # AOM range 200 +- 15 MHz
    "fine": ([10.0, 11.0, 11.5, 12.0, 13.0], [206.0, 207.5, 208.5, 209.5, 211.0]),      # around R1
    "fine2": ([12.5, 13.0, 13.5, 14.0, 15.0], [203.0, 204.5, 206.0, 207.5]),            # fine optimum at its edge
    "near": ([11.5, 12.0, 12.5, 13.0], [206.0, 207.5, 208.5]),                          # near NEW, other beams
    "opt155": ([3.0, 6.0, 9.0, 12.5, 16.0], [192.0, 196.0, 200.0, 204.0, 208.0, 212.0]),  # coarse, at 155 V/cm
    "fine155": ([11.5, 12.5, 13.5], [198.0, 200.0, 202.0, 208.0, 210.0, 212.0]),            # both good f_AOM regions
    "opt217": ([3.0, 6.0, 9.0, 12.5, 16.0], [190.0, 195.0, 200.0, 205.0, 210.0, 215.0]),  # coarse, at 217 V/cm
    "fine217": ([11.5, 12.5, 13.5, 14.5, 16.0, 17.0], [205.0, 207.0, 209.0, 210.0, 211.0, 213.0]),  # around both optima
}
BETA_POINTS = [(12.5, 207.5), (13.0, 206.0)]                       # best proxy points of Steps 4-5
BETAS = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5]                           # one EOM for both beams: beta0 = beta1
NEW = (12.5, 207.5)                                                # best full-grid point for the SPB start (Step 6)
PROXY_SPLITS = [(25, 35), (30, 30), (35, 25), (40, 20), (45, 15)]
CONFIRM_POINTS = [(NEW, 0.0), (NEW, 0.6), ((13.0, 206.0), 0.0), ((13.0, 206.0), 0.6)]   # (set point, vx), y0 = 0
MORE_POWER = [(47, 33), (58, 42), (80, 24)]                        # 80, 100 mW at 35:25, and the old 80 + 24 mW
EOM_FREQS = [13.0, 14.0, 15.0, 15.5, 16.0, 17.0, 17.5, 18.5]       # MHz, stage eomfreq (proxy)
POL_FREQS = [0.1, 0.25, 0.5, 1.0, 2.0]                             # MHz, X/Z square-wave switching
POL_PHASE_FREQS, POL_PHASES = [0.25, 0.5, 1.0], [0.5, 1.0, 1.5]    # MHz, units of pi (phase 0 is stage polfreq)
PROXY_VX, PROXY_Y0 = [-0.6, 0.0, 0.6], [0.0, 6e-3]
CHUNK = 16
info = None


def slot_values(row):
    """(region or (delta0, f_AOM), P0 mW, P1 mW, vx, y0[, beta[, pol MHz[, phase/pi[, EOM MHz]]]]) -> SLOTS."""
    reg, p0, p1, vx, y0, *rest = row
    d0, faom = c2.REGIONS[reg] if isinstance(reg, str) else reg
    d1 = d0 + faom - info["delta_main_MHz"]
    beta = rest[0] if rest else base.BETA
    pol = 2 * np.pi * 1e6 * rest[1] if len(rest) > 1 else base.POL_FREQUENCY
    phase = np.pi * rest[2] if len(rest) > 2 else 0.0
    eom = 2 * np.pi * 1e6 * rest[3] if len(rest) > 3 else base.EOM_FREQUENCY      # one EOM for both beams
    return [2 * np.pi * 1e6 * d0, 2 * np.pi * 1e6 * d1, vx, y0, p0 * 1e-3, p1 * 1e-3, beta, beta, pol, phase, eom, eom]


def run_states(name, jobs, prepared, gidx, full=False):
    """jobs: (level index, row) -> remaining J=2 population starting from that one level (full: final population of
    every level instead); one checkpoint per call."""
    from centrex_tlf import lindblad
    f = OUT / f"{name}.npz"
    if f.exists():
        return np.load(f, allow_pickle=True)["populations" if full else "remaining"]
    dim = len(c2.OBE.QN)
    rho_batch = []
    for i, _ in jobs:
        rho = np.zeros((dim, dim), dtype=np.complex128)
        rho[i, i] = 1.0
        rho_batch.append(prepared.layout.pack(rho))
    t0 = time.perf_counter()
    res = lindblad.solve_lindblad_batch(prepared, np.asarray(rho_batch), (0.0, T_END), parameter_slots=SLOTS,
                                        parameter_batch=np.asarray([slot_values(r) for _, r in jobs],
                                                                   dtype=np.complex128),
                                        threads=THREADS, **{**SCAN_SOLVER, "maxiters": 1_000_000})
    pops = np.asarray(res.values).real
    rem = pops[:, gidx].sum(axis=1)
    np.savez(f, levels=np.array([i for i, _ in jobs]), rows=np.array([repr(r) for _, r in jobs]), remaining=rem,
             wall=time.perf_counter() - t0, loss=LOSS, **({"populations": pops} if full else {}))
    log(f"{name}: {len(jobs)} trajectories in {time.perf_counter() - t0:.0f} s")
    return pops if full else rem


def run(name, rows, prepared, gidx, rho0, solver=SCAN_SOLVER):
    """rows: (region or (delta0, f_AOM) in MHz, P0 mW, P1 mW, vx, y0[, beta[, pol MHz[, phase/pi]]]) -> remaining J=2
    population; one checkpoint per call."""
    from centrex_tlf import lindblad
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / f"{name}.npz"
    if f.exists():
        return np.load(f, allow_pickle=True)["remaining"]
    batch = [slot_values(r) for r in rows]
    t0 = time.perf_counter()
    res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=SLOTS,
                                  parameter_batch=np.asarray(batch, dtype=np.complex128), threads=THREADS,
                                  **solver)
    rem = np.asarray(res.values).real[:, gidx].sum(axis=1)
    np.savez(f, rows=np.array([repr(r) for r in rows]), remaining=rem, wall=time.perf_counter() - t0, loss=LOSS)
    log(f"{name}: {len(rows)} trajectories in {time.perf_counter() - t0:.0f} s")
    return rem


def split_rows():
    return [(reg, p0, p1, 0.0, 0.0) for reg in ("R1", "R2") for p0, p1 in SPLITS + REFERENCE]


def split2_rows():
    return [(reg, p0, p1, 0.0, 0.0) for reg in ("R1", "R2") for p0, p1 in SPLITS2]


def acc_rows(reg, p0, p1):
    return [(reg, p0, p1, v, y) for v in VX for y in Y0]


def opt_rows(grid, p0, p1):
    d0s, aoms = GRIDS[grid]
    return [((d0, fa), p0, p1, v, y) for d0 in d0s for fa in aoms for v in PROXY_VX for y in PROXY_Y0]


def beta_rows():
    return [(pt, 35, 25, v, y, b) for pt in BETA_POINTS for b in BETAS for v in PROXY_VX for y in PROXY_Y0]


def splitproxy_rows():
    return [(NEW, p0, p1, v, y) for p0, p1 in PROXY_SPLITS for v in PROXY_VX for y in PROXY_Y0]


def confirm_rows():
    return [(pt, 35, 25, v, 0.0) for pt, v in CONFIRM_POINTS]


def main_line_shift_MHz(E):
    """Shift of the F'=3 main line (X J=2 F1=5/2 F=3 mF=2 -> B J'=3 F1'=7/2 F'=3 mF'=3) from 200 V/cm to E."""
    import warnings
    from centrex_tlf import hamiltonian, states
    warnings.simplefilter("ignore")
    X, B = states.ElectronicState.X, states.ElectronicState.B
    g = states.CoupledBasisState(J=2, F1=5/2, F=3, mF=2, I1=1/2, I2=1/2, Ω=0, P=1, electronic_state=X)
    e = states.CoupledBasisState(J=3, F1=7/2, F=3, mF=3, I1=1/2, I2=1/2, Ω=1, P=-1, electronic_state=B)

    def line(field):
        hx = hamiltonian.generate_reduced_X_hamiltonian([g], E=np.array([0, 0, field]), B=base.B_FIELD, Jmin=0, Jmax=5)
        hb = hamiltonian.generate_reduced_B_hamiltonian([e], E=np.array([0, 0, field]), B=base.B_FIELD, Jmin=1, Jmax=5)
        return (np.real(hb.H[0, 0]) - np.real(hx.H[0, 0])) / (2 * np.pi * 1e6)
    return line(E) - line(200.0)


def efield_rows(E):
    """Laser fixed while the field drifts (delta = omega_L - omega_line, so delta0 -> delta0 - shift), then re-tuned."""
    s0 = main_line_shift_MHz(E)
    return [((NEW[0] - s0, NEW[1]), 35, 25, v, y) for v in PROXY_VX for y in PROXY_Y0] +            [(NEW, 35, 25, v, y) for v in PROXY_VX for y in PROXY_Y0]


def proxy(rem):
    """rem over (setting, vx, y0) of the proxy points -> proxy for the remaining fraction in the acceptance."""
    w = np.exp(-0.5 * (np.asarray(PROXY_VX) / 0.65) ** 2); w /= w.sum()
    return np.einsum("v,svy->s", w, rem.reshape(-1, len(PROXY_VX), len(PROXY_Y0))) / len(PROXY_Y0)


def pt_name(pt):
    return f"d{pt[0]:g}_f{pt[1]:g}".replace(".", "p")


def parse_pt(s):
    return tuple(float(x) for x in s.split(","))


def perstate(stem, pt, p0, p1, prepared, gidx):
    """Remaining J=2 from each J=2 level alone at one set point, on the reduced acceptance grid (VX24 x Y0).
    The OBEs are linear in rho0, so any incoherent distribution over these levels is a weighted sum of the results.
    Levels with mF < 0 are taken from their mF > 0 partner (reflection through a plane containing z maps mF -> -mF
    and leaves the X/Z light and k*vx unchanged; only the 1 mG field breaks it); checked at vx = -0.6, 0.6, y0 = 0."""
    from scripts.cleanup_spb_doppler import describe
    from scripts.plot_cleanup_maps_current import SY, acceptance_mean
    desc = {int(i): describe(c2.OBE.QN[i]) for i in gidx}
    lab = {i: (float(d["F1"]), float(d["F"]), int(round(d["mF"]))) for i, d in desc.items()}
    assert len(set(lab.values())) == len(lab), "J=2 labels are not unique"
    pos = [i for i in lab if lab[i][2] >= 0]
    partner = {i: next(j for j in pos if lab[j] == (lab[i][0], lab[i][1], -lab[i][2])) for i in lab if lab[i][2] < 0}
    check_vx = (-0.6, 0.6)
    r_check = run_states(f"{stem}_symcheck", [(i, (pt, p0, p1, v, 0.0)) for i in partner for v in check_vx],
                         prepared, gidx)
    jobs = [(i, (pt, p0, p1, v, y)) for i in pos for v in VX24 for y in Y0]
    r = np.concatenate([run_states(f"{stem}_{k:03d}", jobs[k:k + CHUNK], prepared, gidx)
                        for k in range(0, len(jobs), CHUNK)])
    n = len(VX24) * len(Y0)
    g = {i: r[k * n:(k + 1) * n].reshape(len(VX24), len(Y0)) for k, i in enumerate(pos)}
    for k, (i, v) in enumerate((i, v) for i in partner for v in check_vx):
        ref = g[partner[i]][int(np.argmin(abs(VX24 - v))), 0]
        log(f"symmetry {lab[i]} vs {lab[partner[i]]}, vx {v:+.1f}: remaining {r_check[k]:.4e} vs {ref:.4e} "
            f"({100 * (r_check[k] - ref):+.4f} % of J=2)")
    for i, j in partner.items():
        g[i] = g[j]
    levels = sorted(lab, key=lambda i: lab[i])
    acc = {i: acceptance_mean(g[i], VX24, Y0, SY["CU2"]) for i in levels}
    for i in levels:
        log(f"   F1={lab[i][0]:.1f} F={lab[i][1]:.0f} mF={lab[i][2]:+d}  P(mJ=0) {desc[i]['p_mj'].get(0, 0):.2f}: "
            f"on axis {100 * (1 - g[i][len(VX24) // 2, 0]):.2f} % (vx 0.2), acceptance {100 * (1 - acc[i]):.2f} %")
    w_spb = {i: (2 / 3 if desc[i]["spin"]["singlet"] > 0.5 else 1 / 9) for i in levels
             if desc[i]["p_mj"].get(0, 0) > 0.5}
    log(f"{stem}: SPB start {100 * (1 - sum(w * acc[i] for i, w in w_spb.items())):.2f} %, "
        f"uniform start {100 * (1 - np.mean([acc[i] for i in levels])):.2f} % (weighted sums, reduced grid)")
    np.savez(OUT / f"{stem}_result.npz", level=np.array(levels), F1=np.array([lab[i][0] for i in levels]),
             F=np.array([lab[i][1] for i in levels]), mF=np.array([lab[i][2] for i in levels]),
             p_mJ0=np.array([desc[i]["p_mj"].get(0, 0) for i in levels]),
             singlet=np.array([desc[i]["spin"]["singlet"] for i in levels]),
             from_partner=np.array([i in partner for i in levels]), vx=VX24, y0=Y0,
             remaining=np.array([g[i] for i in levels]), acceptance_remaining=np.array([acc[i] for i in levels]),
             set_point_MHz=np.array(pt), powers_mW=np.array([p0, p1]), E_V_cm=float(base.E_FIELD[2]), loss=LOSS)


def chunked(stem, rows, prepared, gidx, rho0):
    return np.concatenate([run(f"{stem}_{k:03d}", rows[k:k + CHUNK], prepared, gidx, rho0)
                           for k in range(0, len(rows), CHUNK)])


def print_rows(rows, rem, init):
    for (reg, p0, p1, _, _), r in zip(rows, rem):
        log(f"on axis, loss {LOSS:.0%}, {init:7s} {reg}, {p0:2d} + {p1:2d} mW: {100 * (1 - r):.3f} %  "
            f"(remaining {r:.2e})")


def setup():
    global info
    prepared, gidx, info, rhos = c2.setup()
    log(f"loss per pass {LOSS:.0%}; relative pass powers {np.round(info['pass_power'], 3).tolist()}")
    return prepared, gidx, rhos


def main(argv):
    stage = argv[0]
    if stage == "efield":
        base.E_FIELD = np.array([0.0, 0.0, float(argv[1])])                  # read by base.build()
    global LOSS, T_END, OUT
    if stage == "loss":
        LOSS = base.LOSS = float(argv[1]) / 100                              # read by base.build()
    if stage == "vz":
        base.VZ = float(argv[1])                                             # read by build() and doppler_sign()
        T_END = base.T_END = 0.05 / base.VZ
    vzt = next((float(a[2:]) for a in argv[1:] if a.startswith("vz") and a[2:].replace(".", "").isdigit()), None)
    argv = [argv[0]] + [a for a in argv[1:] if not (a.startswith("vz") and a[2:].replace(".", "").isdigit())]
    if vzt is not None:            # any stage: forward velocity (m/s); results go to tmp/cu2_power/vz<vz>/
        base.VZ = vzt
        T_END = base.T_END = 0.05 / base.VZ
    sy = next((float(a[2:]) for a in argv if a.startswith("sy")), None)      # beam sigma_y in mm
    argv = [a for a in argv if not a.startswith("sy")]
    if sy is not None:
        c2.SIGMA_Y_CU2 = sy * 1e-3                                           # read by c2.setup()
    tag = "" if sy is None else f"_sy{sy:g}".replace(".", "p")
    ef = next((float(a[1:]) for a in argv[1:] if a[0] == "e" and a[1:].replace(".", "").isdigit()), None)
    argv = [argv[0]] + [a for a in argv[1:] if not (a[0] == "e" and a[1:].replace(".", "").isdigit())]
    if ef is not None:                                                       # model at another field, V/cm
        base.E_FIELD = np.array([0.0, 0.0, ef])                              # read by base.build()
        tag += f"_e{ef:g}".replace(".", "p")
    if vzt is not None:
        OUT = vz_paths.vz_dir(OUT, vzt)                                      # same names, own folder
        OUT.mkdir(parents=True, exist_ok=True)
    inits = [a for a in argv[1:] if a in c2.INITS] or list(c2.INITS)
    argv = [a for a in argv if a not in c2.INITS]
    lp = next((float(a[4:]) for a in argv[1:] if a.startswith("loss")), None)   # e.g. "loss10" for the grid stages
    argv = [a for a in argv if not (a.startswith("loss") and a != "loss")]
    if lp is not None:
        LOSS = base.LOSS = lp / 100
    prepared, gidx, rhos = setup()
    if stage == "split":
        rows = split_rows()
        for init in c2.INITS:
            print_rows(rows, run(f"split_loss5_{init}", rows, prepared, gidx, rhos[init]), init)
    elif stage == "split2":
        rows = split2_rows()
        for init in c2.INITS:
            print_rows(rows, run(f"split2_loss5_{init}", rows, prepared, gidx, rhos[init]), init)
    elif stage == "acceptance":
        reg, p0, p1 = argv[1], int(argv[2]), int(argv[3])
        if reg not in c2.REGIONS:
            reg = tuple(float(x) for x in reg.split(","))
        name = reg if isinstance(reg, str) else f"d{reg[0]:g}_f{reg[1]:g}".replace(".", "p")
        for init in inits:
            rows = acc_rows(reg, p0, p1)
            for k in range(0, len(rows), 22):          # two checkpoints of 22 points per start
                run(f"acc_loss5_{name}_{p0}_{p1}{tag}_{init}_{k:03d}", rows[k:k + 22], prepared, gidx, rhos[init])
    elif stage == "efield":
        E = float(argv[1])
        log(f"E = {E:g} V/cm: F'=3 main line shifted by {main_line_shift_MHz(E):+.3f} MHz from 200 V/cm; "
            f"Delta_main {info['delta_main_MHz']:.3f} MHz")
        for init in inits:
            run(f"efield_loss5_{E:g}_{init}", efield_rows(E), prepared, gidx, rhos[init])
    elif stage == "loss":
        for init in inits:
            run(f"loss_{argv[1]}_{init}", [(NEW, 35, 25, v, y) for v in PROXY_VX for y in PROXY_Y0],
                prepared, gidx, rhos[init])
    elif stage == "power":
        for init in inits:
            run(f"power_loss5_{init}", [(NEW, a, b, v, y) for a, b in MORE_POWER for v in PROXY_VX for y in PROXY_Y0],
                prepared, gidx, rhos[init])
    elif stage == "corrdet":        # measured-scan analogue: AOM fixed, laser scanned, triplet start, F=3 left (R(2) F'=4)
        # corrdet 209 20 30 -90 150 4 e217 vx0   (vx0: on axis; vxpm: vx = -0.6, +0.6 m/s; proxy: the 6 proxy points)
        from scripts.cleanup_spb_doppler import describe
        from centrex_tlf import lindblad
        f_aom, p0, p1 = float(argv[1]), int(argv[2]), int(argv[3])
        lo, hi, step = (float(a) for a in argv[4:7])
        vxs = [0.0] if "vx0" in argv else (list(PROXY_VX) if "proxy" in argv else [-0.6, 0.6])
        ys = list(PROXY_Y0) if "proxy" in argv else [0.0]
        # "sing": start in the singlet F1=5/2 F=2 mF=0 instead (R(2) F1'=7/2 F'=3); new runs also keep every J=2 level
        d0 = np.round(np.arange(lo, hi + step / 2, step), 6)
        desc = {int(i): describe(c2.OBE.QN[i]) for i in gidx}
        lab = {i: (float(desc[i]["F1"]), float(desc[i]["F"]), int(round(desc[i]["mF"]))) for i in gidx}
        levels = sorted(gidx, key=lambda i: lab[i])
        start_tag = "sing" if "sing" in argv else "trip"
        Fstart = 2 if start_tag == "sing" else 3
        trip = [i for i in gidx if desc[i]["p_mj"].get(0, 0) > 0.5 and int(desc[i]["F"]) == Fstart]
        f3 = [i for i in gidx if int(desc[i]["F"]) == 3 and float(desc[i]["F1"]) == 2.5]
        assert len(trip) == (3 if start_tag == "trip" else 1) and len(f3) == 7
        dim = len(c2.OBE.QN)
        rho0 = np.zeros((dim, dim), dtype=np.complex128)
        for i in trip:
            rho0[i, i] = 1 / len(trip)
        solver = {**SCAN_SOLVER, "maxiters": 1_000_000}
        vtag = "vx0" if vxs == [0.0] else ("proxy" if "proxy" in argv else "vxpm")
        stem = (f"corrdet_loss{100 * LOSS:g}_f{f_aom:g}_{p0}_{p1}{tag}_{start_tag}_{vtag}_d{lo:g}_{hi:g}_{step:g}"
                .replace(".", "p"))
        beta = next((float(a[4:]) for a in argv[1:] if a.startswith("beta")), None)    # e.g. beta0: no EOM
        if beta is not None:
            stem += f"_beta{beta:g}".replace(".", "p")
        rows = [((d, f_aom), p0, p1, v, y) + (() if beta is None else (beta,)) for v in vxs for y in ys for d in d0]
        out = {"F3": [], "J2": [], "L": []}
        for k in range(0, len(rows), CHUNK):
            f = OUT / f"{stem}_{k:03d}.npz"
            if not f.exists():
                t0 = time.perf_counter()
                res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=SLOTS,
                                              parameter_batch=np.asarray([slot_values(r) for r in rows[k:k + CHUNK]],
                                                                         dtype=np.complex128),
                                              threads=THREADS, **solver)
                pops = np.asarray(res.values).real
                np.savez(f, rows=np.array([repr(r) for r in rows[k:k + CHUNK]]), F3=pops[:, f3].sum(1),
                         J2=pops[:, gidx].sum(1), L=pops[:, levels], wall=time.perf_counter() - t0)
                log(f"{f.stem}: {len(rows[k:k + CHUNK])} trajectories in {time.perf_counter() - t0:.0f} s")
            z = np.load(f)
            out["F3"] += list(z["F3"]); out["J2"] += list(z["J2"])
            out["L"] = out["L"] + list(z["L"]) if ("L" in z.files and out["L"] is not None) else None
        extra = {} if out["L"] is None else {
            "levels": np.array(out["L"]).reshape(len(vxs), len(ys), len(d0), len(levels)).mean(1),
            "F1": [lab[i][0] for i in levels], "F": [lab[i][1] for i in levels], "mF": [lab[i][2] for i in levels]}
        np.savez(OUT / f"{stem}_result.npz", delta0_MHz=d0, vx=np.array(vxs), y0=np.array(ys), start=start_tag,
                 F3=np.array(out["F3"]).reshape(len(vxs), len(ys), len(d0)).mean(1),
                 J2=np.array(out["J2"]).reshape(len(vxs), len(ys), len(d0)).mean(1), **extra,
                 f_AOM_MHz=f_aom, powers_mW=(p0, p1), E_V_cm=float(base.E_FIELD[2]), loss=LOSS,
                 delta_main_MHz=info["delta_main_MHz"])
        log(f"{stem}: done, {len(rows)} trajectories")
    elif stage == "eomfreq":       # EOM frequency at one set point, proxy: eomfreq 14.5,209 32 43 e217
        pt, p0, p1 = parse_pt(argv[1]), int(argv[2]), int(argv[3])
        for init in inits:
            rem = chunked(f"eomfreq_loss{100 * LOSS:g}_{pt_name(pt)}_{p0}_{p1}{tag}_{init}",
                          [(pt, p0, p1, v, y, base.BETA, base.POL_FREQUENCY / (2 * np.pi * 1e6), 0.0, f)
                           for f in EOM_FREQS for v in PROXY_VX for y in PROXY_Y0], prepared, gidx, rhos[init])
            log(f"eomfreq {init} at {pt}: proxy " + ", ".join(f"{f:g} MHz {100 * (1 - x):.2f} %"
                                                           for f, x in zip(EOM_FREQS, proxy(rem))))
    elif stage == "polfreq":
        for init in inits:
            chunked(f"polfreq_loss5_35_25_{init}", [(NEW, 35, 25, v, y, base.BETA, f) for f in POL_FREQS
                                                     for v in PROXY_VX for y in PROXY_Y0], prepared, gidx, rhos[init])
    elif stage == "polphase":
        for init in inits:
            chunked(f"polphase_loss5_35_25_{init}", [(NEW, 35, 25, v, y, base.BETA, f, ph) for f in POL_PHASE_FREQS
                                                      for ph in POL_PHASES for v in PROXY_VX for y in PROXY_Y0],
                    prepared, gidx, rhos[init])
    elif stage == "vz":
        log(f"vz = {base.VZ:g} m/s, T_END = {T_END * 1e6:.1f} us")
        for init in inits:
            run(f"vz_loss5_{base.VZ:g}_{init}", [(NEW, 35, 25, v, y) for v in PROXY_VX for y in PROXY_Y0],
                prepared, gidx, rhos[init])
    elif stage == "yfill":
        for init in inits:
            run(f"yfill_loss5_d12p5_f207p5_35_25{tag}_{init}", [(NEW, 35, 25, v, 7.5e-3) for v in VX],
                prepared, gidx, rhos[init])
    elif stage == "splitproxy":
        for init in inits:
            chunked(f"splitproxy_loss5_{init}", splitproxy_rows(), prepared, gidx, rhos[init])
    elif stage == "confirm":
        for init in inits:
            run(f"confirm_loss5_{init}", confirm_rows(), prepared, gidx, rhos[init], solver=CONFIRM_SOLVER)
    elif stage == "beta":
        for init in inits:
            chunked(f"beta_loss5_35_25_{init}", beta_rows(), prepared, gidx, rhos[init])
    elif stage == "acc24":                                                  # candidates, both starts
        p0, p1 = int(argv[1]), int(argv[2])
        pts = [parse_pt(a) for a in argv[3:]]
        for init in inits:
            chunked(f"acc24_loss{100 * LOSS:g}_{p0}_{p1}{tag}_{init}_" + "+".join(pt_name(p) for p in pts),
                    [(pt, p0, p1, v, y) for pt in pts for v in VX24 for y in Y0], prepared, gidx, rhos[init])
    elif stage == "setpoints":
        p0, p1 = int(argv[1]), int(argv[2])
        e0 = next((float(a[4:]) for a in argv[3:] if a.startswith("from")), None)
        pts = [parse_pt(a) for a in argv[3:] if not a.startswith("from")]
        if e0 is not None:          # laser and AOM fixed while the field drifts from e0: delta0 -> delta0 - shift
            shift = main_line_shift_MHz(float(base.E_FIELD[2])) - main_line_shift_MHz(e0)
            log(f"field {base.E_FIELD[2]:g} V/cm, set points of {e0:g} V/cm: F'=3 main line moved {shift:+.3f} MHz")
            tag += f"_from{e0:g}".replace(".", "p")
            rows_pts = [(d0 - shift, fa) for d0, fa in pts]
        else:
            rows_pts = pts
        for init in inits:
            chunked(f"setpoints_loss{100 * LOSS:g}_{p0}_{p1}{tag}_{init}_" + "+".join(pt_name(p) for p in pts),
                    [(pt, p0, p1, v, y) for pt in rows_pts for v in PROXY_VX for y in PROXY_Y0], prepared, gidx, rhos[init])
    elif stage == "perstate":
        pt, p0, p1 = parse_pt(argv[1]), int(argv[2]), int(argv[3])
        perstate(f"perstate_loss{100 * LOSS:g}_{pt_name(pt)}_{p0}_{p1}{tag}", pt, p0, p1, prepared, gidx)
    elif stage == "map":            # vx x y0 map (0.2 m/s x 3 mm, y0 0...15 mm) at one set point, reusing acc24 points
        pt, p0, p1 = parse_pt(argv[1]), int(argv[2]), int(argv[3])
        ym = np.round(np.linspace(0.0, 15e-3, 6), 9)
        solver = {**SCAN_SOLVER, "maxiters": 1_000_000}
        for init in inits:
            known = {}
            for f in sorted(OUT.glob(f"acc24_loss{100 * LOSS:g}_{p0}_{p1}{tag}_{init}_*.npz")):
                z = np.load(f, allow_pickle=True)
                for row, r in zip(z["rows"], z["remaining"]):
                    q = eval(str(row), {"np": np})
                    if tuple(float(x) for x in q[0]) == pt:
                        known[(round(float(q[3]), 6), round(float(q[4]), 9))] = float(r)
            todo = [(pt, p0, p1, float(v), float(y)) for v in VX for y in ym
                    if (round(float(v), 6), round(float(y), 9)) not in known]
            stem = f"map_loss{100 * LOSS:g}_{pt_name(pt)}_{p0}_{p1}{tag}_{init}"
            log(f"{stem}: {len(known)} grid points from acc24, {len(todo)} to run")
            rem = np.concatenate([run(f"{stem}_{k:03d}", todo[k:k + CHUNK], prepared, gidx, rhos[init], solver)
                                  for k in range(0, len(todo), CHUNK)]) if todo else []
            known.update({(round(r[3], 6), round(r[4], 9)): float(x) for r, x in zip(todo, rem)})
            grid = np.array([[known[(round(float(v), 6), round(float(y), 9))] for y in ym] for v in VX])
            np.savez(OUT / f"{stem}_grid.npz", vx=VX, y0=ym, remaining=grid, set_point_MHz=pt, powers_mW=(p0, p1),
                     E_V_cm=float(base.E_FIELD[2]), loss=LOSS, init=init)
            log(f"{stem}: on axis {100 * (1 - grid[len(VX) // 2, 0]):.3f} %")
    elif stage == "cu2final":       # final population of every J=2 level, from each lens (m_J ~ 0) level alone
        from scripts.cleanup_spb_doppler import describe
        pt, p0, p1 = parse_pt(argv[1]), int(argv[2]), int(argv[3])
        desc = {int(i): describe(c2.OBE.QN[i]) for i in gidx}
        lab = {i: (float(d["F1"]), float(d["F"]), int(round(d["mF"]))) for i, d in desc.items()}
        starts = [i for i in gidx if desc[i]["p_mj"].get(0, 0) > 0.5 and lab[i][2] >= 0]   # mF < 0 by reflection
        jobs = [(i, (pt, p0, p1, v, y)) for i in starts for v in VX24 for y in Y0]
        stem = f"cu2final_loss{100 * LOSS:g}_{pt_name(pt)}_{p0}_{p1}{tag}"
        pops = np.concatenate([run_states(f"{stem}_{k:03d}", jobs[k:k + CHUNK], prepared, gidx, full=True)
                               for k in range(0, len(jobs), CHUNK)])
        levels = sorted(gidx, key=lambda i: lab[i])
        fin = pops[:, levels].reshape(len(starts), len(VX24), len(Y0), len(levels))
        np.savez(OUT / f"{stem}_result.npz", start_F1=[lab[i][0] for i in starts], start_F=[lab[i][1] for i in starts],
                 start_mF=[lab[i][2] for i in starts], F1=[lab[i][0] for i in levels], F=[lab[i][1] for i in levels],
                 mF=[lab[i][2] for i in levels], p_mJ0=[desc[i]["p_mj"].get(0, 0) for i in levels],
                 vx=VX24, y0=Y0, final=fin, set_point_MHz=pt, powers_mW=(p0, p1), E_V_cm=float(base.E_FIELD[2]),
                 loss=LOSS)
        f3 = [k for k, i in enumerate(levels) if lab[i][1] == 3]
        for s, i in enumerate(starts):
            log(f"{stem}: start {lab[i]}: on axis (vx 0.2) J=2 left {100 * fin[s, len(VX24) // 2, 0].sum():.3f} %, "
                f"F=3 left {100 * fin[s, len(VX24) // 2, 0, f3].sum():.3f} %")
    elif stage == "transfer":       # every J=2 level (mF >= 0) alone at one set point, proxy points, all final levels
        # transfer 13,209 20 30 e217 [acc]  -> final[start, vx, y0, level]; mF < 0 starts follow by reflection
        from scripts.cleanup_spb_doppler import describe
        pt, p0, p1 = parse_pt(argv[1]), int(argv[2]), int(argv[3])
        desc = {int(i): describe(c2.OBE.QN[i]) for i in gidx}
        lab = {i: (float(d["F1"]), float(d["F"]), int(round(d["mF"]))) for i, d in desc.items()}
        levels = sorted(gidx, key=lambda i: lab[i])
        acc = "acc" in argv          # acc: 24-point acceptance grid, only the starts the on/off comparison uses
        tvx, ty0 = (VX24, Y0) if acc else (np.array(PROXY_VX), np.array(PROXY_Y0))
        starts = [i for i in levels if lab[i][2] >= 0 and (not acc or lab[i][0] == 1.5
                  or (lab[i][1] == 3 and lab[i][2] <= 1) or (lab[i][1] == 2 and lab[i][2] == 0))]
        jobs = [(i, (pt, p0, p1, v, y)) for i in starts for v in tvx for y in ty0]
        stem = f"transfer_loss{100 * LOSS:g}_{pt_name(pt)}_{p0}_{p1}{tag}" + ("_acc" if acc else "")
        pops = np.concatenate([run_states(f"{stem}_{k:03d}", jobs[k:k + CHUNK], prepared, gidx, full=True)
                               for k in range(0, len(jobs), CHUNK)])
        fin = pops[:, levels].reshape(len(starts), len(tvx), len(ty0), len(levels))
        np.savez(OUT / f"{stem}_result.npz", start_F1=[lab[i][0] for i in starts], start_F=[lab[i][1] for i in starts],
                 start_mF=[lab[i][2] for i in starts], F1=[lab[i][0] for i in levels], F=[lab[i][1] for i in levels],
                 mF=[lab[i][2] for i in levels], p_mJ0=[desc[i]["p_mj"].get(0, 0) for i in levels],
                 vx=tvx, y0=ty0, final=fin, set_point_MHz=pt, powers_mW=(p0, p1),
                 E_V_cm=float(base.E_FIELD[2]), loss=LOSS)
        for s, i in enumerate(starts):
            log(f"{stem}: start {lab[i]}: vx {tvx[len(tvx) // 2]:g}, y0 0: J=2 left {100 * fin[s, len(tvx) // 2, 0].sum():.3f} %")
    elif stage == "corr":           # AOM fixed, laser scanned: delta0 and delta1 move together, on axis
        f_aom, p0, p1 = float(argv[1]), int(argv[2]), int(argv[3])
        lo, hi, step = (float(a) for a in argv[4:7])
        d0 = np.round(np.arange(lo, hi + step / 2, step), 6)
        log(f"correlated scan, f_AOM {f_aom:g} MHz: delta1 = delta0 {f_aom - info['delta_main_MHz']:+.3f} MHz")
        solver = {**SCAN_SOLVER, "maxiters": 1_000_000}
        for init in inits:
            stem = f"corr_loss{100 * LOSS:g}_f{f_aom:g}_{p0}_{p1}{tag}_{init}_d{lo:g}_{hi:g}_{step:g}".replace(".", "p")
            rows = [((d, f_aom), p0, p1, 0.0, 0.0) for d in d0]
            rem = np.concatenate([run(f"{stem}_{k:03d}", rows[k:k + CHUNK], prepared, gidx, rhos[init], solver)
                                  for k in range(0, len(rows), CHUNK)])
            np.savez(OUT / f"{stem}_result.npz", delta0_MHz=d0, delta1_MHz=d0 + f_aom - info["delta_main_MHz"],
                     remaining=rem, f_AOM_MHz=f_aom, powers_mW=(p0, p1), E_V_cm=float(base.E_FIELD[2]), loss=LOSS,
                     init=init)
            k = int(np.argmin(rem))
            log(f"{stem}: best {100 * (1 - rem[k]):.3f} % at delta0 {d0[k]:g} MHz")
    elif stage in GRIDS:
        p0, p1 = int(argv[1]), int(argv[2])
        for init in inits:
            chunked(f"{stage}_loss{100 * LOSS:g}_{p0}_{p1}{tag}_{init}", opt_rows(stage, p0, p1), prepared, gidx, rhos[init])
    summary()


def summary():
    from scripts.plot_cleanup_maps_current import SY, acceptance_mean
    for init in c2.INITS:
        f = OUT / f"split_loss5_{init}.npz"
        if f.exists():
            print_rows(split_rows(), np.load(f, allow_pickle=True)["remaining"], init)
        f = OUT / f"split2_loss5_{init}.npz"
        if f.exists():
            print_rows(split2_rows(), np.load(f, allow_pickle=True)["remaining"], init)
    for f in sorted(OUT.glob("acc_loss5_*_000.npz")):
        stem = f.stem[:-4]
        parts = sorted(OUT.glob(f"{stem}_*.npz"))
        rem = np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in parts])
        if rem.size < len(VX) * len(Y0):
            log(f"{stem}: {rem.size} of {len(VX) * len(Y0)} points"); continue
        g = rem.reshape(len(VX), len(Y0))
        log(f"{stem}: on axis {100 * (1 - g[len(VX) // 2, 0]):.2f} %, "
            f"mean in acceptance {100 * (1 - acceptance_mean(g, VX, Y0, SY['CU2'])):.2f} %")
    for f in sorted(OUT.glob("acc24_*_000.npz")):
        stem = f.stem[:-4]
        rem = np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in sorted(OUT.glob(f"{stem}_[0-9]*.npz"))])
        names = stem.split("_spb_" if "_spb_" in stem else "_uniform_")[1].split("+")
        n = len(VX24) * len(Y0)
        for k, name in enumerate(names):
            if rem.size >= (k + 1) * n:
                g = rem[k * n:(k + 1) * n].reshape(len(VX24), len(Y0))
                log(f"{stem.split('_d')[0]} {name}: mean in acceptance (reduced grid) "
                    f"{100 * (1 - acceptance_mean(g, VX24, Y0, SY['CU2'])):.2f} %")
    for f in sorted(OUT.glob("setpoints_*_000.npz")):
        stem = f.stem[:-4]
        rem = np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in sorted(OUT.glob(f"{stem}_[0-9]*.npz"))])
        names = stem.split("_spb_" if "_spb_" in stem else "_uniform_")[1].split("+")
        pr = 100 * (1 - proxy(rem[:rem.size // 6 * 6]))
        log(f"{stem.split('_d')[0]}: proxy " + ", ".join(f"{nm} {x:.2f} %" for nm, x in zip(names, pr)))
    for f in sorted(f for g in GRIDS for f in OUT.glob(f"{g}_loss*_000.npz")):
        stem = f.stem[:-4]
        d0s, aoms = GRIDS[stem.split("_")[0]]
        parts = sorted(OUT.glob(f"{stem}_*.npz"))
        rem = np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in parts])
        n = rem.size // (len(PROXY_VX) * len(PROXY_Y0))
        if n == 0:
            continue
        pr = proxy(rem[:n * len(PROXY_VX) * len(PROXY_Y0)])
        log(f"{stem}: proxy efficiency (%), {n} of {len(d0s) * len(aoms)} settings; rows delta0, cols f_AOM {aoms}")
        tab = np.full(len(d0s) * len(aoms), np.nan); tab[:n] = 100 * (1 - pr)
        for d0, r in zip(d0s, tab.reshape(len(d0s), len(aoms))):
            log(f"   {d0:5.1f}: " + "  ".join(f"{x:6.2f}" for x in r))

    for init in c2.INITS:
        parts = sorted(OUT.glob(f"beta_loss5_35_25_{init}_*.npz"))
        if not parts:
            continue
        rem = np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in parts])
        n = rem.size // (len(PROXY_VX) * len(PROXY_Y0))
        tab = np.full(len(BETA_POINTS) * len(BETAS), np.nan)
        tab[:n] = 100 * (1 - proxy(rem[:n * len(PROXY_VX) * len(PROXY_Y0)]))
        log(f"beta_loss5_35_25_{init}: proxy efficiency (%), rows (delta0, f_AOM), cols beta {BETAS}")
        for pt, r in zip(BETA_POINTS, tab.reshape(len(BETA_POINTS), len(BETAS))):
            log(f"   {pt}: " + "  ".join(f"{x:6.2f}" for x in r))

    for init in c2.INITS:
        parts = sorted(OUT.glob(f"polfreq_loss5_35_25_{init}_*.npz"))
        if parts:
            rem = np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in parts])
            n = rem.size // (len(PROXY_VX) * len(PROXY_Y0))
            pr = 100 * (1 - proxy(rem[:n * len(PROXY_VX) * len(PROXY_Y0)]))
            log(f"polfreq_loss5_35_25_{init} at {NEW}: " + ", ".join(f"{f:g} MHz {x:.2f} %"
                                                                    for f, x in zip(POL_FREQS, pr)))

    for init in c2.INITS:
        parts = sorted(OUT.glob(f"polphase_loss5_35_25_{init}_*.npz"))
        p0 = sorted(OUT.glob(f"polfreq_loss5_35_25_{init}_*.npz"))
        if parts and p0:
            rem = np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in parts])
            n = rem.size // (len(PROXY_VX) * len(PROXY_Y0))
            if n == len(POL_PHASE_FREQS) * len(POL_PHASES):
                tab = 100 * (1 - proxy(rem)).reshape(len(POL_PHASE_FREQS), len(POL_PHASES))
                ph0 = 100 * (1 - proxy(np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in p0])))
                log(f"polphase_loss5_35_25_{init}: proxy (%), rows MHz, cols phase/pi [0] + {POL_PHASES}, mean")
                for f, r in zip(POL_PHASE_FREQS, tab):
                    row = [ph0[POL_FREQS.index(f)], *r]
                    log(f"   {f:4g}: " + "  ".join(f"{x:6.2f}" for x in row) + f"   mean {np.mean(row):6.2f}")

    for init in c2.INITS:
        parts = sorted(OUT.glob(f"splitproxy_loss5_{init}_*.npz"))
        if parts:
            rem = np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in parts])
            n = rem.size // (len(PROXY_VX) * len(PROXY_Y0))
            pr = 100 * (1 - proxy(rem[:n * len(PROXY_VX) * len(PROXY_Y0)]))
            log(f"splitproxy_loss5_{init} at {NEW}: " + ", ".join(f"{a}+{b} mW {x:.2f} %"
                                                               for (a, b), x in zip(PROXY_SPLITS, pr)))
        f = OUT / f"confirm_loss5_{init}.npz"
        if f.exists():
            for (pt, v), r in zip(CONFIRM_POINTS, np.load(f, allow_pickle=True)["remaining"]):
                log(f"confirm (tight tolerances) {init} {pt} vx {v}: {100 * (1 - r):.3f} %")

    for init in c2.INITS:
        for f in sorted(OUT.glob(f"loss_*_{init}.npz")):
            log(f"{f.stem} %: proxy {100 * (1 - proxy(np.load(f, allow_pickle=True)['remaining'])[0]):.2f} %")
        f = OUT / f"power_loss5_{init}.npz"
        if f.exists():
            pr = 100 * (1 - proxy(np.load(f, allow_pickle=True)["remaining"]))
            log(f"power_loss5_{init} at {NEW}: " + ", ".join(f"{a}+{b} mW {x:.2f} %" for (a, b), x in zip(MORE_POWER, pr)))
        for f in sorted(OUT.glob(f"vz_loss5_*_{init}.npz")):
            log(f"{f.stem}: proxy {100 * (1 - proxy(np.load(f, allow_pickle=True)['remaining'])[0]):.2f} %")
        for f in sorted(OUT.glob(f"efield_loss5_*_{init}.npz")):
            pr = 100 * (1 - proxy(np.load(f, allow_pickle=True)["remaining"]))
            log(f"{f.stem}: proxy, laser fixed {pr[0]:.2f} %, re-tuned {pr[1]:.2f} %")


if __name__ == "__main__":
    if sys.argv[1] == "summary":
        summary()
    else:
        main(sys.argv[1:])
