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
    python scripts/cu2_power_check.py polfreq       # X/Z switching frequency at NEW, 35 + 25 mW, proxy
    python scripts/cu2_power_check.py polphase      # ... and its phase at the molecule's arrival (random in the lab)
    python scripts/cu2_power_check.py summary

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
from scripts.cu2_laser_aom_scan import CONFIRM_SOLVER, SCAN_SOLVER, T_END, log  # noqa: E402

LOSS = 0.05
base.LOSS = LOSS                                            # read by base.build() when the model is prepared
OUT = ROOT / "tmp" / "cu2_power"
SLOTS = ["detuning0", "detuning1", "vx", "y0", "power0", "power1", "beta0", "beta1", "pol_frequency", "pol_phase"]
TOTAL = 60
SPLITS = [(TOTAL - p1, p1) for p1 in (10, 15, 20, 25, 30, 35, 40)]   # mW, (F'=3, F'=4)
SPLITS2 = [(15, 45), (10, 50)]                                     # extension: uniform start still rising
REFERENCE = [(80, 24)]
VX = np.round(np.arange(-1.0, 1.0 + 1e-9, 0.2), 6)
Y0 = np.array([0.0, 3e-3, 6e-3, 9e-3])
THREADS = 16
GRIDS = {                                                          # (delta0/2pi, f_AOM) in MHz
    "opt": ([5.0, 8.0, 11.5, 15.0, 18.5], [186.0, 191.0, 197.0, 203.0, 208.5, 214.0]),   # AOM range 200 +- 15 MHz
    "fine": ([10.0, 11.0, 11.5, 12.0, 13.0], [206.0, 207.5, 208.5, 209.5, 211.0]),      # around R1
    "fine2": ([12.5, 13.0, 13.5, 14.0, 15.0], [203.0, 204.5, 206.0, 207.5]),            # fine optimum at its edge
    "near": ([11.5, 12.0, 12.5, 13.0], [206.0, 207.5, 208.5]),                          # near NEW, other beams
}
BETA_POINTS = [(12.5, 207.5), (13.0, 206.0)]                       # best proxy points of Steps 4-5
BETAS = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5]                           # one EOM for both beams: beta0 = beta1
NEW = (12.5, 207.5)                                                # best full-grid point for the SPB start (Step 6)
PROXY_SPLITS = [(25, 35), (30, 30), (35, 25), (40, 20), (45, 15)]
CONFIRM_POINTS = [(NEW, 0.0), (NEW, 0.6), ((13.0, 206.0), 0.0), ((13.0, 206.0), 0.6)]   # (set point, vx), y0 = 0
MORE_POWER = [(47, 33), (58, 42), (80, 24)]                        # 80, 100 mW at 35:25, and the old 80 + 24 mW
POL_FREQS = [0.1, 0.25, 0.5, 1.0, 2.0]                             # MHz, X/Z square-wave switching
POL_PHASE_FREQS, POL_PHASES = [0.25, 0.5, 1.0], [0.5, 1.0, 1.5]    # MHz, units of pi (phase 0 is stage polfreq)
PROXY_VX, PROXY_Y0 = [-0.6, 0.0, 0.6], [0.0, 6e-3]
CHUNK = 16
info = None


def run(name, rows, prepared, gidx, rho0, solver=SCAN_SOLVER):
    """rows: (region or (delta0, f_AOM) in MHz, P0 mW, P1 mW, vx, y0[, beta[, pol MHz[, phase/pi]]]) -> remaining J=2
    population; one checkpoint per call."""
    from centrex_tlf import lindblad
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / f"{name}.npz"
    if f.exists():
        return np.load(f, allow_pickle=True)["remaining"]
    batch = []
    for reg, p0, p1, vx, y0, *rest in rows:
        d0, faom = c2.REGIONS[reg] if isinstance(reg, str) else reg
        d1 = d0 + faom - info["delta_main_MHz"]
        beta = rest[0] if rest else base.BETA
        pol = 2 * np.pi * 1e6 * rest[1] if len(rest) > 1 else base.POL_FREQUENCY
        phase = np.pi * rest[2] if len(rest) > 2 else 0.0
        batch.append([2 * np.pi * 1e6 * d0, 2 * np.pi * 1e6 * d1, vx, y0, p0 * 1e-3, p1 * 1e-3, beta, beta, pol,
                      phase])
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
    global LOSS, T_END
    if stage == "loss":
        LOSS = base.LOSS = float(argv[1]) / 100                              # read by base.build()
    if stage == "vz":
        base.VZ = float(argv[1])                                             # read by build() and doppler_sign()
        T_END = base.T_END = 0.05 / base.VZ
    sy = next((float(a[2:]) for a in argv if a.startswith("sy")), None)      # beam sigma_y in mm
    argv = [a for a in argv if not a.startswith("sy")]
    if sy is not None:
        c2.SIGMA_Y_CU2 = sy * 1e-3                                           # read by c2.setup()
    tag = "" if sy is None else f"_sy{sy:g}".replace(".", "p")
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
