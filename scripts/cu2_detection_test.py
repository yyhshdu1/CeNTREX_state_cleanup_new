"""Expected CU2 on/off ratio of the R(2) F1'=7/2 F'=4 LIF signal in the first CU2 test (CU1 and SPB off).

    python scripts/cu2_detection_test.py yields [line]       # probe OBE: photons per molecule per J=2 level (~10 s)
    python scripts/cu2_detection_test.py plot 12.5,200 32 43 e155    # slide schematic of the on/off simulation
    python scripts/cu2_detection_test.py ratio 12.5,200 32 43 e155 [line]   # combine with the `cu2final` stage of
                                                                     # scripts/cu2_power_check.py at that set point
    line: F4 (default, R(2) F1'=7/2 F'=4, the test line), F1p5F2 / F1p5F3 (R(2) F1'=5/2 F'=2 / F'=3, which also
    see the F1 = 3/2 levels), R7o2F3 (R(2) F1'=7/2 F'=3), Q3o2F1 / Q3o2F2 / Q5o2F2 / Q5o2F3 (Q(2) F1'=3/2 F'=1, 2
    and F1'=5/2 F'=2, 3)

Probe (user, 2026-11-03): locked to the R(2) F1'=7/2 F'=4 main line, no EOM, sigma+ circular polarization
(DET2; the DET field plates are off for the CU2 test, so the axis is the earth field, perpendicular to k), flat-top 1.5 cm (along the molecular beam) x 1.5 cm (height), 15 mW, single pass,
same direction as the CU2 beams (Doppler shift k*vx), after CU2. E = 0 in the detection region, earth field taken as
0.5 G along the laser z axis (vertical, the CU2 field direction), every fluorescence photon collected, no background.

Default treatment (2026-11-04, `pol k`): quantization axis along k, where the light is pure sigma+, with the CU2
populations rotated onto k (`level_photons`). centrex's decay operators are pairwise (no coherence transfer), which is
adequate only with the axis along the light's symmetry axis; `yields_coherent` (`pol spc`) is the exact check with the
earth field perpendicular to k and coherent decay. The older `pol sp` numbers (vertical axis) are wrong.

    python scripts/cu2_detection_test.py yields Q5o2F3 polspk      # centrex probe run, sigma+ with the axis along k
    python scripts/cu2_detection_test.py yields_coherent Q5o2F3    # exact check (numpy, ~20 s per line)

Assumptions: the field ramp from CU2 to zero is adiabatic, so each CU2 eigenstate arrives as its zero-field parent
(F1, F, mF); a molecule keeps its height y between CU2 and the probe, so the probe sees |y| <= 7.5 mm; vx and y
distributions as for CU2 (sigma_vx 0.65 m/s cut at |vx| <= 1 m/s, sigma_y 10.35 mm). Coherences are not carried.
"""
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import vz_paths  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CACHE_DIR = ROOT / "tmp" / "detection"
LINES = {"F4": ("R2_F1_7o2_F4", dict(F1=2.5, F=3)), "F1p5F2": ("R2_F1_5o2_F2", dict(F1=1.5, F=2)),
         "F1p5F3": ("R2_F1_5o2_F3", dict(F1=1.5, F=2)), "R7o2F3": ("R2_F1_7o2_F3", dict(F1=2.5, F=2)),
         "Q3o2F1": ("Q2_F1_3o2_F1", dict(F1=1.5, F=1)), "Q3o2F2": ("Q2_F1_3o2_F2", dict(F1=1.5, F=2)),
         "Q5o2F2": ("Q2_F1_5o2_F2", dict(F1=2.5, F=2)), "Q5o2F3": ("Q2_F1_5o2_F3", dict(F1=2.5, F=3)),
         "P1o2F1": ("P2_F1_1o2_F1", dict(F1=1.5, F=2)), "P1o2F0": ("P2_F1_1o2_F0", dict(F1=1.5, F=1, mshift=-1)),
         "P3o2F1": ("P2_F1_3o2_F1", dict(F1=2.5, F=2)), "P3o2F2": ("P2_F1_3o2_F2", dict(F1=2.5, F=3)),
         "Q1o2F0": ("Q1_F1_1o2_F0", dict(J=1, F1=0.5, F=1, mshift=-1))}
# transition, main ground level (mF = 0): the J=2 level with the largest |<e|d.X|g>| (tmp/q2_me_table.py).
# P(2) F1'=1/2 lines only couple to F1=3/2. mshift moves the main pair (ground mF, excited mF') = (0, 1) by
# that amount: F'=0 has no mF'=1, so its main pair is (-1, 0). J (default 2) is the ground J of the line:
# Q(1) F1'=1/2 F'=0 is a J=1 line; its yields are then per J=1 level.


def ground_J(line):
    return LINES[line][1].get("J", 2)


def cache(line, pol=None):
    pol = pol or {"k": "spk", "kiso": "spk"}.get(POL, POL)
    sfx = "" if pol == "x" else f"_pol{pol}"
    return CACHE_DIR / (f"cu2_r2f4_yields{sfx}.json" if line == "F4" else f"cu2_r2_{line}_yields{sfx}.json")


def wigner_d2(F, beta=np.pi / 2):
    """|d^F_{m'm}(beta)|^2, rows m' and columns m both ordered -F ... +F (integer F)."""
    from scipy.linalg import expm
    m = np.arange(-F, F + 1)
    jp = np.zeros((2 * F + 1, 2 * F + 1))
    for k in range(2 * F):                    # J+ |m> = sqrt(F(F+1) - m(m+1)) |m+1>
        jp[k + 1, k] = np.sqrt(F * (F + 1) - m[k] * (m[k] + 1))
    jy = (jp - jp.T) / 2j
    return np.abs(expm(-1j * beta * jy)) ** 2


def level_photons(line, levels, vxs, power=15e-3):
    """Photons per molecule [vx, level] for J=2 levels (F1, F, mF) with mF about the vertical CU2 axis.

    POL "k" (default): quantization axis along the probe k, where the light is pure sigma+ and centrex's pairwise
    decay is adequate; each vertical-axis level is rotated onto k with |d^F_{m'm}(pi/2)|^2 and its photons are
    sum_m' |d|^2 n_k(m'). Exact if the probe dynamics are symmetric about k (coherences Delta m' != 0 created by the
    rotation then do not contribute); the 0.5 G earth field perpendicular to k is neglected during the probe.
    POL "kiso": alignment scrambled between CU2 and the probe (mF-averaged n_k). POL "spc": exact check, vertical
    axis, earth field 0.5 G vertical, coherence-preserving decay (numpy, `yields_coherent`). Other POL values read
    centrex runs in the vertical basis, which are wrong for light that is not pi/sigma about that axis.
    """
    d = json.loads(vz_paths.find(cache(line), CACHE_DIR).read_text())
    kP, lab = d["powers_W"].index(power), [tuple(x) for x in d["labels"]]
    ph = np.array(d["photons"])[kP]                                        # [vx grid, label]
    raw = np.array([[np.interp(v, d["vx_m_s"], ph[:, k]) for k in range(len(lab))] for v in vxs])
    n = np.zeros((len(vxs), len(levels)))
    for j, lv in enumerate(levels):
        if POL not in ("k", "kiso"):
            n[:, j] = raw[:, lab.index(lv)] if lv in lab else 0.0
            continue
        F1, F, mF = lv
        idx = [lab.index((F1, F, mp)) if (F1, F, mp) in lab else None for mp in range(-F, F + 1)]
        nk = np.array([raw[:, i] if i is not None else np.zeros(len(vxs)) for i in idx])   # [m', vx]
        w = np.full(2 * F + 1, 1 / (2 * F + 1)) if POL == "kiso" else wigner_d2(F)[:, mF + F]
        n[:, j] = w @ nk
    return n


OUT = ROOT / "tmp" / "cu2_power"

VZ = 156.0                                   # m/s, measured (TOF); 184 before 2026-11-09 (tmp/*/vz184/)
WIDTH, HEIGHT = 0.015, 0.015                 # flat top, along the beam and in height (m)
T_PROBE = WIDTH / VZ
POWERS = [7.5e-3, 15e-3, 24e-3, 30e-3]       # W; 24 mW measured in the CU2 test (2026-11-09), the others show saturation
VXS = np.round(np.arange(-1.0, 1.0 + 1e-9, 0.2), 6)
K_VX = 2 * np.pi * 3.68e6                    # k*vx per m/s at 271.7 nm, rad/s
B_EARTH = 0.5                                # G
THREADS = 4                                  # Rust batch threads of the probe OBE
SOLVER = "dopri5"                            # probe ODE solver; "tsit5" / "fixed_dopri5" if dopri5 reports
DT = 2e-9                                    # stiffness; DT = initial step (adaptive) or the step (fixed)
E_DET = 0.0                                  # V/cm, DET plates field along the probe z axis = k (0: plates off).
                                             # E != 0 also keeps the opposite-parity B levels (Stark mixed).
POL = "k"                                    # probe treatment, set with "pol<name>". The DET laser is sigma+ about
                                             # its direction k (y); the DET plates are off for the CU2 test.
                                             # "k" (default): axis along k, CU2 populations rotated onto k (see
                                             # level_photons); "kiso": same, alignment scrambled; "spc": exact
                                             # check (vertical axis, earth field 0.5 G vertical, coherent decay).
                                             # Probe runs: "spk" = sigma+ about the run's z axis with B along it
                                             # (used by "k"/"kiso"); "sp" = sigma+ about y in the vertical basis
                                             # with centrex's pairwise decay: WRONG (no coherence transfer in
                                             # decay), kept only for comparison; "x" / "z": linear, same caveat
                                             # for "x".
POL_VEC = {"x": [1.0, 0, 0], "z": [0, 0, 1.0], "sp": [1j / 2 ** 0.5, 0, 1 / 2 ** 0.5],
           "spk": [1 / 2 ** 0.5, -1j / 2 ** 0.5, 0], "spc": [1j / 2 ** 0.5, 0, 1 / 2 ** 0.5],
           "k": None, "kiso": None}             # centrex: x - iy drives Delta mF = +1
POLNAME = {"x": "X pol", "z": "Z pol", "sp": r"$\sigma^+$ (pairwise decay, wrong)", "spk": r"$\sigma^+$ along $B$",
           "spc": r"$\sigma^+$ about $k$, exact", "k": r"$\sigma^+$ about $k$", "kiso": r"$\sigma^+$ about $k$"}
SVX, VCUT, SY = 0.65, 1.0, 10.35e-3


def build(line="F4"):
    from centrex_tlf import couplings, lindblad, states, transitions
    from centrex_tlf.couplings.polarization import Polarization
    warnings.simplefilter("ignore")
    X, B = states.ElectronicState.X, states.ElectronicState.B
    tr = getattr(transitions, LINES[line][0])
    pol = Polarization(np.array(POL_VEC[POL], dtype=complex), name=POL.upper())
    qn = dict(LINES[line][1])
    shift = qn.pop("mshift", 0)
    J = qn.pop("J", 2)
    g = 1 * next(iter(states.generate_coupled_states_X(states.QuantumSelector(
        electronic=X, J=J, mF=(1 if POL == "z" else 0) + shift, P=tr.P_ground, **qn))))
    e = 1 * next(iter(states.generate_coupled_states_B(states.QuantumSelector(       # X light: Delta mF = +-1
        electronic=B, J=tr.J_excited, F1=tr.F1_excited, F=tr.F_excited, mF=1 + shift, P=tr.P_excited))))
    sel = couplings.generate_transition_selectors(transitions=[tr], polarizations=[[pol]], ground_mains=[g],
                                                  excited_mains=[e])
    obe = lindblad.generate_OBE_system_transitions([tr], sel, qn_compact=True, E=np.array([0.0, 0.0, E_DET]),
                                                   B=np.array([0.0, 0.0, B_EARTH]),
                                                   retain_opposite_parity_levels=E_DET != 0.0)
    p = lindblad.LindbladParameters()
    rabi, det = p.real("rabi", 0.0), p.real("det", 0.0)
    for s in obe.H_symbolic.free_symbols:
        if s in obe.coupling_symbols:
            p.bind(s, rabi, finalize=False)
        elif str(s) == str(sel[0].δ):
            p.bind(s, det, finalize=False)
        else:
            p.real(str(s), 0.0)
    for grp in obe.polarization_symbols:
        for s in (grp if isinstance(grp, (list, tuple)) else [grp]):
            p.bind(s, 1.0, finalize=False)
    p._finalize()
    prepared = lindblad.prepare_lindblad_problem(obe, p, backend="rust", hamiltonian_representation="decomposed")
    return prepared, obe


def yields(line="F4"):
    from centrex_tlf import hamiltonian, lindblad, states
    from centrex_tlf.utils.rabi import power_to_rabi_rectangular_beam
    gamma = float(getattr(hamiltonian, "Γ"))
    X, B = states.ElectronicState.X, states.ElectronicState.B
    t0 = time.time()
    prepared, obe = build(line)
    QN, dim = obe.QN, len(obe.QN)
    mc = float(abs(obe.couplings[0].main_coupling))
    g2 = [i for i, s in enumerate(QN) if s.largest.electronic_state == X and s.largest.J == ground_J(line)]
    exc = [i for i, s in enumerate(QN) if s.largest.electronic_state == B]
    labels = [[float(QN[i].largest.F1), int(QN[i].largest.F), int(QN[i].largest.mF)] for i in g2]
    print(f"probe OBE: {dim} levels, {len(g2)} J={ground_J(line)}, {len(exc)} excited, built in {time.time() - t0:.0f} s", flush=True)
    rhos, rows = [], []
    for P in POWERS:
        for vx in VXS:
            for i in g2:
                rho = np.zeros((dim, dim), dtype=complex); rho[i, i] = 1.0
                rhos.append(prepared.layout.pack(rho))
                rows.append([power_to_rabi_rectangular_beam(P, mc, WIDTH, HEIGHT), K_VX * vx])
    res = lindblad.solve_lindblad_batch(prepared, np.asarray(rhos), (0.0, T_PROBE), parameter_slots=["rabi", "det"],
                                        parameter_batch=np.asarray(rows, dtype=np.complex128), solver=SOLVER,
                                        execution_mode="expanded_sparse", dt=DT, reltol=1e-7, abstol=1e-9,
                                        maxiters=2_000_000, output="photon_integral",
                                        integral_weights=[(i, gamma) for i in exc], output_when="final", threads=THREADS)
    ph = np.asarray(res.values).real.reshape(len(POWERS), len(VXS), len(g2))
    out = dict(probe=dict(line=LINES[line][0], polarization=POL.upper(), width_m=WIDTH, height_m=HEIGHT, vz=VZ,
                          t_end_s=T_PROBE, field=f"E = {E_DET:g} V/cm, B = {B_EARTH} G along z", profile="flat top, single pass"),
               powers_W=POWERS, vx_m_s=VXS.tolist(), labels=labels, photons=ph.tolist(),
               levels=dim, wall_s=time.time() - t0)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache(line).write_text(json.dumps(out, indent=1))
    k15, k0 = (POWERS.index(15e-3) if 15e-3 in POWERS else 0), int(np.argmin(abs(VXS)))
    for lab, n in zip(labels, ph[k15, k0]):
        if n > 1e-3:
            print(f"  F1={lab[0]} F={lab[1]} mF={lab[2]:+d}: {n:.3f} photons at vx 0 (15 mW); "
                  f"vx +-1: {ph[k15, 0, labels.index(lab)]:.3f}, 7.5 / 30 mW: "
                  f"{ph[0, k0, labels.index(lab)]:.3f} / {ph[-1, k0, labels.index(lab)]:.3f}")
    print("wrote", cache(line), f"({time.time() - t0:.0f} s)")


def yields_coherent(line="F4"):
    """Exact probe check ("spc"): vertical axis, earth field 0.5 G vertical (perpendicular to k), sigma+ about k,
    and coherence-preserving decay operators C_q (q = -1, 0, +1) inside J=2, integrated with a dense matrix
    exponential (the probe is static: flat top, one frequency). centrex's own collapse operators are pairwise
    (one |g><e| per pair) and its prepared solvers reject anything else, so decay there never transfers
    coherences; that is only adequate with the axis along the light's symmetry axis."""
    import sympy as smp
    from scipy.linalg import expm
    from centrex_tlf import couplings, hamiltonian, states
    from centrex_tlf.utils.rabi import power_to_rabi_rectangular_beam
    global POL
    gamma = float(getattr(hamiltonian, "Γ"))
    X, B = states.ElectronicState.X, states.ElectronicState.B
    POL = "spc"
    t0 = time.time()
    _, obe = build(line)
    QN, dim = obe.QN, len(obe.QN)
    gi = [i for i, s in enumerate(QN) if s.largest.electronic_state == X and s.largest.J == 2]
    ei = [i for i, s in enumerate(QN) if s.largest.electronic_state == B]
    labels = [[float(QN[i].largest.F1), int(QN[i].largest.F), int(QN[i].largest.mF)] for i in gi]
    mc = float(abs(obe.couplings[0].main_coupling))
    rabi_s, det_s = smp.Symbol("rabi_"), smp.Symbol("det_")
    sub = {}
    for s in obe.H_symbolic.free_symbols:
        sub[s] = rabi_s if s in obe.coupling_symbols else (det_s if str(s).startswith("δ") else 0.0)
    for grp in obe.polarization_symbols:
        for s in (grp if isinstance(grp, (list, tuple)) else [grp]):
            sub[s] = 1.0
    assert sum(1 for s in obe.H_symbolic.free_symbols if str(s).startswith("δ")) == 1
    Hf = smp.lambdify((rabi_s, det_s), obe.H_symbolic.subs(sub), "numpy")
    C0 = np.asarray(obe.C_array, dtype=complex)
    keep = [C for C in C0 if not any(abs(C[i, j]) > 0 for i in gi for j in ei)]   # out of J=2: pairwise is fine
    target = np.array([sum(abs(C[i, j]) ** 2 for C in C0 for i in gi) for j in ei])
    sub_qn = [QN[i] for i in gi + ei]
    Ms = [couplings.generate_coupling_matrix(sub_qn, sub_qn[:len(gi)], sub_qn[len(gi):],
                                             pol_vec=np.asarray(v, dtype=complex) / np.linalg.norm(v),
                                             normalize_pol=False)[:len(gi), len(gi):]
          for v in ([1, -1j, 0], [1, 1j, 0], [0, 0, 1])]                       # q = +1, -1, 0
    scale = np.sqrt(target / sum(abs(M) ** 2 for M in Ms).sum(axis=0))
    Cs = keep
    for M in Ms:
        C = np.zeros((dim, dim), complex)
        C[np.ix_(gi, ei)] = M * scale[None, :]
        Cs.append(C)
    I = np.eye(dim)
    D = sum(np.kron(C, C.conj()) - 0.5 * np.kron(C.conj().T @ C, I) - 0.5 * np.kron(I, (C.conj().T @ C).T)
            for C in Cs)
    w = np.zeros(dim * dim)
    w[[e * dim + e for e in ei]] = gamma
    V0 = np.zeros((dim * dim + 1, len(gi)), complex)
    for k, i in enumerate(gi):
        V0[i * dim + i, k] = 1.0
    ph = np.zeros((len(POWERS), len(VXS), len(gi)))
    for a, P in enumerate(POWERS):
        rabi = power_to_rabi_rectangular_beam(P, mc, WIDTH, HEIGHT)
        for b, vx in enumerate(VXS):
            H = np.asarray(Hf(rabi, K_VX * vx), dtype=complex)
            L = -1j * (np.kron(H, I) - np.kron(I, H.T)) + D
            A = np.zeros((dim * dim + 1, dim * dim + 1), complex)
            A[:-1, :-1], A[-1, :-1] = L, w
            out = expm(A * T_PROBE) @ V0
            ph[a, b] = out[-1].real
            tr = np.array([np.trace(out[:-1, k].reshape(dim, dim)).real for k in range(len(gi))])
            assert np.allclose(tr, 1.0, atol=1e-6), tr
    out = dict(probe=dict(line=LINES[line][0], polarization="SPC (sigma+ about k, coherent decay, numpy)",
                          width_m=WIDTH, height_m=HEIGHT, vz=VZ, t_end_s=T_PROBE,
                          field=f"E = 0, B = {B_EARTH} G along z (vertical, perpendicular to k)",
                          profile="flat top, single pass"),
               powers_W=POWERS, vx_m_s=VXS.tolist(), labels=labels, photons=ph.tolist(), levels=dim,
               wall_s=time.time() - t0)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache(line).write_text(json.dumps(out, indent=1))
    print("wrote", cache(line), f"({time.time() - t0:.0f} s)")


class Model:
    """CU2 final populations (cu2final stage) for the lens starts, and the acceptance weighting of the probe."""

    def __init__(self, pt, p0, p1, tag):
        stem = f"cu2final_loss5_d{pt[0]:g}_f{pt[1]:g}".replace(".", "p") + f"_{p0}_{p1}{tag}"
        z = np.load(vz_paths.find(OUT / f"{stem}_result.npz"), allow_pickle=True)
        self.z, self.pt, self.p = z, pt, (p0, p1)
        self.levels = [(float(a), int(b), int(c)) for a, b, c in zip(z["F1"], z["F"], z["mF"])]
        starts = [(float(a), int(b), int(c)) for a, b, c in zip(z["start_F1"], z["start_F"], z["start_mF"])]
        self.vx, self.y0, fin = z["vx"], z["y0"], z["final"]                 # fin[start, vx, y0, level]
        flip = [self.levels.index((f1, f, -m)) for f1, f, m in self.levels]   # mF < 0 starts by reflection
        self.fin = {s: fin[k] for k, s in enumerate(starts)}
        for k, (f1, f, m) in enumerate(starts):
            if m > 0:
                self.fin[(f1, f, -m)] = fin[k][..., flip]
        self.wv = np.exp(-0.5 * (self.vx / SVX) ** 2) * (abs(self.vx) <= VCUT + 1e-12)
        self.wv /= self.wv.sum()
        self.ym = np.linspace(0.0, HEIGHT / 2, 76)
        self.wy = np.exp(-0.5 * (self.ym / SY) ** 2)
        self.wy[[0, -1]] *= 0.5
        self.wy /= self.wy.sum()
        trip = [s for s in self.fin if s[1] == 3]
        sing = [s for s in self.fin if s[1] == 2]
        self.weights = {"equal": {s: 1 / 4 for s in self.fin},
                        "2/3": {**{s: 2 / 3 for s in sing}, **{s: 1 / 9 for s in trip}}}

    def mean(self, a):
        """a[vx, y0] -> mean over the probed acceptance."""
        return self.wv @ np.array([np.interp(self.ym, self.y0, r) for r in a]) @ self.wy

    def photons(self, line):
        """Photons per molecule [vx, level] at 15 mW (see level_photons for the POL modes)."""
        return level_photons(line, self.levels, self.vx)

    def off(self, s):
        """CU2 off: the population stays in the start level."""
        return np.broadcast_to(np.eye(len(self.levels))[self.levels.index(s)], self.fin[s].shape)

    def signal(self, pop, n):
        """pop[vx, y0, level] -> mean photons per molecule."""
        return self.mean(np.einsum("vyl,vl->vy", pop, n))

    def readings(self, line, w):
        n, w = self.photons(line), self.weights[w]
        on = sum(w[s] * self.signal(self.fin[s], n) for s in self.fin)
        off = sum(w[s] * self.signal(self.off(s), n) for s in self.fin)
        return off, on

    def residual(self, w):
        """Mean final population per level and the J=2 total."""
        w = self.weights[w]
        per = np.array([sum(w[s] * self.mean(self.fin[s][..., k]) for s in self.fin)
                        for k in range(len(self.levels))])
        return per, per.sum()


def ratio(pt, p0, p1, tag, line="F4"):
    """CU2 on/off ratio of the probe signal for lens starts (weights over the four m_J ~ 0 levels)."""
    m = Model(pt, p0, p1, tag)
    n = m.photons(line)
    print(f"CU2 {float(m.z['E_V_cm']):g} V/cm, set point {tuple(float(x) for x in pt)} MHz, {p0} + {p1} mW; "
          f"probe {LINES[line][0]}, 15 mW, X pol, flat top 1.5 x 1.5 cm")
    for s in sorted(m.fin):
        on, off = m.signal(m.fin[s], n), m.signal(m.off(s), n)
        print(f"  start {s}: photons off {off:.3f}, on {on:.2e}, ratio on/off {on / off:.2e}")
    for w, name in (("equal", "equal over the 4 lens levels"), ("2/3", "singlet 2/3, triplets 1/9 each")):
        off, on = m.readings(line, w)
        print(f"  {name}: off {off:.3f}, on {on:.2e} photons per molecule -> on/off {on / off:.2e} "
              f"(J=2 left {100 * m.residual(w)[1]:.3f} %)")


def plot(pt, p0, p1, tag):
    """Slide schematic: how the CU2-test on/off ratio is simulated, and what the three recommended lines read."""
    global POL
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                         "mathtext.fontset": "cm", "hatch.color": "#b0b0b0", "hatch.linewidth": 0.8})
    INK, GRAY, BLUE, ORANGE, GREEN, RED = "#111111", "#555555", "#1537a6", "#d4652f", "#2e8b57", "#b8312f"
    YELLOW, RESID = "#f3d36b", "#7a7a7a"
    m = Model(pt, p0, p1, tag)
    E = float(m.z["E_V_cm"])
    init = {s: 0.0 for s in m.levels}
    init.update(m.weights["equal"])
    per, total = m.residual("equal")
    res = {s: v / total for s, v in zip(m.levels, per)}
    used = ("F4", "Q5o2F3", "Q3o2F2")
    shown = {line: dict(zip(m.levels, level_photons(line, m.levels, [0.0])[0])) for line in used}   # vx = 0
    read = {line: {w: m.readings(line, w) for w in m.weights} for line in used}
    f13 = sum(v for s, v in res.items() if s[0] == 1.5)
    f3 = sum(v for s, v in res.items() if s[1] == 3)

    groups = [(1.5, 1), (1.5, 2), (2.5, 2), (2.5, 3)]
    cols = [s for g in groups for s in sorted((s for s in m.levels if s[:2] == g), key=lambda s: s[2])]
    xs, x, prev = [], 0.0, None
    for c in cols:
        if prev is not None and c[:2] != prev:
            x += 0.55
        xs.append(x)
        x += 1.0
        prev = c[:2]
    xs = np.asarray(xs)
    x_end = xs[-1] + 1.0

    fig = plt.figure(figsize=(23.5, 15.6))
    ax = fig.add_axes([0.005, 0.005, 0.99, 0.965])
    ax.set_axis_off()
    ax.set_xlim(-7.6, 43.0)
    ax.set_ylim(16.85, -0.4)
    fig.suptitle(r"First CU2 test: how the on/off ratio of the $J=2$ detection signal is simulated",
                 fontsize=22, fontweight="bold", y=0.99)

    # simulation chain
    boxes = [
        ("Lens", r"$J=2$, the four $m_J\approx 0$" + "\nlevels at the CU2 field\n(CU1 and SPB off)", YELLOW),
        ("CU2 OBE  (on / off)", f"56 levels, {E:g} V/cm, 13 passes\n" + r"$\delta_0/2\pi$" + f" = {pt[0]:g} MHz, "
         + r"$f_{\mathrm{AOM}}$" + f" = {pt[1]:g} MHz\n{p0} + {p1} mW, 5 % loss / pass", "#dfe6f6"),
        ("Field ramp to 0", "adiabatic: each CU2 eigenstate → its\nzero-field " + r"$(F_1, F, m_F)$, $m_F$ about the"
         + "\nvertical CU2 field; turned onto " + r"$k$" + " with " + r"$d^F(\pi/2)$", "#ececec"),
        ("Probe OBE", r"one R(2) or Q(2) line, $E=0$" + "\n15 mW, " + r"$\sigma^+$, axis along $k$" + ", no EOM\n"
         "flat top 1.5 × 1.5 cm, single pass",
         "#dfe6f6"),
        ("Signal", r"$S=\langle\,\sum_l P_l\, n_l(v_x)\,\rangle$" + "\n" + r"$|v_x|\leq 1$ m/s, $|y|\leq 7.5$ mm"
         + "\n" + r"reading $S^{\mathrm{on}}/S^{\mathrm{off}}$", "#e4f1e8"),
    ]
    bw, gap, y_b, bh, x_b = 8.0, 1.85, 0.3, 2.1, -7.3
    for k, (title, body, color) in enumerate(boxes):
        x0 = x_b + k * (bw + gap)
        ax.add_patch(FancyBboxPatch((x0, y_b), bw, bh, boxstyle="round,pad=0.02,rounding_size=0.25",
                                    facecolor=color, edgecolor=GRAY, lw=1.0))
        ax.text(x0 + bw / 2, y_b + 0.38, title, ha="center", va="center", fontsize=14, fontweight="bold")
        ax.text(x0 + bw / 2, y_b + 1.3, body, ha="center", va="center", fontsize=11.5,
                linespacing=1.75 if "sum" in body else 1.35)
        if k:
            ax.add_patch(FancyArrowPatch((x0 - gap + 0.12, y_b + bh / 2), (x0 - 0.12, y_b + bh / 2),
                                         arrowstyle="-|>", mutation_scale=22, color=GRAY, lw=1.6))
    ax.text(x_b + 1.5 * (bw + gap) + bw / 2 - 2.0, y_b + bh + 0.32, r"$P_l$ = population of $J=2$ level $l$ at "
            r"the probe: the CU2 final state (on) or the start level (off)", ha="center", va="center",
            fontsize=11, color=GRAY)
    ax.text(x_b + 3.5 * (bw + gap) + bw / 2, y_b + bh + 0.32, r"$n_l$ = photons per molecule from level $l$",
            ha="center", va="center", fontsize=11, color=GRAY)

    # header: hyperfine groups and m_F
    yh = 3.5
    names = {(1.5, 1): r"$F_1=3/2,\ F=1$", (1.5, 2): r"$F_1=3/2,\ F=2$", (2.5, 2): r"$F_1=5/2,\ F=2$",
             (2.5, 3): r"$F_1=5/2,\ F=3$"}
    for g, name in names.items():
        sel = [x_ for x_, c in zip(xs, cols) if c[:2] == g]
        ax.text(np.mean(sel) + 0.5, yh, name, ha="center", va="center", fontsize=15)
        ax.plot([min(sel) + 0.05, max(sel) + 0.95], [yh + 0.27, yh + 0.27], color=GRAY, lw=1.2)
    for x_, c in zip(xs, cols):
        ax.text(x_ + 0.5, yh + 0.6, f"{c[2]:+d}".replace("-", "−") if c[2] else "0", ha="center", va="center",
                fontsize=12, color=GRAY)
    ax.text(-0.3, yh + 0.6, "$m_F$", ha="right", va="center", fontsize=13, color=GRAY)

    # (a) populations
    ax.text(-7.5, 4.55, r"(a) where the $J=2$ molecules are", fontsize=15, fontweight="bold", va="center")
    rows_a = [(5.9, "CU2 off\n(lens levels, equal)", init, YELLOW),
              (7.25, f"CU2 on: residual\n{100 * total:.2f} % of " + r"$J=2$", res, RESID)]
    for base, label, frac, color in rows_a:
        ax.text(-0.3, base - 0.45, label, ha="right", va="center", fontsize=12.5, linespacing=1.3)
        ax.plot([xs[0], x_end], [base, base], color="#cfcfcf", lw=1.0, zorder=1)
        top = max(frac.values())
        for x_, c in zip(xs, cols):
            v = frac[c]
            if v > 0:
                h = 0.95 * v / top
                ax.add_patch(Rectangle((x_ + 0.15, base - h), 0.7, h, facecolor=color, edgecolor="none", zorder=2))
                if v >= 0.045:
                    ax.text(x_ + 0.5, base - h - 0.1, f"{100 * v:.0f} %", ha="center", va="bottom", fontsize=10.5)
    ax.text(x_end + 0.4, 5.45, "% of the " + r"$J=2$" + " population; the other lens\nweighting (singlet 2/3, "
            "triplets 1/9 each)\nis the second column of (c)", ha="left", va="center", fontsize=11.5, color=GRAY,
            linespacing=1.4)
    ax.text(x_end + 0.4, 6.8, f"% of the residual: {100 * f13:.0f} % in " + r"$F_1=3/2$" + ",\nonly "
            f"{100 * f3:.0f} % in " + r"$F=3$" + " (decay from B feeds " + r"$F_1=3/2$" + ",\nwhich CU2 pumps weakly)",
            ha="left", va="center", fontsize=11.5, color=GRAY, linespacing=1.4)

    # (b) photons per molecule, (c) readings
    ax.text(-7.5, 8.1, r"(b) photons per molecule from each sublevel ($m_F$ about the vertical, $v_x=0$, 15 mW, "
            + POLNAME[POL] + ")",
            fontsize=15,
            fontweight="bold", va="center")
    ax.text(x_end + 0.4, 8.1, r"(c) predicted $S^{\mathrm{on}}/S^{\mathrm{off}}$ (true $J=2$ left "
            f"{100 * total:.2f} %)", fontsize=15, fontweight="bold", va="center")
    ax.text(x_end + 0.6, 8.6, "lens: equal     2/3, 1/9", fontsize=11, color=GRAY, va="center")
    norm, cmap = Normalize(0.0, 3.0), plt.get_cmap("Blues")
    low = sorted(m.residual(w)[1] * read["F4"][w][0] / read["F4"][w][1] for w in m.weights)
    F13, F53 = r"$F_1=3/2$", r"$F_1=5/2$"
    true13 = {w: sum(v for s, v in zip(m.levels, m.residual(w)[0]) if s[0] == 1.5) for w in m.weights}
    per_res = [read["Q3o2F2"][w][1] / true13[w] for w in m.weights]
    lensw = [read["F4"][w][0] / read["Q5o2F3"][w][0] for w in ("equal", "2/3")]
    lens53 = [shown["Q5o2F3"][s_] for s_ in m.fin]
    rows_b = [
        (8.9, "F4", "1  check\nR(2) " + r"$F_1'=7/2,\ F'=4$", ORANGE,
         r"Only $F=3$ (" + f"{100 * f3:.0f} % of the residual), saturated, no dark states:\n"
         f"checks that CU2 empties " + r"$F=3$" + f"; reads {low[0]:.0f}–{low[1]:.0f}× below the true residual."),
        (10.15, "Q5o2F3", "2  normalization\nQ(2) " + r"$F_1'=5/2,\ F'=3$", BLUE,
         f"All four lens levels ({min(lens53):.1f}–{max(lens53):.1f} photons): the CU2-off signal = the\nnumber of " + r"$J=2$" + " molecules. Its "
         r"$S^{\mathrm{off}}$" + " over that of R(2) " + r"$F'=4$" + " gives the lens\nweighting: "
         f"{lensw[0]:.2f} (equal) or {lensw[1]:.2f} (2/3, 1/9)."),
        (11.4, "Q3o2F2", "3  residual\nQ(2) " + r"$F_1'=3/2,\ F'=2$", GREEN,
         f"All 8 {F13} levels: the CU2-on signal is the residual "
         + (f"({min(per_res):.2f}" if f"{min(per_res):.2f}" == f"{max(per_res):.2f}"
            else f"({min(per_res):.2f}–{max(per_res):.2f}") + "\nphotons per residual molecule), almost background-free. "
         "The ratio is\nlarge only because its CU2-off signal is tiny."),
    ]
    for top, line, label, color, tell in rows_b:
        ax.text(-0.3, top + 0.42, label, ha="right", va="center", fontsize=12.5, linespacing=1.3)
        for x_, c in zip(xs, cols):
            v = shown[line][c]
            dark = v < 0.05
            ax.add_patch(Rectangle((x_ + 0.04, top), 0.92, 0.84, facecolor="white" if dark else cmap(norm(v)),
                                   edgecolor="#b0b0b0", lw=0.8, hatch="///" if dark else None, zorder=2))
            ax.text(x_ + 0.5, top + 0.42, "0" if dark else f"{v:.1f}", ha="center", va="center", fontsize=11,
                    color=GRAY if dark else ("white" if norm(v) > 0.6 else INK), zorder=3)
        r_eq = read[line]["equal"][1] / read[line]["equal"][0]
        r_23 = read[line]["2/3"][1] / read[line]["2/3"][0]
        ax.text(x_end + 0.6, top + 0.42, f"{100 * r_eq:.2f} %   {100 * r_23:.2f} %", ha="left", va="center",
                fontsize=15, fontweight="bold", color=color)
        ax.text(x_end + 7.6, top + 0.42, tell, ha="left", va="center", fontsize=11.5, color=color,
                linespacing=1.35)

    # two lines together: residual on Q(2) F1'=3/2 F'=2 (CU2 on), normalization on Q(2) F1'=5/2 F'=3 (CU2 off)
    pol0, two = POL, {}
    pols = ("k", "kiso", "spc") if pol0 == "k" else (pol0,)
    for p in pols:
        POL = p
        mp = Model(pt, p0, p1, tag)
        for w in mp.weights:
            two[p, w] = mp.readings("Q3o2F2", w)[1] / mp.readings("Q5o2F3", w)[0]
    POL = pol0
    fac = [true13[w] / two[p, w] for p in pols for w in m.weights]
    yb = 12.75
    ax.add_patch(Rectangle((-7.4, yb - 0.1), 50.2, 2.35, facecolor="#f2f5fb", edgecolor=BLUE, lw=1.2, zorder=0))
    ax.text(-7.2, yb + 0.15, "Lines 2 and 3 together", fontsize=14, fontweight="bold", color=BLUE, va="top")
    ax.text(-7.1, yb + 1.05, r"$\frac{S^{\mathrm{on}}\,[\mathrm{Q}(2)\ F_1'=3/2,\ F'=2]}"
            r"{S^{\mathrm{off}}\,[\mathrm{Q}(2)\ F_1'=5/2,\ F'=3]}$", fontsize=21, color=INK, va="center")
    ax.text(6.3, yb + 0.75, f"= {100 * two[pol0, 'equal']:.2f} % / {100 * two[pol0, '2/3']:.2f} % "
            + ("(axis along k),  " + f"{100 * two['kiso', 'equal']:.2f} % / {100 * two['kiso', '2/3']:.2f} % "
               "(alignment scrambled),  " + f"{100 * two['spc', 'equal']:.2f} % / {100 * two['spc', '2/3']:.2f} % "
               + r"(exact, earth field $\perp k$)" if pol0 == "k" else ""),
            fontsize=13, color=INK, va="center", fontweight="bold")
    ax.text(6.3, yb + 1.45, f"true {F13} residual {100 * true13['equal']:.2f} % / {100 * true13['2/3']:.2f} %  "
            f"→  one model factor {min(fac):.1f}–{max(fac):.1f} for both weightings and all three treatments "
            f"(lens equal / 2/3, 1/9)", fontsize=12.5, color=INK, va="center")
    ax.text(-7.2, yb + 1.95, r"Axis along $k$: the light is pure $\sigma^+$; on Q(2) ($\Delta F=0$) the stretched "
            r"$m_F'=+F$ along $k$ is dark, R(2) $F'=4$ ($F\to F+1$) has none. The 0.5 G earth field ($\perp k$, "
            r"plates off) turns a molecule by $\leq 0.4$ rad in the 81 µs probe.",
            fontsize=11.5, color=INK, va="center")

    ax.text(-7.5, 15.95,
            f"CU2 at {E:g} V/cm, set point ({pt[0]:g}, {pt[1]:g}) MHz, {p0} + {p1} mW. Lens start = the four "
            r"$J=2$ levels with $P(m_J=0)>0.5$ at the CU2 field ($F_1=5/2$: $F=2,\ m_F=0$ and $F=3,\ m_F=0,\pm1$)."
            + "\n" + r"Means over $\sigma_{v_x}=0.65$ m/s cut at $|v_x|\leq 1$ m/s and $\sigma_y=10.35$ mm cut at the "
            r"probe height $|y|\leq 7.5$ mm (CU2 grid 6 $v_x$ × 4 $y_0$); (b) is $v_x=0$, (c) includes the Doppler "
            r"shift $k v_x$." + "\n" + r"Only the probed line is in the probe model; no background. centrex decay "
            r"operators are pairwise (no coherence transfer), adequate only with the axis along the light's symmetry "
            r"axis; the exact check (vertical axis, coherent decay) agrees to $\leq 3$ %.",
            fontsize=11, color=GRAY, va="center", linespacing=1.6)
    fig_dir = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
    fig_dir.mkdir(parents=True, exist_ok=True)
    name = f"CU2_test_onoff_obe_e{E:g}_d{pt[0]:g}_f{pt[1]:g}".replace(".", "p")
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"{name}.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", fig_dir / name)


if __name__ == "__main__":
    POL = next((a[3:] for a in sys.argv[1:] if a.startswith("pol") and a[3:] in POL_VEC), "k")
    args = [a for a in sys.argv[1:] if a not in LINES and not a.startswith("pol")]
    line = next((a for a in sys.argv[1:] if a in LINES), "F4")
    if args[0] == "yields":
        yields(line)
    elif args[0] == "yields_coherent":
        yields_coherent(line)
    elif args[0] in ("ratio", "plot"):
        pt, p0, p1 = tuple(float(x) for x in args[1].split(",")), int(args[2]), int(args[3])
        tag = f"_{args[4]}" if len(args) > 4 else ""
        if args[0] == "ratio":
            ratio(pt, p0, p1, tag, line)
        else:
            plot(pt, p0, p1, tag)
