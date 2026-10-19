"""Level diagrams of the two J=1 detection lines, with the photon yield of every ground sublevel.

    python scripts/detection_level_diagrams.py compute   # probe OBE: photons per molecule per J=1 sublevel (~4 min)
    python scripts/detection_level_diagrams.py plot      # diagrams (slides style) from the cached yields
    python scripts/detection_level_diagrams.py info      # slide schematic: what each line tells about CU1

Lines: R(1) F1'=5/2 F'=3 (the AP-paper detection line, sees F=2 only) and Q(1) F1'=1/2 F'=1 (sees all of J=1).
Probe model: the present detection laser, 2.5 cm along the molecular beam x 1.7 cm high, taken as a Gaussian with
these 1/e^2 diameters; single frequency on resonance (detuning = Doppler shift k*vx), fixed linear polarization,
no phase modulation; zero field; molecule at 184 m/s. Output is the number of photons scattered by a molecule
that starts in each J=1 sublevel, on a grid of peak intensity, vx and height y.

Polarization "par" is pi light along the quantization axis (the axis the molecules keep from the last electric
field), "perp" is linear light perpendicular to it. The solver has one collapse operator per (e, g) pair, i.e.
no decay-induced ground-state coherence; this is exact for "par" with molecules in m_F eigenstates, approximate
for "perp".

The diagrams show the whole excited rotational level (all F1', F', both parities) and the 12 J=1 sublevels,
in energy order but not to scale; the label columns sit at the height of the levels they name.
"""
import itertools
import json
import sys
import textwrap
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "tmp" / "detection" / "probe_yields.json"
FIG_DIR = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
HDF_J1 = ROOT / "results" / "SPB_J1_cleanup_scans.h5"

VZ = 184.0                                   # m/s
SIGMA_Z, SIGMA_Y = 6.25e-3, 4.25e-3          # probe sigma along the beam and in height (1/e^2 diameters / 4)
T_PROBE = 6 * SIGMA_Z / VZ                   # +-3 sigma
K_VX = 2 * np.pi * 3.68e6                    # k*vx per m/s at 271.7 nm, rad/s
INTENSITIES = [20.0, 100.0, 500.0]           # peak, W/m^2 (2, 10, 50 mW/cm^2)
VXS = [0.0, 0.5, 1.0]                        # m/s
YS = [0.0, 8.5e-3]                           # m, centre and edge of the detection acceptance
SHOWN = (500.0, 0.0, 0.0)                    # grid point drawn on the diagrams (saturated, centre, vx = 0)
LINES = {
    "R1F3": dict(name="R(1) F1'=5/2, F'=3", main_g=dict(J=1, F1=1.5, F=2, mF=0), main_e=dict(J=2, F1=2.5, F=3, P=1)),
    "Q1F1": dict(name="Q(1) F1'=1/2, F'=1", main_g=dict(J=1, F1=0.5, F=0, mF=0), main_e=dict(J=1, F1=0.5, F=1, P=1)),
}


def build(key, pol):
    from centrex_tlf import couplings, lindblad, states, transitions
    from centrex_tlf.couplings.polarization import Polarization
    warnings.simplefilter("ignore")
    X, B = states.ElectronicState.X, states.ElectronicState.B
    trans = {"R1F3": transitions.R1_F1_5o2_F3, "Q1F1": transitions.Q1_F1_1o2_F1}[key]
    line = LINES[key]
    z = Polarization(np.array([0, 0, 1.0], dtype=complex), name="Z")
    x = Polarization(np.array([1.0, 0, 0], dtype=complex), name="X")
    # main pair: Delta m_F = 0 for pi light, +1 for perpendicular light, so main_coupling is never zero
    me = dict(line["main_e"], mF=0 if pol == "par" else 1)
    g = 1 * next(iter(states.generate_coupled_states_X(states.QuantumSelector(electronic=X, P=-1, **line["main_g"]))))
    e = 1 * next(iter(states.generate_coupled_states_B(states.QuantumSelector(electronic=B, **me))))
    sel = couplings.generate_transition_selectors(transitions=[trans], polarizations=[[z, x] if pol == "par" else [x, z]],
                                                  ground_mains=[g], excited_mains=[e])
    obe = lindblad.generate_OBE_system_transitions([trans], sel, qn_compact=True, E=np.array([0.0, 0.0, 0.0]),
                                                   B=np.array([0.0, 0.0, 1e-3]), retain_opposite_parity_levels=False)
    p = lindblad.LindbladParameters()
    t = p.time()
    I0, det, y = p.real("I0", 100.0), p.real("det", 0.0), p.real("y", 0.0)
    mc = float(abs(obe.couplings[0].main_coupling))
    p.bind(sel[0].Ω, lindblad.gaussian_beam_rabi(VZ * t, y, I0, 3 * SIGMA_Z, 0.0, SIGMA_Z, SIGMA_Y, mc), finalize=False)
    p.bind(sel[0].δ, det, finalize=False)
    a, b = sel[0].polarization_symbols          # first polarization of the selector is the one used
    p.bind(a, 1.0, finalize=False); p.bind(b, 0.0, finalize=False)
    p._finalize()
    prepared = lindblad.prepare_lindblad_problem(obe, p, backend="rust", hamiltonian_representation="decomposed")
    return prepared, obe


def compute():
    from centrex_tlf import hamiltonian, lindblad, states
    gamma = float(getattr(hamiltonian, "Γ"))
    X, B = states.ElectronicState.X, states.ElectronicState.B
    combos = list(itertools.product(INTENSITIES, VXS, YS))
    out = dict(grid=dict(intensity_W_m2=INTENSITIES, vx_m_s=VXS, y_m=YS), probe=dict(
        sigma_z_m=SIGMA_Z, sigma_y_m=SIGMA_Y, vz_m_s=VZ, t_end_s=T_PROBE, frequency="single, on resonance",
        polarization="fixed linear", field="E = 0, B = 1e-3 G along z"), runs={})
    for key, pol in itertools.product(LINES, ("par", "perp")):
        t0 = time.time()
        prepared, obe = build(key, pol)
        QN, dim = obe.QN, len(obe.QN)
        gidx = [i for i, s in enumerate(QN) if s.largest.electronic_state == X and s.largest.J == 1]
        exc = [i for i, s in enumerate(QN) if s.largest.electronic_state == B]
        rhos, rows = [], []
        for I0, vx, y in combos:
            for i in gidx:
                rho = np.zeros((dim, dim), dtype=complex)
                rho[i, i] = 1.0
                rhos.append(prepared.layout.pack(rho)); rows.append([I0, K_VX * vx, y])
        # fixed step: the adaptive solvers stop with a stiffness error in the dark tail of the probe
        res = lindblad.solve_lindblad_batch(prepared, np.asarray(rhos), (0.0, T_PROBE), parameter_slots=["I0", "det", "y"],
                                            parameter_batch=np.asarray(rows, dtype=np.complex128), solver="fixed_dopri5",
                                            execution_mode="expanded_sparse", dt=2e-9, maxiters=5_000_000,
                                            output="photon_integral", integral_weights=[(i, gamma) for i in exc],
                                            output_when="final", parallel=True)
        ph = np.asarray(res.values).real.reshape(len(combos), len(gidx))
        out["runs"][f"{key}_{pol}"] = dict(
            labels=[[float(QN[i].largest.F1), int(QN[i].largest.F), int(QN[i].largest.mF)] for i in gidx],
            combos=[list(c) for c in combos], photons=ph.tolist())
        print(f"{key} {pol}: {time.time() - t0:.0f} s", flush=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(out, indent=1))
    print("wrote", CACHE)


def yields_at(data, run, point=SHOWN):
    r = data["runs"][run]
    k = [tuple(c) for c in r["combos"]].index(tuple(point))
    return {tuple(lab): v for lab, v in zip(r["labels"], r["photons"][k])}


def cu1_residual_split():
    """Where the CU1 residual sits at +19.0 MHz (4.8 cm, on axis) for the AP1 triplet start."""
    import h5py
    with h5py.File(HDF_J1, "r") as h5:
        g = h5["runs"]["trajectories_cu1_vx0_19p0_20p2MHz_4p8cm"]
        r = g["results"][0, :, -1, 1:13].real
        labels = json.loads(g.attrs["fixed_parameters_json"])["level_labels"][1:13]
    trip = [k for k, lab in enumerate(labels) if "F1=1.5, F=2" in lab and abs(int(lab[-2:])) <= 1]
    fin = r[trip].mean(axis=0)
    f0 = sum(fin[k] for k, lab in enumerate(labels) if "F=0" in lab)
    f22 = sum(fin[k] for k, lab in enumerate(labels) if "F1=1.5, F=2" in lab and abs(int(lab[-2:])) == 2)
    return fin.sum(), f0 / fin.sum(), f22 / fin.sum()


def excited_manifold(key):
    """Energies (MHz) at E = 0 of every hyperfine level of the excited rotational level, both parities,
    and of the four J=1 ground hyperfine levels. Same Hamiltonian and origin for all calls."""
    from centrex_tlf.utils import plotting as lp
    Je, branch = {"R1F3": (2, "R"), "Q1F1": (1, "Q")}[key]
    exc, ground = {}, None
    for F1 in (Je - 0.5, Je + 0.5):
        for F in (int(F1 - 0.5), int(F1 + 0.5)):
            s = lp.calculate_transition_level_structure(E=0.0, B=0.0, J_ground=1, branch=branch, F1_excited=F1,
                                                        F_excited=F)
            ground = s.ground
            for lv in s.excited:                  # all m_F' of one (F1', F', P') are degenerate at E = 0
                exc[(lv.F1, lv.F, lv.P)] = lv.energy_MHz
    gl = {(lv.F1, lv.F): lv.energy_MHz for lv in ground}
    return Je, exc, gl


def plot():
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, Rectangle

    data = json.loads(CACHE.read_text())
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    YELLOW, GREY, INK, BLUE = "#f3d36b", "#b9b9b9", "#111111", "#1537a6"
    PCOL = {-1: "#1f77b4", 1: "#ff7f0e"}                                 # excited parity P' = -, +
    GCOL = {(1.5, 2): "#2ca02c", (1.5, 1): "#d62728", (0.5, 1): "#9467bd", (0.5, 0): "#8c564b"}
    width, pad, lw = 0.34, 0.03, 4.2
    total, frac_f0, frac_f22 = cu1_residual_split()
    cu1_off = {(1.5, 2, -1), (1.5, 2, 0), (1.5, 2, 1), (1.5, 1, 0)}       # filled by AP1 (triplet, singlet)
    cu1_res = {(0.5, 0, 0), (1.5, 2, -2), (1.5, 2, 2)}                     # where CU1 leaves its residual
    target = {"R1F3": (2.5, 3, 1), "Q1F1": (0.5, 1, 1)}
    I0 = SHOWN[0]

    def fmt(v):
        return "0" if v < 0.005 else f"{v:.1f}"

    def half(v):
        return f"{int(round(2 * v))}/2" if abs(v - round(v)) > 0.25 else f"{v:.0f}"

    def figure(key, fname):
        par, perp = yields_at(data, f"{key}_par"), yields_at(data, f"{key}_perp")
        Je, exc, gl = excited_manifold(key)
        tF1, tF, tP = target[key]
        e0 = exc[target[key]]
        mfmax = int(max(F for _, F, _ in exc))
        fig = plt.figure(figsize=(16.5, 12.0))
        ax = fig.add_axes([0.02, 0.24, 0.96, 0.69])
        ax.set_axis_off()
        xP = mfmax + 0.75                                                  # label columns right of the levels
        xF, xF1, xO = xP + 0.6, xP + 1.25, xP + 1.9
        ax.set_xlim(-mfmax - 3.3, xO + 3.7)
        ax.set_ylim(0.0, 1.0)

        # ---- excited: hyperfine groups in energy order (not to scale), doublet splitting exaggerated -----------
        groups = sorted({(F1, F) for F1, F, _ in exc}, key=lambda g: exc[(g[0], g[1], 1)] + exc[(g[0], g[1], -1)])
        y_exc = {}
        for (F1, F), yc in zip(groups, np.linspace(0.50, 0.87, len(groups))):
            lo, hi = sorted((-1, 1), key=lambda p: exc[(F1, F, p)])
            y_exc[(F1, F, lo)], y_exc[(F1, F, hi)] = yc - 0.018, yc + 0.018
            ax.text(xF, yc, f"{F}", ha="center", va="center", fontsize=13)      # F', F1' at the doublet centre
            ax.text(xF1, yc, half(F1), ha="center", va="center", fontsize=13)
        for (F1, F, P), yy in y_exc.items():
            for m in range(-F, F + 1):
                ax.plot([m - width / 2, m + width / 2], [yy, yy], color=PCOL[P], lw=lw, solid_capstyle="butt", zorder=2)
            ax.text(xP, yy, "+" if P > 0 else "−", ha="center", va="center", fontsize=13, color=PCOL[P],
                    fontweight="bold")                                     # P' at the height of its own line
            is_t = (F1, F, P) == target[key]
            off = f"{exc[(F1, F, P)] - e0:+,.1f}".replace(",", " ").replace("-", "−")
            ax.text(xO, yy, "0   ← probe on resonance" if is_t else off, ha="left", va="center",
                    fontsize=11.5 if is_t else 11, color=BLUE if is_t else INK, fontweight="bold" if is_t else "normal")
        ytop = max(y_exc.values()) + 0.055
        for x_, t_ in ((xP, "P'"), (xF, "F'"), (xF1, "F₁'")):
            ax.text(x_, ytop, t_, ha="center", va="center", fontsize=13, style="italic")
        ax.text(xO, ytop, "offset from the probed\nlevel (MHz)", ha="left", va="center", fontsize=11)
        ax.text(-mfmax - 1.0, 0.975, r"$m_F'$", ha="center", va="center", fontsize=14)
        for m in range(-mfmax, mfmax + 1):
            ax.text(m, 0.975, f"{m:+d}".replace("-", "−") if m else "0", ha="center", va="center", fontsize=13)
        ax.text(-mfmax - 3.2, 0.69, rf"$B\,^3\Pi_1\ (v'=0)$" + "\n" + rf"$J'={Je}$", ha="left", va="center",
                fontsize=16, linespacing=1.6)

        # probe: dashed line through the probed level, outlines of the driven sublevels
        yt = y_exc[target[key]]
        ax.plot([-mfmax - 0.55, mfmax + 0.45], [yt, yt], ls=(0, (6, 3)), lw=1.8, color=BLUE, zorder=1)
        for m in range(-tF, tF + 1):
            ax.add_patch(Rectangle((m - width / 2 - pad, yt - 0.011), width + 2 * pad, 0.022, facecolor="none",
                                   edgecolor=BLUE, lw=2.2, zorder=3))

        # ---- ground J = 1: hyperfine levels in energy order, the two F1 pairs apart (not to scale) ------------
        gorder = sorted(gl, key=lambda g: gl[g])                          # (1/2,0), (1/2,1), (3/2,1), (3/2,2)
        y_g = dict(zip(gorder, [0.085, 0.135, 0.255, 0.305]))
        same = all(abs(par[k] - perp[k]) < 0.05 for k in par)
        for (F1, F), yy in y_g.items():
            for m in range(-F, F + 1):
                k = (F1, F, m)
                band = YELLOW if k in cu1_off else GREY if k in cu1_res else None
                if band:                                                   # larger than the outline, so it shows
                    ax.add_patch(Rectangle((m - width / 2 - pad - 0.05, yy - 0.016), width + 2 * pad + 0.10, 0.032,
                                           facecolor=band, edgecolor="none", zorder=0))
                ax.plot([m - width / 2, m + width / 2], [yy, yy], color=GCOL[(F1, F)], lw=lw, solid_capstyle="butt",
                        zorder=2)
                lo, hi = min(par[k], perp[k]), max(par[k], perp[k])
                if hi > 0.05:                                              # detected in at least one polarization
                    ax.add_patch(Rectangle((m - width / 2 - pad, yy - 0.010), width + 2 * pad, 0.020,
                                           facecolor="none", edgecolor=BLUE, lw=2.0,
                                           ls="solid" if lo > 0.05 else (0, (3, 2)), zorder=3))
                lab = fmt(par[k]) if same else f"{fmt(par[k])}|{fmt(perp[k])}"
                ax.text(m + width / 2 + pad + 0.07, yy, lab, ha="left", va="center", fontsize=10.5,
                        color=INK if hi > 0.05 else "0.55", fontweight="bold" if hi > 0.05 else "normal")
            ax.text(xF, yy, f"{F}", ha="center", va="center", fontsize=13, color=GCOL[(F1, F)])
            ax.text(xO, yy, f"{1e3 * (gl[(F1, F)] - gl[gorder[-1]]):+.1f}".replace("-", "−").replace("+0.0", "0"),
                    ha="left", va="center", fontsize=11)
        for F1 in (0.5, 1.5):                                              # F1 bracket at the centre of its pair
            ys_ = [y_g[g] for g in y_g if g[0] == F1]
            ax.plot([xF1 - 0.3, xF1 - 0.3], [min(ys_) - 0.012, max(ys_) + 0.012], color=INK, lw=1.4)
            ax.text(xF1, np.mean(ys_), half(F1), ha="center", va="center", fontsize=13)
        yg_top = max(y_g.values()) + 0.06
        for x_, t_ in ((xF, "F"), (xF1, "F₁")):
            ax.text(x_, yg_top, t_, ha="center", va="center", fontsize=13, style="italic")
        ax.text(xO, yg_top, "offset from F₁=3/2, F=2\n(kHz)", ha="left", va="center", fontsize=11)
        ax.text(-mfmax - 1.0, 0.025, r"$m_F$", ha="center", va="center", fontsize=14)
        for m in range(-2, 3):
            ax.text(m, 0.025, f"{m:+d}".replace("-", "−") if m else "0", ha="center", va="center", fontsize=13)
        ax.text(-mfmax - 3.2, 0.20, r"$X\,^1\Sigma^+\ (v=0)$" + "\n" + r"$J=1$", ha="left", va="center", fontsize=16,
                linespacing=1.6)
        bx0, bx1, by0, by1 = -2 - width / 2 - 0.2, 2 + width / 2 + 0.85, 0.055, 0.335
        ax.add_patch(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fill=False, ls=(0, (5, 3)), lw=1.5,
                               edgecolor="0.3", zorder=4))
        # arrow between the m_F columns, so it does not run over the bars it passes
        ax.add_patch(FancyArrowPatch((0.5, by1), (0.5, yt), arrowstyle="-|>", mutation_scale=18, lw=2.0,
                                     color=BLUE, shrinkA=0, shrinkB=0, zorder=2.5))
        fig.suptitle(f"Detection  |  {LINES[key]['name'].replace('F1', 'F₁')}  |  E = 0  |  probe 2.5 cm × 1.7 cm, "
                     "single frequency", fontsize=19, y=0.985)

        # ---- key -------------------------------------------------------------------------------------------
        p = f"{I0 / 10:.0f} mW/cm² peak, probe centre, vx = 0"
        if same:
            num = (f"Number at each J=1 sublevel: photons per molecule, saturated ({p}), any polarization. "
                   "Blue box: detected; grey: not detected.")
            how = "F'=3 decays inside J=1 only to F=2: no dark state, 2.1 photons per F=2 molecule."
        else:
            num = (f"Number at each J=1 sublevel: photons per molecule, saturated ({p}), polarization ∥ | ⊥ to the "
                   "molecules' quantization axis (⊥ approximate). Solid box: detected for both; dashed: only for ⊥.")
            how = "3 excited vs 12 ground sublevels: fixed polarization leaves dark states. F=0 is never dark."
        key_text = [
            f"Blue dashed line: probe on resonance with |F₁'={half(tF1)}, F'={tF}, P'=+⟩, single frequency, fixed "
            "linear polarization; blue boxes on B: driven.",
            num,
            how,
            f"Yellow: populated with CU1 off (AP1).  Grey: where the CU1 residual sits ({100 * frac_f0:.0f} % F=0, "
            f"{100 * frac_f22:.0f} % F=2, m_F=±2; +19.0 MHz).",
            "Levels in energy order, not to scale; right columns: offsets in MHz (B) and kHz (X).",
        ]
        wrapped = [w for par_ in key_text for w in textwrap.wrap(par_, 140)]
        for k, line in enumerate(wrapped):
            fig.text(0.02, 0.215 - 0.024 * k, line, ha="left", va="top", fontsize=15, color=INK)
        for ext in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"{fname}.{ext}", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("wrote", fname)

    figure("R1F3", "detection_level_diagram_R1_F3_0Vcm")
    figure("Q1F1", "detection_level_diagram_Q1_F1_0Vcm")


def cu1_residual_by_level():
    """J=1 population per sublevel after CU1 (+19.0 MHz, 4.8 cm, on axis) for the AP1 triplet start
    (1/3 in each of |F1=3/2, F=2, m_F=-1,0,+1>)."""
    import h5py
    with h5py.File(HDF_J1, "r") as h5:
        g = h5["runs"]["trajectories_cu1_vx0_19p0_20p2MHz_4p8cm"]
        r = g["results"][0, :, -1, 1:13].real
        labels = json.loads(g.attrs["fixed_parameters_json"])["level_labels"][1:13]
    keys = []
    for lab in labels:                                   # "X J=1, F1=0.5, F=0, mF=+0"
        f1, f, mf = (part.split("=")[1] for part in lab.split(", ")[1:])
        keys.append((float(f1), int(f), int(mf)))
    trip = [k for k, key in enumerate(keys) if key[:2] == (1.5, 2) and abs(key[2]) <= 1]
    fin = r[trip].mean(axis=0)
    return dict(zip(keys, fin))


def plot_info():
    """Slide schematic: what each detection configuration tells about CU1."""
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle

    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                         "mathtext.fontset": "cm", "hatch.color": "#b0b0b0", "hatch.linewidth": 0.8})
    INK, GRAY, BLUE, ORANGE, GREEN, RED = "#111111", "#555555", "#1537a6", "#d4652f", "#2e8b57", "#b8312f"
    YELLOW, RESID = "#f3d36b", "#7a7a7a"
    data = json.loads(CACHE.read_text())
    yR, yP, yX = (yields_at(data, r) for r in ("R1F3_par", "Q1F1_par", "Q1F1_perp"))
    fin = cu1_residual_by_level()
    total = sum(fin.values())
    init = {k: (1 / 3 if k[:2] == (1.5, 2) and abs(k[2]) <= 1 else 0.0) for k in fin}
    groups = {}
    for k in fin:
        groups.setdefault(k[:2], []).append(k)
    yI = {k: np.mean([yP[j] for j in groups[k[:2]]]) for k in fin}       # alignment scrambled (exact pi average)

    def ratio(y):
        return sum(fin[k] * y[k] for k in fin) / sum(init[k] * y[k] for k in fin)

    nF2 = sum(fin[k] * yR[k] for k in fin) / 2.10
    yF0 = yP[(0.5, 0, 0)]
    bound = [nF2 + sum(fin[k] * y[k] for k in fin) / yF0 for y in (yP, yX, yI)]
    rR, rP, rX = ratio(yR), ratio(yP), ratio(yX)
    f0, f22 = fin[(0.5, 0, 0)] / total, (fin[(1.5, 2, -2)] + fin[(1.5, 2, 2)]) / total

    cols = ([(1.5, 2, m) for m in range(-2, 3)] + [(1.5, 1, m) for m in (-1, 0, 1)]
            + [(0.5, 1, m) for m in (-1, 0, 1)] + [(0.5, 0, 0)])
    xs, x, prev = [], 0.0, None
    for c in cols:
        if prev is not None and c[:2] != prev:
            x += 0.55
        xs.append(x)
        x += 1.0
        prev = c[:2]
    xs = np.asarray(xs)
    x_end = xs[-1] + 1.0

    fig = plt.figure(figsize=(18.5, 10.6))
    ax = fig.add_axes([0.01, 0.01, 0.98, 0.92])
    ax.set_axis_off()
    ax.set_xlim(-7.4, 31.5)
    ax.set_ylim(11.9, -0.4)
    fig.suptitle(r"What each $J=1$ detection line tells about CU1", fontsize=21, fontweight="bold", y=0.985)

    # header: hyperfine groups and m_F
    names = {(1.5, 2): r"$F_1=3/2,\ F=2$", (1.5, 1): r"$F_1=3/2,\ F=1$", (0.5, 1): r"$F_1=1/2,\ F=1$",
             (0.5, 0): r"$F=0$"}
    for g, name in names.items():
        sel = [x_ for x_, c in zip(xs, cols) if c[:2] == g]
        ax.text(np.mean(sel) + 0.5, 0.15, name, ha="center", va="center", fontsize=15)
        ax.plot([min(sel) + 0.05, max(sel) + 0.95], [0.42, 0.42], color=GRAY, lw=1.2)
    for x_, c in zip(xs, cols):
        ax.text(x_ + 0.5, 0.75, f"{c[2]:+d}".replace("-", "−") if c[2] else "0", ha="center", va="center",
                fontsize=12, color=GRAY)
    ax.text(-0.3, 0.75, "$m_F$", ha="right", va="center", fontsize=13, color=GRAY)

    # (a) populations
    ax.text(-7.3, 1.25, r"(a) where the $J=1$ molecules are", fontsize=15, fontweight="bold", va="center")
    rows_a = [(2.55, "CU1 off (after AP1,\ntriplet)", init, YELLOW, 1.0),
              (3.85, f"CU1 on: residual\n{100 * total:.2f} % of " + r"$J=1$", {k: v / total for k, v in fin.items()},
               RESID, 1.0)]
    for base, label, frac, color, _ in rows_a:
        ax.text(-0.3, base - 0.45, label, ha="right", va="center", fontsize=12.5, linespacing=1.3)
        ax.plot([xs[0], x_end], [base, base], color="#cfcfcf", lw=1.0, zorder=1)
        for x_, c in zip(xs, cols):
            v = frac[c]
            if v > 0:
                h = 0.95 * v / max(frac.values())
                ax.add_patch(Rectangle((x_ + 0.15, base - h), 0.7, h, facecolor=color, edgecolor="none", zorder=2))
                if v >= 0.05:
                    ax.text(x_ + 0.5, base - h - 0.12, f"{100 * v:.0f} %", ha="center", va="bottom", fontsize=11)
    ax.text(x_end + 0.3, 2.1, "% of the " + r"$J=1$" + " population\n(singlet → " + r"$F_1=3/2,\ F=1,\ m_F=0$"
            + " if\n" + r"$J=0,\ F=0$" + " was filled)", ha="left", va="center", fontsize=11.5, color=GRAY,
            linespacing=1.4)
    ax.text(x_end + 0.3, 3.4, f"% of the residual; the other\n{100 * (1 - f0 - f22):.0f} % are spread thinly",
            ha="left", va="center", fontsize=11, color=GRAY, linespacing=1.3)

    # (b) photons per molecule; (c) reading and meaning
    ax.text(-7.3, 4.65, "(b) photons per molecule from each sublevel", fontsize=15, fontweight="bold", va="center")
    ax.text(x_end + 0.3, 4.65, r"(c) predicted $S^{\mathrm{on}}/S^{\mathrm{off}}$, and what it tells about CU1",
            fontsize=15, fontweight="bold", va="center")
    norm = Normalize(0.0, 4.6)
    cmap = plt.get_cmap("Blues")
    F0, F2, m2 = r"$F=0$", r"$F=2$", r"$F=2,\ m_F=\pm 2$"
    rows_b = [
        (5.05, "R(1) " + r"$F_1'=5/2,\ F'=3$" + "\nany polarization", yR, rR, GREEN,
         f"{F2} only, no dark states → the {F2} part of the residual,\n"
         f"exactly. Blind to {F0} ({100 * f0:.0f} % of the residual) and to the singlet."),
        (6.3, "Q(1) " + r"$F_1'=1/2,\ F'=1$" + "\npolarization " + r"$\parallel$" + " axis", yP, rP, ORANGE,
         f"{F0} is never dark (4.4), but {m2} is: misses {100 * f22:.0f} %\n"
         f"of the residual. Gives an upper bound on the {F0} part."),
        (7.55, "Q(1) " + r"$F_1'=1/2,\ F'=1$" + "\npolarization " + r"$\perp$" + " axis", yX, rX, RED,
         "Sees all 12, but with 0.9–3.7 photons: the CU1-off levels\n"
         f"weigh little, the residual a lot → {rX / total:.1f}× too high (approx.)."),
    ]
    for top, label, y, r_, color, tell in rows_b:
        ax.text(-0.3, top + 0.42, label, ha="right", va="center", fontsize=12.5, linespacing=1.3)
        for x_, c in zip(xs, cols):
            v = y[c]
            dark = v < 0.05
            ax.add_patch(Rectangle((x_ + 0.04, top), 0.92, 0.84, facecolor="white" if dark else cmap(norm(v)),
                                   edgecolor="#b0b0b0", lw=0.8, hatch="///" if dark else None, zorder=2))
            ax.text(x_ + 0.5, top + 0.42, "0" if dark else f"{v:.1f}", ha="center", va="center", fontsize=12,
                    color=GRAY if dark else ("white" if norm(v) > 0.55 else INK), zorder=3)
        ax.text(x_end + 0.6, top + 0.42, f"{100 * r_:.2f} %", ha="left", va="center", fontsize=16,
                fontweight="bold", color=color)
        ax.text(x_end + 4.6, top + 0.42, tell, ha="left", va="center", fontsize=12, color=color, linespacing=1.35)

    # combined
    yb = 8.85
    ax.add_patch(Rectangle((-7.2, yb - 0.1), 38.4, 2.25, facecolor="#f2f5fb", edgecolor=BLUE, lw=1.2, zorder=0))
    ax.text(-7.0, yb + 0.12, "Both lines together", fontsize=14, fontweight="bold", color=BLUE, va="top")
    formula = (r"$\frac{N_{J=1}^{\mathrm{on}}}{N_{J=1}^{\mathrm{off}}}\ \lesssim\ "
               r"\frac{S_R^{\mathrm{on}}/2.1\ +\ S_Q^{\mathrm{on}}/4.4}{S_R^{\mathrm{off}}/2.1}$")
    ax.text(-6.9, yb + 0.95, formula, fontsize=22, color=INK, va="center")
    ax.text(0.7, yb + 0.95, f"= {100 * min(bound):.2f}–{100 * max(bound):.2f} % for any alignment "
            f"(true {100 * total:.2f} %)", fontsize=14, color=INK, va="center", fontweight="bold")
    ax.text(13.6, yb + 0.95,
            r"$S_R$, $S_Q$: signals on R(1) $F'=3$ and Q(1) $F'=1$;  2.1, 4.4: photons per" + "\n"
            + r"$F=2$ molecule (R) and per $F=0$ molecule (Q).  $F=2$ exact, $F=0$ bounded," + "\n"
            + rf"the $F=1$ levels ({100 * (1 - f0 - f22):.0f} %) from the model.",
            fontsize=12.5, color=INK, va="center", linespacing=1.5)
    ax.text(-7.0, yb + 1.85,
            r"Q(1) at $\parallel$ and $\perp$:  $S_Q^{\mathrm{on}}/S_Q^{\mathrm{off}}$"
            + f" = {100 * rP:.2f} % vs {100 * rX:.2f} % if the molecules keep the CU1 field axis, the same if their "
            "alignment is scrambled → tests the alignment and the CU1 model.", fontsize=12.5, color=INK, va="center")

    ax.text(-7.3, 11.6,
            f"Readings predicted for CU1 at +19.0 MHz, 4.8 cm beam, on axis, AP1 triplet start (true " + r"$J=1$"
            + f" residual {100 * total:.2f} %); probe saturated ({SHOWN[0] / 10:.0f} mW/cm², centre, "
            + r"$v_x=0$" + "), zero field.\n" + r"$\parallel$ / $\perp$" + ": probe polarization relative to the "
            "molecules' quantization axis; " + r"$\perp$" + " approximate (no decay-induced ground coherence in the "
            "solver).", fontsize=11, color=GRAY, va="center", linespacing=1.5)
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"detection_cu1_information.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote detection_cu1_information; readings R, par, perp:", rR, rP, rX, "bound", bound)


if __name__ == "__main__":
    {"compute": compute, "plot": plot, "info": plot_info}[sys.argv[1]]()
