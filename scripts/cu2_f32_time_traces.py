"""m_F-resolved population of X J=2 F1=3/2 F=2 along the CU2 interaction region (same model as the transfer stage
of scripts/cu2_power_check.py: 217 V/cm, (delta0, f_AOM) = (13, 209) MHz, 20 + 30 mW, 5 % loss per pass), at the measured forward velocity VZ = 156 m/s (184 m/s with `vz184`).

    python scripts/cu2_f32_time_traces.py compute [vz184]   # 8 trajectories with saved populations (~20 min)
    python scripts/cu2_f32_time_traces.py plot [vz184]      # -> docs/state_cleanup/figures/slides/CU2_F1eq3o2_F2_vs_z_vz156.{png,pdf} (or _vz184)

Starts: F1=3/2 F=2 mF = 0, +1, +2 (mF < 0 mirror these) and the singlet F1=5/2 F=2 mF=0 (it feeds F1=3/2 F=2);
vx = 0.2 and 0.6 m/s, y0 = 0. z = vz * t from the start of the 5 cm window; pass centres at 3.5 mm + 3.5 mm * i.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "tmp" / "cu2_power"
VZ = 184.0 if "vz184" in sys.argv else 156.0           # m/s; 156 = measured (TOF), 184 = old model default
TAG = f"_vz{VZ:g}"                                     # figure names always carry vz
from scripts import vz_paths  # noqa: E402
CACHE = vz_paths.vz_dir(OUT, VZ) / "f32_traces_loss5_d13_f209_20_30_e217.npz"   # 184 m/s: tmp/cu2_power/vz184/
FIG = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
POINT, P0, P1 = (13.0, 209.0), 20, 30
VXS, Y0 = [0.2, 0.6], 0.0
STARTS = [(1.5, 2.0, 0), (1.5, 2.0, 1), (1.5, 2.0, 2), (2.5, 2.0, 0)]
NSAVE = 401


def compute():
    from centrex_tlf import lindblad
    from scripts import cu2_power_check as pc
    from scripts.cleanup_spb_doppler import describe
    pc.base.E_FIELD = np.array([0.0, 0.0, 217.0])
    pc.base.VZ = VZ                                       # read by build() and doppler_sign()
    pc.T_END = pc.base.T_END = 0.05 / VZ
    prepared, gidx, _ = pc.setup()
    QN = pc.c2.OBE.QN
    desc = {int(i): describe(QN[i]) for i in gidx}
    lab = {i: (float(d["F1"]), float(d["F"]), int(round(d["mF"]))) for i, d in desc.items()}
    idx = {v: k for k, v in lab.items()}
    jobs = [(idx[s], (POINT, P0, P1, v, Y0)) for s in STARTS for v in VXS]
    dim = len(QN)
    rhos = []
    for i, _ in jobs:
        r = np.zeros((dim, dim), dtype=np.complex128); r[i, i] = 1.0
        rhos.append(prepared.layout.pack(r))
    t_end = pc.T_END
    saveat = np.linspace(0.0, t_end, NSAVE)
    solver = {**pc.SCAN_SOLVER, "maxiters": 1_000_000, "output_when": "saveat"}
    res = lindblad.solve_lindblad_batch(prepared, np.asarray(rhos), (0.0, t_end), parameter_slots=pc.SLOTS,
                                        parameter_batch=np.asarray([pc.slot_values(r) for _, r in jobs],
                                                                   dtype=np.complex128),
                                        saveat=saveat, threads=8, **solver)
    pops = np.asarray(res.values).real.reshape(len(jobs), NSAVE, dim)
    levels = sorted(gidx, key=lambda i: lab[i])
    np.savez(CACHE, t=saveat, z=pc.base.VZ * saveat, pops=pops[:, :, levels].reshape(len(STARTS), len(VXS), NSAVE, -1),
             F1=[lab[i][0] for i in levels], F=[lab[i][1] for i in levels], mF=[lab[i][2] for i in levels],
             starts=np.array(STARTS), vx=np.array(VXS), y0=Y0, vz=pc.base.VZ, first_pass=pc.base.FIRST_PASS,
             pass_spacing=pc.base.PASS_SPACING, sigma_z=pc.base.SIGMA_Z, n_pass=pc.base.N_PASS,
             pass_power=np.array(pc.info["pass_power"]))
    print("wrote", CACHE)


def plot():
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    INK, GRAY, GRID = "#111111", "#555555", "#e6e6e6"
    MCOL = {-2: "#7a3db8", -1: "#2a78d6", 0: "#111111", 1: "#e8673a", 2: "#2e8b57"}
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "mathtext.fontset": "dejavusans",
        "font.size": 8, "axes.linewidth": 0.7, "axes.edgecolor": GRAY, "xtick.color": GRAY, "ytick.color": GRAY,
        "xtick.labelcolor": INK, "ytick.labelcolor": INK, "axes.labelcolor": INK, "axes.spines.top": False,
        "axes.spines.right": False, "savefig.facecolor": "white"})
    z = np.load(CACHE, allow_pickle=True)
    zz = 1e3 * z["z"]
    lab = list(zip(z["F1"], z["F"], z["mF"]))
    sel = {m: lab.index((1.5, 2.0, m)) for m in range(-2, 3)}
    J2 = np.arange(len(lab))
    # laser intensity along z (both beams share the multipass), relative units
    xs = float(z["first_pass"]) + float(z["pass_spacing"]) * np.arange(int(z["n_pass"]))
    prof = sum(p * np.exp(-0.5 * ((z["z"] - x) / float(z["sigma_z"])) ** 2) for p, x in zip(z["pass_power"], xs))
    titles = ["start F₁=3/2, F=2, m$_F$ = 0", "start F₁=3/2, F=2, m$_F$ = +1", "start F₁=3/2, F=2, m$_F$ = +2",
              "start singlet F₁=5/2, F=2, m$_F$ = 0"]
    fig = plt.figure(figsize=(8.6, 9.0))
    gs = fig.add_gridspec(5, 2, height_ratios=[0.35, 1, 1, 1, 1], hspace=0.42, wspace=0.12, left=0.08, right=0.985,
                          top=0.925, bottom=0.10)
    for c, vx in enumerate(z["vx"]):
        axp = fig.add_subplot(gs[0, c])
        axp.fill_between(zz, prof, color="#f2c14e", lw=0)
        axp.set_xlim(0, zz[-1]); axp.set_yticks([]); axp.set_ylim(0, 1.1)
        axp.spines["left"].set_visible(False)
        axp.set_title(f"({'ab'[c]}) v$_x$ = {vx:g} m/s:  Doppler ±{3.68 * vx:.1f} MHz", loc="left", fontsize=8.5,
                      weight="bold", color=INK, pad=4)
        axp.text(1.0, 1.02, "laser intensity, 13 passes", transform=axp.transAxes, fontsize=6.5, color=GRAY,
                 va="bottom", ha="right")
        axp.tick_params(labelbottom=False)
        for r, title in enumerate(titles):
            ax = fig.add_subplot(gs[r + 1, c], sharex=axp)
            P = z["pops"][r, c]
            for m in range(-2, 3):
                ax.plot(zz, P[:, sel[m]], color=MCOL[m], lw=1.3, label=f"m$_F$ = {m:+d}".replace("+0", "0"))
            ax.plot(zz, P[:, sel[-2]:sel[2] + 1].sum(1), color=GRAY, lw=1.0, ls="--", label="F₁=3/2, F=2 total")
            ax.plot(zz, P[:, J2].sum(1), color=GRAY, lw=0.8, ls=":", label="all J = 2")
            for x in xs:
                ax.axvline(1e3 * x, color=GRID, lw=0.6, zorder=0)
            ax.set_yscale("log"); ax.set_ylim(1e-4, 1.5)
            ax.set_xlim(0, zz[-1])
            if c == 0:
                ax.set_ylabel("population")
            else:
                ax.tick_params(labelleft=False)
            ax.text(0.01, 0.04, title, transform=ax.transAxes, fontsize=7.5, color=INK, va="bottom",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
            fin = P[-1, sel[-2]:sel[2] + 1].sum()
            ax.text(0.99, 0.96, f"F₁=3/2, F=2 at the end: {100 * fin:.2f} %", transform=ax.transAxes, fontsize=7,
                    color=INK, ha="right", va="top")
            if r == 3:
                ax.set_xlabel(f"z along CU2 (mm)  =  v$_z$ t,  v$_z$ = {float(z['vz']):g} m/s")
            else:
                ax.tick_params(labelbottom=False)
    h, l = fig.axes[1].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, ncol=7, loc="upper left", bbox_to_anchor=(0.06, 0.975))
    fig.suptitle("CU2 model: m$_F$ populations of X J=2 F₁=3/2, F=2 along the interaction region  "
                 "(217 V/cm, δ₀ 13 MHz, f$_\\mathrm{AOM}$ 209 MHz, 20 + 30 mW, v$_z$ " + f"{float(z['vz']):g} m/s)", x=0.01, y=1.0, ha="left", fontsize=9.5,
                 weight="bold", color=INK)
    fig.text(0.01, 0.045, "y₀ = 0; the Doppler sign alternates pass by pass.  m$_F$ < 0 starts mirror these (m$_F$ → −m$_F$).  Gray lines: pass centres.  "
             "Polarization switches X ↔ Z at 0.5 MHz (≈ 0.37 mm per half period), EOM 17.5 MHz, β = 1.5.\n"
             "Populations are of the field-dressed levels at 217 V/cm, labelled by their zero-field parent "
             "(F₁, F, m$_F$); m$_F$ = 0, ±1 are |m$_J$| = 1, m$_F$ = ±2 are |m$_J$| = 2.",
             fontsize=6.8, color=GRAY, va="top")
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"CU2_F1eq3o2_F2_vs_z{TAG}.{ext}", dpi=200, bbox_inches="tight", pad_inches=0.05)
    print("wrote", FIG / f"CU2_F1eq3o2_F2_vs_z{TAG}.png")


if __name__ == "__main__":
    {"compute": compute, "plot": plot}[sys.argv[1]]()
