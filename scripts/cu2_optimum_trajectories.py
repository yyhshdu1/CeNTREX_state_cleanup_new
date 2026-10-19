"""Population trajectories at the two CU2 optima at 200 V/cm (regions R1 and R2 of the 2D scan).

    python scripts/cu2_optimum_trajectories.py compute   # ~5 min, cached in tmp/cu2_laser_aom_scan/
    python scripts/cu2_optimum_trajectories.py plot

Writes docs/state_cleanup/figures/CU2_200V_optimum_populations.{pdf,png}  (rotational manifolds)
and   docs/state_cleanup/figures/CU2_200V_optimum_J2_sublevels.{pdf,png}  (the 20 J=2 sublevels).
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cu2_laser_aom_scan import (CHECKPOINT_DIR, CONFIRM_SOLVER, FIRST_PASS, N_PASS,  # noqa: E402
                                        PASS_SPACING, T_END, VZ, build)

CACHE = CHECKPOINT_DIR / "optimum_trajectories.npz"
OUT = ROOT / "docs" / "state_cleanup" / "figures"
POINTS = {"R1": (11.5, 208.5), "R2": (11.5, 191.0)}   # (delta0/2pi, f_AOM) in MHz
SAVE_STEP = 1e-7                                       # s


def compute():
    from centrex_tlf import hamiltonian, lindblad, states
    t0 = time.perf_counter()
    prepared, rho0, ground_idx, info, obe = build(return_obe=True)
    QN = obe.QN
    X, B = states.ElectronicState.X, states.ElectronicState.B
    excited = np.asarray(states.QuantumSelector(electronic=B).get_indices(QN), dtype=int).ravel()
    other_J = {J: np.asarray(states.QuantumSelector(J=J, electronic=X).get_indices(QN), dtype=int).ravel()
               for J in (1, 3, 4, 5)}
    labels = np.array([(QN[i].largest.F1, QN[i].largest.F, QN[i].largest.mF) for i in ground_idx])
    pts = np.array(list(POINTS.values()))
    d1 = pts[:, 0] + pts[:, 1] - info["delta_main_MHz"]
    batch = 2 * np.pi * 1e6 * np.column_stack((pts[:, 0], d1))
    print(f"model ready ({time.perf_counter() - t0:.0f} s); delta1 = {d1}", flush=True)

    saveat = np.append(np.arange(0.0, T_END, SAVE_STEP), T_END)
    solver = {**CONFIRM_SOLVER, "output": "populations", "output_when": "saveat", "saveat": saveat}
    t0 = time.perf_counter()
    res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=["detuning0", "detuning1"],
                                  parameter_batch=batch, threads=len(pts), **solver)
    pops = np.asarray(res.values).real
    print(f"trajectories {time.perf_counter() - t0:.0f} s, shape {pops.shape}", flush=True)

    gamma = float(getattr(hamiltonian, "Γ"))
    solver = {**CONFIRM_SOLVER, "output": "photon_integral", "output_when": "final",
              "integral_weights": [(int(i), gamma) for i in excited]}
    t0 = time.perf_counter()
    ph = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=["detuning0", "detuning1"],
                                 parameter_batch=batch, threads=len(pts), **solver)
    photons = np.asarray(ph.values).real.reshape(len(pts), -1)[:, -1]
    print(f"photon integral {time.perf_counter() - t0:.0f} s: {photons}", flush=True)

    np.savez(CACHE, names=np.array(list(POINTS)), points=pts, detuning1=d1, t=saveat, pops=pops,
             photons=photons, ground_idx=ground_idx, excited=excited, labels=labels,
             **{f"J{J}": v for J, v in other_J.items()}, delta_main=info["delta_main_MHz"])
    for k, name in enumerate(POINTS):
        fin = pops[k, -1]
        print(f"{name}: cleanup {100*(1-fin[ground_idx].sum()):.4f}%, photons {photons[k]:.3f}, "
              + ", ".join(f"J={J} {fin[v].sum():.4f}" for J, v in other_J.items()), flush=True)


def plot():
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    INK, MUTED, GRID, SURFACE = "#262624", "#6b6b68", "#e4e4e1", "#ffffff"
    SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]  # categorical, fixed order
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm", "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
        "axes.labelcolor": INK, "text.color": INK, "savefig.facecolor": SURFACE, "legend.frameon": False,
    })
    z = np.load(CACHE)
    names, pts, d1 = z["names"], z["points"], z["detuning1"]
    s_mm = z["t"] * VZ * 1e3
    pops, photons, gidx, exc = z["pops"], z["photons"], z["ground_idx"], z["excited"]
    labels = z["labels"]
    pass_mm = (FIRST_PASS + PASS_SPACING * np.arange(N_PASS)) * 1e3

    def frame(ax):
        for x in pass_mm:
            ax.axvline(x, color=GRID, lw=0.6, zorder=0)
        ax.set_xlim(0, 50); ax.set_yscale("log")
        ax.set_xlabel("distance along the molecular beam (mm)")
        ax.tick_params(length=2.5, width=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    def title(ax, k):
        rem = pops[k, -1, gidx].sum()
        ax.set_title(f"{names[k]}: $f_\\mathrm{{AOM}}={pts[k,1]:.1f}$ MHz, $\\delta_0/2\\pi={pts[k,0]:.1f}$ MHz\n"
                     f"$\\eta_\\mathrm{{c}}={100*(1-rem):.2f}$ %, $N_\\gamma={photons[k]:.2f}$",
                     loc="left", color=INK)

    # ---- Figure A: rotational manifolds ----------------------------------------------------------
    series = [("$X$, $J=2$ (addressed)", gidx, SLOTS[0], 1.8),
              ("$X$, $J=3$", z["J3"], SLOTS[1], 1.4),
              ("$X$, $J=4$", z["J4"], SLOTS[2], 1.4),
              ("$X$, $J=5$", z["J5"], SLOTS[3], 1.4),
              ("$B$ (excited)", exc, SLOTS[4], 0.8),
              ("$X$, $J=1$", z["J1"], SLOTS[5], 1.4)]
    fig, axes = plt.subplots(1, 2, figsize=(6.1, 3.3), sharey=True, constrained_layout=True)
    for k, ax in enumerate(axes):
        frame(ax)
        for name, idx, color, lw in series:
            y = pops[k][:, idx].sum(axis=1)
            ax.plot(s_mm, np.clip(y, 1e-9, None), color=color, lw=lw, label=name)
        # selective direct labels: final value of the addressed level and the two main destinations
        for (name, idx, color, lw), dy in zip(series[:3], (2.2, 1.35, 0.72)):
            yf = pops[k, -1, idx].sum()
            ax.annotate(f"{yf:.3g}" if yf < 0.01 else f"{yf:.3f}", xy=(50, yf), xytext=(49, yf * dy),
                        ha="right", va="center", fontsize=6.5, color=INK)
        ax.set_ylim(1e-5, 1.6)
        title(ax, k)
    axes[0].set_ylabel("population")
    handles = [Line2D([], [], color=c, lw=max(lw, 1.4)) for _, _, c, lw in series]
    fig.legend(handles, [s[0] for s in series], loc="outside lower center", ncol=6, fontsize=7,
               handlelength=1.6, columnspacing=1.1)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"CU2_200V_optimum_populations.{ext}", dpi=220)
    plt.close(fig)

    # ---- Figure B: the 20 J=2 sublevels, grouped by hyperfine level ---------------------------------
    groups = [(1.5, 1.0), (1.5, 2.0), (2.5, 2.0), (2.5, 3.0)]
    gname = {g: f"$F_1={'3/2' if g[0] == 1.5 else '5/2'}$, $F={int(g[1])}$" for g in groups}
    fig, axes = plt.subplots(1, 2, figsize=(6.1, 3.3), sharey=True, constrained_layout=True)
    for k, ax in enumerate(axes):
        frame(ax)
        for col, (F1, F, mF) in zip(gidx, labels):
            c = SLOTS[groups.index((F1, F))]
            ax.plot(s_mm, np.clip(pops[k][:, col], 1e-12, None), color=c, lw=0.9)
        # label the three sublevels with the most population left
        fin = pops[k, -1, gidx]
        order = np.argsort(fin)[::-1][:3]
        # stacked in the empty upper-right corner, leader lines to the line ends
        for r, ytxt in zip(order, (2.6e-2, 1.25e-2, 6.0e-3)):
            yf = fin[r]
            F1, F, mF = labels[r]
            ax.annotate(f"{gname[(F1, F)]}, $m_F={int(mF):+d}$" if mF else f"{gname[(F1, F)]}, $m_F=0$",
                        xy=(50, yf), xytext=(46.5, ytxt), ha="right", va="center", fontsize=6.5,
                        color=INK, arrowprops=dict(arrowstyle="-", lw=0.5, color=MUTED, shrinkA=1, shrinkB=1,
                                                   relpos=(1.0, 0.5)))
        ax.set_ylim(1e-6, 0.08)
        title(ax, k)
    axes[0].set_ylabel("population of each $J=2$ sublevel")
    handles = [Line2D([], [], color=SLOTS[i], lw=1.6) for i in range(len(groups))]
    fig.legend(handles, [f"$J=2$, {gname[g]}" for g in groups], loc="outside lower center", ncol=4,
               fontsize=7, handlelength=1.6, columnspacing=1.4)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"CU2_200V_optimum_J2_sublevels.{ext}", dpi=220)
    plt.close(fig)
    print("written:", *sorted(p.name for p in OUT.glob("CU2_200V_optimum_*")))


if __name__ == "__main__":
    {"compute": compute, "plot": plot}[sys.argv[1]]()
