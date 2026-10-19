"""CU1 single trajectories on the beam axis (vx = 0, y0 = 0): J=1 sublevel populations versus time.

    python scripts/cu1_trajectories.py compute     # one trajectory per starting J=1 sublevel, 19.0 and 20.2 MHz
    python scripts/cu1_trajectories.py plot        # docs/state_cleanup/figures/CU1_19MHz_J1_sublevels.{pdf,png}

The OBE is linear in rho0, so the per-start trajectories give any incoherent initial mixture; the figure uses
the 8 SPB-related levels (|m_J|=1 content > 0.5, equal weights) and the uniform 12-level start.
Model: scripts/cleanup_spb_doppler.build_cu1 (80 mW, 4.8 cm beam, 13 passes, 10 % loss, X/Z at 0.5 MHz).
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cleanup_spb_doppler import CU1_SOLVER, build_cu1, describe  # noqa: E402
from scripts.cu2_laser_aom_scan import FIRST_PASS, N_PASS, PASS_SPACING, T_END, VZ, log  # noqa: E402

CACHE = ROOT / "tmp" / "cu1_trajectories" / "cu1_trajectories.npz"
OUT = ROOT / "docs" / "state_cleanup" / "figures"
DETUNINGS = (19.0, 20.2)            # MHz: chosen operating point, previous on-axis optimum
SAVE_STEP = 0.1e-6                  # s


def compute():
    from centrex_tlf import hamiltonian, lindblad, states
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    prepared, gidx, obe = build_cu1(return_obe=True)
    QN, dim = obe.QN, len(obe.QN)
    X, B = states.ElectronicState.X, states.ElectronicState.B
    excited = np.asarray(states.QuantumSelector(electronic=B).get_indices(QN), dtype=int).ravel()
    other = {J: np.asarray(states.QuantumSelector(J=J, electronic=X).get_indices(QN), dtype=int).ravel()
             for J in (0, 2, 3, 4)}
    desc = [describe(QN[i]) for i in gidx]
    labels = np.array([(d["F1"], d["F"], d["mF"]) for d in desc])
    spb = np.array([d["p_mj"].get(1, 0) + d["p_mj"].get(-1, 0) > 0.5 for d in desc])
    log(f"model ready ({time.perf_counter() - t0:.0f} s, {dim} levels)")

    saveat = np.append(np.arange(0.0, T_END, SAVE_STEP), T_END)
    solver = {**CU1_SOLVER, "output": "populations", "output_when": "saveat", "saveat": saveat}
    rho, rows = [], []
    for d in DETUNINGS:
        for i in gidx:
            r = np.zeros((dim, dim), dtype=np.complex128); r[i, i] = 1.0
            rho.append(prepared.layout.pack(r)); rows.append([2 * np.pi * 1e6 * d])
    t0 = time.perf_counter()
    res = lindblad.solve_lindblad_batch(prepared, np.asarray(rho), (0.0, T_END), parameter_slots=["laser_detuning"],
                                        parameter_batch=np.asarray(rows, dtype=np.complex128), threads=12, **solver)
    pops = np.asarray(res.values).real.reshape(len(DETUNINGS), len(gidx), len(saveat), dim)
    log(f"{len(rows)} trajectories in {time.perf_counter() - t0:.0f} s; trace drift "
        f"{np.abs(pops.sum(axis=-1) - 1).max():.1e}")

    # photons per molecule for the two initial mixtures, from the running integral
    gamma = float(getattr(hamiltonian, "Γ"))
    ph_solver = {**CU1_SOLVER, "output": "photon_integral", "output_when": "final",
                 "integral_weights": [(int(i), gamma) for i in excited]}
    photons = {}
    for name, sel in (("spb", spb), ("uniform", np.ones(len(gidx), bool))):
        rho0 = np.zeros((dim, dim), dtype=np.complex128)
        rho0[gidx[sel], gidx[sel]] = 1.0 / sel.sum()
        ph = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=["laser_detuning"],
                                     parameter_batch=np.array([[2 * np.pi * 1e6 * d] for d in DETUNINGS],
                                                              dtype=np.complex128), threads=2, **ph_solver)
        photons[name] = np.asarray(ph.values).real.reshape(len(DETUNINGS), -1)[:, -1]
    log(f"photons per molecule: {photons}")
    np.savez(CACHE, detunings=np.array(DETUNINGS), t=saveat, pops=pops, ground_idx=gidx, excited=excited,
             labels=labels, spb=spb, photons_spb=photons["spb"], photons_uniform=photons["uniform"],
             **{f"J{J}": v for J, v in other.items()})
    fin = pops[:, :, -1, :][:, :, gidx].sum(axis=-1)
    for k, d in enumerate(DETUNINGS):
        log(f"{d} MHz: cleanup SPB {100 * (1 - fin[k, spb].mean()):.3f} %, uniform {100 * (1 - fin[k].mean()):.3f} %")
        for j in np.argsort(-fin[k]):
            F1, F, mF = labels[j]
            log(f"   start F1={F1} F={F:.0f} mF={mF:+.0f}{' (SPB)' if spb[j] else ''}: left {fin[k, j]:.2e}")


def plot():
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from scripts.plot_cu1_doppler_map import GRID, INK, MUTED, SERIES  # noqa: F401  (also sets rcParams)

    z = np.load(CACHE)
    t, pops, gidx, labels, spb = z["t"], z["pops"], z["ground_idx"], z["labels"], z["spb"]
    dets = list(z["detunings"])
    k19, k20 = dets.index(19.0), dets.index(20.2)
    x = VZ * t * 1e3                                  # position along the interaction region (mm)
    passes = (FIRST_PASS + PASS_SPACING * np.arange(N_PASS)) * 1e3
    fam_color = {(1.5, 2): SERIES[0], (1.5, 1): SERIES[1], (0.5, 1): SERIES[2], (0.5, 0): "#4a3aa7"}
    fam_name = {(1.5, 2): "$F_1=3/2,\\ F=2$", (1.5, 1): "$F_1=3/2,\\ F=1$", (0.5, 1): "$F_1=1/2,\\ F=1$",
                (0.5, 0): "$F_1=1/2,\\ F=0$"}
    dash = {0: "-", 1: (0, (5, 2)), 2: (0, (1.5, 1.5))}

    def mixture(k, sel):
        return pops[k][sel].mean(axis=0)              # (time, level)

    fig = plt.figure(figsize=(6.3, 8.8))
    gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 0.22, 0.8], hspace=0.38)
    axes = [fig.add_subplot(gs[0])]
    axes += [fig.add_subplot(gs[1], sharex=axes[0]), None, fig.add_subplot(gs[3], sharex=axes[0])]
    leg_ax = fig.add_subplot(gs[2]); leg_ax.axis("off")
    starts = [("spb", spb, "(a) start: 8 SPB-related levels, $\\delta/2\\pi=19.0$ MHz",
               z["photons_spb"][k19]),
              ("uniform", np.ones(len(gidx), bool), "(b) start: all 12 $J=1$ sublevels, $\\delta/2\\pi=19.0$ MHz",
               z["photons_uniform"][k19])]
    for ax, (key, sel, title, nph) in zip(axes[:2], starts):
        p19, p20 = mixture(k19, sel), mixture(k20, sel)
        done = set()
        for j, (F1, F, mF) in enumerate(labels):
            fam, m = (F1, int(F)), int(abs(mF))
            if (fam, m) in done:
                continue
            done.add((fam, m))
            idx = [gidx[q] for q, (a, b, c) in enumerate(labels) if (a, int(b)) == fam and abs(c) == m]
            y = p19[:, idx].mean(axis=1)              # one curve per +-mF pair (equal by symmetry)
            ax.plot(x, y, color=fam_color[fam], ls=dash[m], lw=1.3,
                    label=f"{fam_name[fam]}, $|m_F|={m}$" + (" (each)" if m else ""))
        ax.plot(x, p19[:, gidx].sum(axis=1), color=INK, lw=2.0, label="total $J=1$")
        ax.plot(x, p20[:, gidx].sum(axis=1), color=MUTED, lw=1.4, ls=(0, (4, 2)),
                label="total $J=1$ at $\\delta/2\\pi=20.2$ MHz")
        fin19, fin20 = p19[-1, gidx].sum(), p20[-1, gidx].sum()
        ax.set_yscale("log"); ax.set_ylim(1e-5, 1.5)
        ax.set_ylabel("population")
        ax.set_title(f"{title}\nleft in $J=1$: {fin19:.2e} (20.2 MHz: {fin20:.2e});  "
                     f"photons per molecule {nph:.2f}", loc="left", fontsize=8.5)
        for p in passes:
            ax.axvline(p, color=GRID, lw=0.8, zorder=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    handles, labs = axes[0].get_legend_handles_labels()
    leg_ax.legend(handles, labs, loc="center", fontsize=6.5, frameon=False, ncol=3, handlelength=2.4,
                  columnspacing=1.2)
    for a in axes[:2]:
        a.tick_params(labelbottom=False)

    ax = axes[3]
    p19 = mixture(k19, spb)
    dest = [("$J=3$", z["J3"], SERIES[0]), ("$J=2$", z["J2"], SERIES[1]), ("$J=4$", z["J4"], SERIES[2]),
            ("$J=0$", z["J0"], "#4a3aa7"), ("excited ($B$)", z["excited"], INK)]
    for name, idx, color in dest:
        ax.plot(x, p19[:, idx].sum(axis=1), color=color, lw=1.4,
                label=f"{name}: {p19[-1, idx].sum():.3g}" if name[0] == "$" else name)
    ax.set_yscale("log"); ax.set_ylim(1e-7, 2.0)
    ax.set_ylabel("population")
    ax.set_xlabel("position along the interaction region $v_z t$ (mm)")
    ax.set_title("(c) where the molecules go, start as in (a), $\\delta/2\\pi=19.0$ MHz (final values in legend)",
                 loc="left", fontsize=8.5)
    for p in passes:
        ax.axvline(p, color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(loc="lower center", fontsize=6.5, frameon=False, ncol=3, borderaxespad=0.2)
    fig.text(0.01, 0.005, "Single trajectory on the beam axis, $v_x=0$, $y_0=0$; 80 mW, 4.8 cm beam, 13 passes "
             "(grey lines: pass centres), 10 % loss per pass, X/Z switching at 0.5 MHz. Populations sampled every "
             "0.1 µs.", fontsize=6.5, ha="left", va="top", color=INK, wrap=True)
    fig.savefig(OUT / "CU1_19MHz_J1_sublevels.pdf", bbox_inches="tight")
    fig.savefig(OUT / "CU1_19MHz_J1_sublevels.png", dpi=200, bbox_inches="tight")
    log("wrote CU1_19MHz_J1_sublevels.{pdf,png}")


if __name__ == "__main__":
    {"compute": compute, "plot": plot}[sys.argv[1]]()
