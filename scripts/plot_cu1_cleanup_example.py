"""Slide figure: one CU1 trajectory from the stored run (no new simulation).

    python scripts/plot_cu1_cleanup_example.py

Reads run trajectories_cu1_vx0_19p0_20p2MHz_4p8cm of results/SPB_J1_cleanup_scans.h5 (one trajectory per starting
J=1 sublevel, on the beam axis), averages over the 12 starting sublevels (uniform start) at +19.0 MHz, and plots
(a) the J=1 sublevel populations and (b) the rotational levels the molecules end in, versus distance.
Writes docs/state_cleanup/figures/slides/cu1_cleanup_example.{png,pdf}.
"""
import json
from pathlib import Path

import h5py
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
RUN = "trajectories_cu1_vx0_19p0_20p2MHz_4p8cm"
mpl.rcParams.update({"font.family": ["Calibri", "DejaVu Sans"], "font.size": 14, "axes.spines.top": False,
                     "axes.spines.right": False})
COL = {(0.5, 0): "#4b3d9c", (0.5, 1): "#1baf7a", (1.5, 1): "#eb6834", (1.5, 2): "#2a78d6"}
NAME = {(0.5, 0): "$F_1$=1/2, F=0", (0.5, 1): "$F_1$=1/2, F=1", (1.5, 1): "$F_1$=3/2, F=1", (1.5, 2): "$F_1$=3/2, F=2"}
STYLE = {0: "-", 1: "--", 2: ":"}


def main():
    with h5py.File(ROOT / "results" / "SPB_J1_cleanup_scans.h5", "r") as h:
        g = h["runs"][RUN]
        res = g["results"][0].real                      # +19.0 MHz: (start level, time, level)
        t = g["axes"]["time_us"][:]
        fp = json.loads(g.attrs["fixed_parameters_json"])
    labels = eval(fp["level_labels"]) if isinstance(fp["level_labels"], str) else fp["level_labels"]
    vz = fp["molecular_vz_m_per_s"]
    pop = res.mean(axis=0)                              # uniform start: equal weight on the 12 starting sublevels
    s = t * 1e-6 * vz * 1e3                             # distance in mm

    j1 = [k for k, l in enumerate(labels) if l.startswith("X J=1")]
    fig, (a, b) = plt.subplots(1, 2, figsize=(14.0, 4.8), gridspec_kw={"width_ratios": [1.35, 1]})
    for ax in (a, b):
        for n in range(1, fp["passes"] + 1):
            ax.axvline(3.5 * n, color="#dddddd", lw=0.8, zorder=0)       # pass centres: 3.5 mm + 3.5 mm steps
    done = set()
    for k in j1:
        f1, f, mf = (part.split("=")[1] for part in labels[k].split(", ")[1:])
        key = (float(f1), int(f))
        m = abs(int(mf))
        tag = (key, m)
        lab = None if tag in done else f"{NAME[key]}, |$m_F$|={m}"
        done.add(tag)
        a.plot(s, pop[:, k], color=COL[key], ls=STYLE[m], lw=1.6, label=lab)
    a.set_ylim(0, 0.09)
    a.set_xlim(0, 50)
    a.set_xlabel("distance along the CU1 region (mm)")
    a.set_ylabel("population")
    a.set_title("(a) the 12 J = 1 sublevels (start 1/12 each)", loc="left", fontsize=15)
    a.legend(fontsize=11, frameon=False, loc="upper right", handlelength=2.4)

    rot = {"J = 2": "X J=2 (compacted)", "J = 3": "X J=3 (compacted)", "J = 4": "X J=4 (compacted)",
           "J = 0": "X J=0 (compacted)"}
    cols = {"J = 2": "#eb6834", "J = 3": "#2a78d6", "J = 4": "#1baf7a", "J = 0": "#4b3d9c"}
    for name, lab in rot.items():
        k = labels.index(lab)
        b.plot(s, pop[:, k], color=cols[name], lw=2.0, label=f"{name}: {pop[-1, k]:.3f}" if pop[-1, k] > 1e-3
               else f"{name}: {pop[-1, k]:.1e}")
    b.plot(s, pop[:, j1].sum(axis=1), color="black", lw=2.6, label=f"J = 1: {pop[-1, j1].sum():.3f}")
    exc = [k for k, l in enumerate(labels) if l.startswith("B ")]
    b.plot(s, pop[:, exc].sum(axis=1), color="#888888", lw=1.0, label="excited B")
    b.set_yscale("log")
    b.set_ylim(1e-6, 1.5)
    b.set_xlim(0, 50)
    b.set_xlabel("distance along the CU1 region (mm)")
    b.set_title("(b) total per J (final value in legend)", loc="left", fontsize=15)
    b.legend(fontsize=11.5, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"cu1_cleanup_example.{ext}", dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", OUT / "cu1_cleanup_example.png", "| left in J=1:", pop[-1, j1].sum())


if __name__ == "__main__":
    main()
