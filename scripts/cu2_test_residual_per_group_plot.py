"""Slide figure: J=2 left after CU2 per group (and the J=2 total), measured (line combinations) vs CU2 OBE model.

    python -I scripts/cu2_test_residual_per_group_plot.py
    -> docs/state_cleanup/figures/slides/CU2_test_residual_per_group.{png,pdf}

Numbers from docs/state_cleanup/figures/cu2_residual_lines/numbers.json (scripts/cu2_test_line_combinations.py):
measured = main result (eta_on from the CU2 OBE residual's own m_F and v_x), 1 sigma from the Monte Carlo; model =
CU2 OBE at v_z 156 m/s with the same start split.
"""
import json
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NUM = ROOT / "docs" / "state_cleanup" / "figures" / "cu2_residual_lines" / "numbers.json"
OUT = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
INK, GRAY, GRID = "#111111", "#555555", "#e6e6e6"
MEAS, MODEL = "#1f3a5f", "#d62728"
GROUPS = [("F1=5/2 F=3", "F1=5/2 F=3\n(triplet start)"), ("F1=5/2 F=2", "F1=5/2 F=2\n(singlet start)"),
          ("F1=3/2 F=2", "F1=3/2 F=2\n(pumped in)"), ("F1=3/2 F=1", "F1=3/2 F=1\n(pumped in)"),
          ("total", "J=2 total")]


def main():
    d = json.loads(NUM.read_text())
    main_, model = d["main"], d["model"]
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": GRAY,
                         "xtick.color": GRAY, "ytick.color": GRAY, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                         "savefig.facecolor": "white"})
    fig, ax = plt.subplots(figsize=(6.6, 4.6), constrained_layout=True)
    x = np.arange(len(GROUPS), dtype=float)
    x[-1] += 0.35                                                       # gap before the total
    for k, (key, _) in enumerate(GROUPS):
        m, e = (100 * v for v in main_[key])
        ax.errorbar(x[k] - 0.13, m, e, fmt="o", ms=7, color=MEAS, ecolor=MEAS, elinewidth=1.6, capsize=3, zorder=4)
        ax.annotate(f"{m:.2f}", (x[k] - 0.13, m), xytext=(-7, 0), textcoords="offset points", ha="right",
                    va="center", fontsize=8.5, color=INK)
        mod = 100 * model[key]
        ax.plot(x[k] + 0.13, mod, "D", ms=7, mfc="white", mec=MODEL, mew=1.6, zorder=4)
        ax.annotate(f"{mod:.2f}", (x[k] + 0.13, mod), xytext=(7, 0), textcoords="offset points", ha="left",
                    va="center", fontsize=8.5, color=MODEL)
    ax.axvline((x[3] + x[4]) / 2, color=GRAY, lw=0.7, ls=":")
    ax.axhline(0, color=GRAY, lw=0.7)
    ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
    ax.set_xticks(x, [lab for _, lab in GROUPS], fontsize=9)
    ax.set_xlim(x[0] - 0.6, x[-1] + 0.6)
    ax.set_ylim(-0.7, 1.9)
    ax.set_ylabel("J=2 left after CU2 (% of J=2 before CU2)")
    h = [ax.errorbar([], [], [], fmt="o", ms=6, color=MEAS, capsize=3),
         ax.plot([], [], "D", ms=6, mfc="white", mec=MODEL, mew=1.4)[0]]
    ax.legend(h, ["measured, four detection lines (±1σ)", "CU2 OBE model, v$_z$ 156 m/s"], frameon=False,
              fontsize=8.5, loc="upper left")
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"CU2_test_residual_per_group.{ext}", dpi=220)
    print("wrote", OUT / "CU2_test_residual_per_group.png")


if __name__ == "__main__":
    main()
