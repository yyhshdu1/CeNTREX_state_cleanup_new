"""Slide schematic: where CU1 sits in the CeNTREX beamline and why it is needed.

    python scripts/plot_beamline_cu1_schematic.py

Writes docs/state_cleanup/figures/slides/beamline_cu1_schematic.{png,pdf}. Slide style (not for the report).
"""
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
NAVY, INK, MUTED, RED = "#1F3A5F", "#222222", "#6b6b68", "#c0392b"
FILL, CU_FILL = "#eef2f7", "#fde3df"
mpl.rcParams.update({"font.family": ["Calibri", "DejaVu Sans"], "font.size": 13})

# (name, condition, what happens, state after)
STAGES = [
    ("Rotational\ncooling", r"$\mathcal{E}\approx 0$", "laser + microwaves\npump J = 1, 2 → J = 0", "J = 0"),
    ("State prep. A\n(SPA)", r"$\mathcal{E}\approx 100$ V/cm", "AP1: J = 0 → 1\nAP2: J = 1 → 2", "J = 2, $m_J$ = 0"),
    ("Electrostatic\nlens", "± 30 kV", "focuses\nJ = 2, $m_J$ = 0", "J = 2, $m_J$ = 0"),
    ("CU1", "SPB entrance", "R(1) $F_1'$=5/2, $F'$=2\npumps off J = 1", "J = 2 only"),
    ("State prep. B\n(SPB)", r"$\mathcal{B}\approx 0$", "microwave AP\nJ = 2, $m_J$=0 → J = 1",
     "J = 1 target"),
    ("CU2", "SPB exit", "R(2) $F_1'$=7/2, $F'$=3, 4\npumps off J = 2", "J = 1 target"),
    ("Interaction\nregion", "", "Ramsey measurement\nin J = 1", ""),
]
WIDTH = [1.55, 1.55, 1.55, 1.45, 1.6, 1.45, 1.65]


def main():
    fig = plt.figure(figsize=(13.0, 5.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 13.0)
    ax.set_ylim(0, 5.6)
    ax.set_axis_off()

    # molecular beam arrow behind the boxes
    ax.add_patch(FancyArrowPatch((0.15, 3.55), (12.85, 3.55), arrowstyle="-|>", mutation_scale=22, lw=2.2,
                                 color="#9aa7b5", zorder=0))

    gap, x = 0.18, 0.35
    centres, edges = [], []
    for (name, cond, action, state), w in zip(STAGES, WIDTH):
        cu1 = name == "CU1"
        box = FancyBboxPatch((x, 2.75), w, 1.6, boxstyle="round,pad=0.02,rounding_size=0.12",
                             facecolor=CU_FILL if cu1 else FILL, edgecolor=RED if cu1 else "#7f8c9a",
                             lw=3.0 if cu1 else 1.3, zorder=2)
        ax.add_patch(box)
        cx = x + w / 2
        centres.append(cx)
        edges.append((x, x + w))
        ax.text(cx, 4.02, name, ha="center", va="center", fontsize=14.5, fontweight="bold",
                color=RED if cu1 else NAVY, zorder=3, linespacing=1.0)
        ax.text(cx, 3.25, action, ha="center", va="center", fontsize=11.5, color=INK, zorder=3, linespacing=1.2)
        if state:                                                       # state after the stage, below the box
            ax.text(cx, 2.45, state, ha="center", va="center", fontsize=12, color=NAVY, zorder=3)
        x += w + gap
    for a, b in zip(centres[:-1], centres[1:]):                        # state flow arrows below the boxes
        ax.add_patch(FancyArrowPatch((a + 0.45, 2.45), (b - 0.45, 2.45), arrowstyle="-|>", mutation_scale=12,
                                     lw=1.0, color="#9aa7b5", zorder=1))
    ax.text(0.35, 2.12, "state after each stage", fontsize=10, color=MUTED, ha="left", va="center")

    # goal line above
    ax.text(0.35, 4.80, "Goal of the cleanup: pump off J = 1 before SPB (CU1) and J = 2 after SPB (CU2)",
            fontsize=15, color=INK, ha="left", va="center")

    # why CU1
    why = FancyBboxPatch((0.35, 0.18), 12.3, 1.62, boxstyle="round,pad=0.02,rounding_size=0.10",
                         facecolor="white", edgecolor=RED, lw=1.6)
    ax.add_patch(why)
    ax.text(0.6, 1.52, "Why CU1?", fontsize=15, fontweight="bold", color=RED, ha="left", va="center")
    lines = [
        "• J = 1 molecules still reach SPB: thermal J = 1 left by rotational cooling, and J = 1, $m_J$ = 0 left by an "
        "incomplete AP2 (J = 1 → 2).",
        "• In SPB the adiabatic passage would transfer them back, J = 1 → J = 2, or pass to the interaction region mixed with the "
        "prepared molecules",
        "   → extra noise and systematics downstream. CU1 removes all of J = 1 before SPB.",
    ]
    for k, t in enumerate(lines):
        ax.text(0.6, 1.12 - 0.36 * k, t, fontsize=13, color=INK, ha="left", va="center")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"beamline_cu1_schematic.{ext}", dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", OUT / "beamline_cu1_schematic.png")


if __name__ == "__main__":
    main()
