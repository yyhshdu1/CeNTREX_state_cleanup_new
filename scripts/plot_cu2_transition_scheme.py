"""Main CU2 transitions: X J=2+ -> B J'=3- F1'=7/2 F'=3 and F'=4 (R(2)), slide schematic.

    python scripts/plot_cu2_transition_scheme.py
    -> docs/state_cleanup/figures/slides/CU2_transition_scheme_R2_F3_F4.{png,pdf}

Ground J levels are single bars (hyperfine not drawn, rotational spacings not to scale). The excited levels are
drawn to scale in MHz:
  F'=4 - F'=3 = 179.9 MHz at E = 0 (centrex B Hamiltonian, J'=1-5, both parities; CU2_test_line_frequencies.md).
  Shaded band = spread of the B mF' sublevels at 217 V/cm about their mean (F'=3: -19.5..+24.5 MHz,
  F'=4: -20.8..+24.1 MHz, same Hamiltonian).
  Comb = EOM sidebands, f_EOM = 17.5 MHz, beta = 1.5; tick length ~ J_n(beta)^2 (n = -2..2), centred at the
  carrier detuning delta from the main pair (CU2 test: delta0 = 13 MHz, f_AOM = 209 MHz,
  delta1 = delta0 + f_AOM - Delta_main = 13 + 209 - 214.6 = 7.4 MHz; Delta_main from CU2_test_line_frequencies.md).
Branching at 217 V/cm from `scripts/cu2_decay_branching.py 217` (range over |mF'| >= 1 of F'=3 and F'=4):
J=2+ 0.21-0.32, J=3- 0.38-0.59, J=4+ 0.18-0.29, J=5- 0.012-0.019, J=1- <= 2.6e-5; mF'=0: 0.53 / 0.47 to J=2+ / 4+.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.special import jv

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
                     "mathtext.fontset": "stix"})

RED, BLUE, GREY = "#d1495b", "#1f6fb4", "0.6"
MHZ = 1 / 50.0                                     # drawing units per MHz in the excited manifold
Y_F3 = 9.0
Y_F4 = Y_F3 + 179.9 * MHZ
BAND = {Y_F3: (-19.5, 24.5), Y_F4: (-20.8, 24.1)}
F_EOM, BETA = 17.5, 1.5
D0, F_AOM, D_MAIN = 13.0, 209.0, 214.551
D1 = D0 + F_AOM - D_MAIN

fig, ax = plt.subplots(figsize=(10.8, 7.6))
ax.set_xlim(-1.2, 13.4)
ax.set_ylim(-0.9, 15.8)
ax.axis("off")


def level(x0, x1, y, lw=3.0, c="k"):
    ax.plot([x0, x1], [y, y], color=c, lw=lw, solid_capstyle="butt", zorder=3)


def wavy(x0, y0, x1, y1, c="k", lw=1.1, n=14, amp=0.09):
    t = np.linspace(0, 1, 400)
    dx, dy = x1 - x0, y1 - y0
    L = np.hypot(dx, dy)
    nx, ny = -dy / L, dx / L
    s = amp * np.sin(2 * np.pi * n * t) * np.clip(np.minimum(t, 1 - t) * 12, 0, 1)
    ax.plot(x0 + dx * t + nx * s, y0 + dy * t + ny * s, color=c, lw=lw, zorder=2)
    ax.annotate("", (x1, y1), (x1 - 0.02 * dx, y1 - 0.02 * dy),
                arrowprops=dict(arrowstyle="-|>", color=c, lw=lw, mutation_scale=11))


# ---- ground X levels (schematic) -------------------------------------------------------------------------
X_LEV = {2: (0.0, 4.6, 0.0, "k", 3.2), 3: (5.2, 6.8, 1.3, "k", 3.2), 4: (7.4, 9.4, 2.8, "k", 3.2),
         5: (10.0, 11.6, 4.5, "k", 3.2)}
X_LAB = {2: "$|J=2^{+}\\rangle$", 3: "$|J=3^{-}\\rangle$", 4: "$|J=4^{+}\\rangle$", 5: "$|J=5^{-}\\rangle$"}
for J, (x0, x1, y, c, lw) in X_LEV.items():
    level(x0, x1, y, lw=lw, c=c)
    ax.text((x0 + x1) / 2, y - 0.25, X_LAB[J], ha="center", va="top", fontsize=14, color=c)
ax.text(13.2, 0.0, "$X\\,^1\\Sigma^+$", ha="right", va="center", fontsize=16)

# ---- excited B levels (to scale in MHz) ------------------------------------------------------------------
for y, lab in ((Y_F3, "$F'=3$"), (Y_F4, "$F'=4$")):
    lo, hi = BAND[y]
    ax.add_patch(Rectangle((0.3, y + lo * MHZ), 8.2, (hi - lo) * MHZ, fc="0.9", ec="none", zorder=1))
    level(0.3, 8.5, y)
    ax.text(8.7, y, lab, ha="left", va="center", fontsize=15)
ax.text(4.4, Y_F4 + 0.75, "$B\\,^3\\Pi_1$,  $|J'=3^{-},\\,F_1'=7/2\\rangle$", ha="center", va="bottom",
        fontsize=15)
# hyperfine splitting
ax.annotate("", (10.0, Y_F4), (10.0, Y_F3), arrowprops=dict(arrowstyle="<->", color="k", lw=1.0))
ax.text(10.15, (Y_F3 + Y_F4) / 2, "179.9 MHz\n($E=0$)", ha="left", va="center", fontsize=11.5)

# ---- laser beams ------------------------------------------------------------------------------------------
pn = {n: jv(n, BETA) ** 2 for n in range(-2, 3)}
for x, y_lev, d, col, name in ((1.4, Y_F3, D0, RED, "$F'=3$ beam"), (3.0, Y_F4, D1, BLUE, "$F'=4$ beam")):
    ax.add_patch(FancyArrowPatch((x, 0.05), (x, y_lev), arrowstyle="-|>", mutation_scale=22, lw=3.2,
                                 color=col, zorder=4))
    yc = y_lev + d * MHZ
    for n, p in pn.items():                            # EOM comb next to the arrow tip
        w = 1.6 * p
        ax.plot([x - 0.45 - w, x - 0.45], [yc + n * F_EOM * MHZ] * 2, color=col, lw=2.0, zorder=5)
    ax.text(x + 0.2, 4.6, name, rotation=90, ha="left", va="center", fontsize=12.5, color=col)
# the F'=4 beam is made from the F'=3 beam by the AOM
ax.add_patch(FancyArrowPatch((1.48, 1.7), (2.92, 1.7), arrowstyle="-|>", mutation_scale=16, lw=1.6,
                             color="0.2", connectionstyle="arc3,rad=-0.35", zorder=6))
ax.text(2.2, 2.25, "AOM", ha="center", va="bottom", fontsize=12, fontweight="bold", color="0.2")
ax.text(2.2, 1.45, "+$f_{AOM}$\n(209 MHz)", ha="center", va="top", fontsize=10, color="0.2", linespacing=1.2)

# ---- decays (branching at 217 V/cm, scripts/cu2_decay_branching.py 217) ---------------------------------
for (x0, x1, y1), lab, (tx, ty) in (((4.4, 3.9, 0.12), "0.21–0.32", (4.45, 4.3)),
                                    ((5.7, 6.0, 1.42), "0.38–0.59", (6.05, 5.2)),
                                    ((7.2, 8.4, 2.92), "0.18–0.29", (7.95, 6.1)),
                                    ((8.3, 10.6, 4.62), "0.01–0.02", (9.75, 6.9))):
    wavy(x0, Y_F3 - 0.05, x1, y1)
    ax.text(tx, ty, lab, fontsize=12.5)

# ---- title and beam settings ------------------------------------------------------------------------------
ax.text(-1.1, 15.75, "CU2: R(2) $F_1'=7/2$, $F'=3$ and $F'=4$ at $E$ = 217 V/cm", ha="left", va="top",
        fontsize=16, fontweight="bold")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"CU2_transition_scheme_R2_F3_F4.{ext}", dpi=170)
print(FIG / "CU2_transition_scheme_R2_F3_F4.png", f"delta1 = {D1:.2f} MHz")
