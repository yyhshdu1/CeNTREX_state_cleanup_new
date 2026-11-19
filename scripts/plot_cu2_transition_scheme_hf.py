"""CU2 transitions with the full B J'=3- hyperfine structure, no decays (slide schematic).

    python scripts/plot_cu2_transition_scheme_hf.py
    -> docs/state_cleanup/figures/slides/CU2_transition_scheme_R2_hf.{png,pdf}

Companion of scripts/plot_cu2_transition_scheme.py (same beams, EOM comb and AOM marker). The four hyperfine
levels of B J'=3- are drawn to scale at their E = 0 positions relative to F1'=7/2 F'=3 (centrex B Hamiltonian,
J'=1-5, both parities): F1'=7/2 F'=4 +179.9 MHz, F1'=5/2 F'=2 +577.2 MHz, F1'=5/2 F'=3 +921.5 MHz.
The shaded bands (driven F1'=7/2 levels only; the F1'=5/2 levels also shift by ~+60 MHz) are the mF' sublevel ranges at 217 V/cm on the same axis (relative to the 217 V/cm mean of
F1'=7/2 F'=3): -19.5..+24.5, +159.9..+204.9, +601.0..+665.5, +946.9..+1015.0 MHz.
Driven levels (F1'=7/2) are black, the F1'=5/2 levels grey. Ground J=2+ is one bar (no hyperfine).
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.special import jv

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
                     "mathtext.fontset": "stix"})

RED, BLUE, GREY = "#d1495b", "#1f6fb4", "0.55"
MHZ = 1 / 95.0                                     # drawing units per MHz in the excited manifold
Y0 = 5.0                                          # F1'=7/2 F'=3
# (F1', F', E=0 position MHz, 217 V/cm band MHz, driven)
LEVELS = [("7/2", 3, 0.0, (-19.5, 24.5), True), ("7/2", 4, 179.9, (159.9, 204.9), True),
          ("5/2", 2, 577.2, (601.0, 665.5), False), ("5/2", 3, 921.5, (946.9, 1015.0), False)]
F_EOM, BETA = 17.5, 1.5
D0, F_AOM, D_MAIN = 13.0, 209.0, 214.551
D1 = D0 + F_AOM - D_MAIN

fig, ax = plt.subplots(figsize=(9.0, 8.0))
ax.set_xlim(-0.6, 11.0)
ax.set_ylim(-0.9, Y0 + 921.5 * MHZ + 2.0)
ax.axis("off")

# ground level
ax.plot([0.0, 4.6], [0, 0], color="k", lw=3.2, solid_capstyle="butt")
ax.text(2.3, -0.25, "$|J=2^{+}\\rangle$", ha="center", va="top", fontsize=14)
ax.text(10.8, 0.0, "$X\\,^1\\Sigma^+$", ha="right", va="center", fontsize=16)

# excited hyperfine levels
y_of = {}
for f1, f, pos, (lo, hi), driven in LEVELS:
    y = Y0 + pos * MHZ
    y_of[(f1, f)] = y
    c = "k" if driven else GREY
    if driven:                                     # 217 V/cm mF' spread, driven levels only
        ax.add_patch(Rectangle((0.3, Y0 + lo * MHZ), 6.4, (hi - lo) * MHZ, fc="0.9", ec="none", zorder=1))
    ax.plot([0.3, 6.7], [y, y], color=c, lw=3.0 if driven else 2.0, solid_capstyle="butt", zorder=3)
    ax.text(6.9, y, f"$F_1'={f1},\\ F'={f}$", ha="left", va="center", fontsize=13.5, color=c)
    if pos > 0:
        ax.text(9.0, y, f"+{pos:.1f} MHz", ha="left", va="center", fontsize=11.5, color=c)
ax.text(3.5, Y0 + 921.5 * MHZ + 0.45, "$B\\,^3\\Pi_1$,  $|J'=3^{-}\\rangle$", ha="center", va="bottom",
        fontsize=15)

# laser beams with EOM combs
pn = {n: jv(n, BETA) ** 2 for n in range(-2, 3)}
for x, key, d, col, name in ((1.4, ("7/2", 3), D0, RED, "$F'=3$ beam"), (3.0, ("7/2", 4), D1, BLUE, "$F'=4$ beam")):
    y_lev = y_of[key]
    ax.add_patch(FancyArrowPatch((x, 0.05), (x, y_lev), arrowstyle="-|>", mutation_scale=22, lw=3.2,
                                 color=col, zorder=4))
    yc = y_lev + d * MHZ
    for n, p in pn.items():
        w = 1.6 * p
        ax.plot([x - 0.45 - w, x - 0.45], [yc + n * F_EOM * MHZ] * 2, color=col, lw=1.8, zorder=5)
    ax.text(x + 0.2, 3.5, name, rotation=90, ha="left", va="center", fontsize=12.5, color=col)

# the F'=4 beam is made from the F'=3 beam by the AOM
ax.add_patch(FancyArrowPatch((1.48, 1.5), (2.92, 1.5), arrowstyle="-|>", mutation_scale=16, lw=1.6,
                             color="0.2", connectionstyle="arc3,rad=-0.35", zorder=6))
ax.text(2.2, 2.05, "AOM", ha="center", va="bottom", fontsize=12, fontweight="bold", color="0.2")
ax.text(2.2, 1.25, "+$f_{AOM}$\n(209 MHz)", ha="center", va="top", fontsize=10, color="0.2", linespacing=1.2)

ax.text(-0.5, Y0 + 921.5 * MHZ + 1.95, "CU2: R(2) $F_1'=7/2$, $F'=3$ and $F'=4$ at $E$ = 217 V/cm",
        ha="left", va="top", fontsize=16, fontweight="bold")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"CU2_transition_scheme_R2_hf.{ext}", dpi=170)
print(FIG / "CU2_transition_scheme_R2_hf.png")
