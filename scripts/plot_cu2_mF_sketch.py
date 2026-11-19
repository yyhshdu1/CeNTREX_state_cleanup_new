"""CU2 mF' level sketch at 217 V/cm in the style of the user's hand drawing (slides; not to scale).

    python scripts/plot_cu2_mF_sketch.py
    -> docs/state_cleanup/figures/slides/CU2_mF_sketch_217.{png,pdf}

Only the P'=- (laser-driven) branch of B J'=3 F1'=7/2 F'=3 and F'=4 is drawn, one bar per mF'. Within each
manifold the bar heights follow the 217 V/cm transition frequencies of scripts/plot_cu2_mF_scheme.py (cached in
tmp/cu2_mF_scheme_217/level_scheme.npz), compressed by U units per MHz; the F'=3 - F'=4 gap is drawn much
smaller than to scale and labelled (~200 MHz; 214.6 MHz between the two main lines at 217 V/cm).
F'=3 beam (blue): carrier Delta_laser = delta0 = 13 MHz above the F'=3 main line (mF'=+-3 bar), EOM sidebands
+-Delta_EOM = 17.5 MHz (the +1 order lands on mF'=+-1). F'=4 beam (red): to mF'=+1, delta1 = +7.4 MHz.
Ground: the X J=2 sublevels at 217 V/cm form two clusters about 1.4 MHz apart (0-0.30 and 1.17-1.74 MHz).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
FIG = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
CACHE = ROOT / "tmp" / "cu2_mF_scheme_217" / "level_scheme.npz"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans"})

BLUE, RED, GREY = "#1010e0", "#e0301e", "0.4"
D0, F_AOM, F_EOM = 13.0, 209.0, 17.5
U = 0.055                      # drawing units per MHz inside a manifold (compressed)
Y3, Y4 = 10.0, 16.0            # F'=3 and F'=4 main lines (gap NOT to scale)
W = 1.0                        # bar width = one mF' column

z = np.load(CACHE)
rows, dmain = z["rows"], float(z["delta_main"])
d1 = D0 + F_AOM - dmain
pos = {k: {int(r[1]): r[3] for r in rows if r[0] == k and r[2] < 0} for k in (0, 1)}   # P' = - only
y3 = {m: Y3 + U * (p - pos[0][3]) for m, p in pos[0].items()}     # F'=3 main line = mF'=+3 bar
y4 = {m: Y4 + U * (p - pos[1][1]) for m, p in pos[1].items()}     # F'=4 main line = mF'=+1 bar

fig, ax = plt.subplots(figsize=(13.0, 6.6))
ax.set_xlim(-8.6, 6.4)
ax.set_ylim(-0.9, 19.2)
ax.axis("off")


def bar(m, y, lw=2.2):
    ax.plot([m - W / 2, m + W / 2], [y, y], color="k", lw=lw, solid_capstyle="butt", zorder=3)


for m, y in y3.items():
    bar(m, y)
for m, y in y4.items():
    bar(m, y)
ax.text(-8.5, Y4 + 0.4, "$|B,\\ J'=3,\\ F_1'=7/2,\\ F'=4,\\ -\\rangle$", fontsize=13, va="center")
ax.text(-8.5, Y3 + 0.4, "$|B,\\ J'=3,\\ F_1'=7/2,\\ F'=3,\\ -\\rangle$", fontsize=13, va="center")
for m in range(-4, 5):
    ax.text(m, max(y4.values()) + 0.55, f"{m:+d}" if m else "0", ha="center", va="bottom", fontsize=10.5,
            color=GREY)
ax.text(-4.9, max(y4.values()) + 0.55, "$m_F'$ =", ha="right", va="bottom", fontsize=11, color=GREY)

# ~200 MHz brace between the two manifolds (gap drawn smaller than to scale)
xb = 5.25
ax.plot([xb - 0.12, xb, xb, xb - 0.12], [Y4 + 0.2, Y4 + 0.2, Y3 - 0.2, Y3 - 0.2], color="0.4", lw=1.3)
ax.plot([xb, xb + 0.15], [(Y3 + Y4) / 2] * 2, color="0.4", lw=1.3)
ax.text(xb + 0.25, (Y3 + Y4) / 2, "~200 MHz", fontsize=13, va="center")

# ground: two clusters of X J=2 sublevels
G0, G1 = 0.0, 0.45
for y in (G0, G1):
    ax.plot([-0.5, 0.5], [y, y], color="k", lw=2.2, solid_capstyle="butt")
ax.plot([-0.62, -0.7, -0.7, -0.62], [G1, G1, G0, G0], color=GREY, lw=1.1)
ax.text(-0.85, (G0 + G1) / 2, "1.4 MHz", fontsize=12, color=GREY, ha="right", va="center")
ax.text(0.65, (G0 + G1) / 2, "$|X,\\ J=2\\rangle$", fontsize=12.5, ha="left", va="center")

# F'=3 beam: carrier (dashed), EOM sidebands, detuning
yc = Y3 + U * D0
ax.add_patch(FancyArrowPatch((0.05, G1 + 0.05), (0.05, yc - 0.05), arrowstyle="-|>", mutation_scale=16,
                             lw=1.8, color=BLUE, zorder=4))
ax.plot([-1.45, 1.45], [yc, yc], color=BLUE, lw=1.4, ls=(0, (4, 3)), zorder=4)
ax.add_patch(FancyArrowPatch((-0.5, yc), (-0.5, yc + U * F_EOM), arrowstyle="-|>", mutation_scale=12, lw=1.5,
                             color=BLUE, zorder=5))
ax.add_patch(FancyArrowPatch((-1.42, yc), (-1.42, yc - U * F_EOM), arrowstyle="-|>", mutation_scale=12,
                             lw=1.5, color=BLUE, zorder=5))
ax.text(-0.4, yc + 0.5 * U * F_EOM + 0.1, "$\\Delta_{EOM}$", fontsize=11, color=GREY, va="center")
ax.text(-1.32, yc - 0.5 * U * F_EOM - 0.1, "$\\Delta_{EOM}$", fontsize=11, color=GREY, va="center")
xl = 1.6                                           # Delta_laser: carrier above the F'=3 main line (mF'=+3)
ax.plot([xl - 0.1, xl, xl, xl - 0.1], [yc, yc, y3[3], y3[3]], color=GREY, lw=1.1)
ax.text(xl + 0.12, (yc + y3[3]) / 2, "$\\Delta_{laser}$", fontsize=11, color=GREY, va="center")
ax.plot([xl - 0.1, 2.5], [y3[3], y3[3]], color=GREY, lw=0.6, ls=":")

# F'=4 beam (made from the F'=3 beam with the AOM)
y4c = Y4 + U * d1
ax.add_patch(FancyArrowPatch((0.1, G1 + 0.05), (1.0, y4c - 0.03), arrowstyle="-|>", mutation_scale=16,
                             lw=1.8, color=RED, zorder=4))
ax.plot([0.55, 1.45], [y4c, y4c], color=RED, lw=1.4, ls=(0, (4, 3)), zorder=4)

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"CU2_mF_sketch_217.{ext}", dpi=170)
print(FIG / "CU2_mF_sketch_217.png", {m: round(v, 2) for m, v in y3.items()}, f"d1 = {d1:.2f}")
