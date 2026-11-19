"""CU2 mF' level sketch at 217 V/cm in the style of the user's hand drawing (slides; not to scale).

    python scripts/plot_cu2_mF_sketch.py [delta0 f_AOM]      # default: 217 V/cm simulation set point 14.5 209
    -> docs/state_cleanup/figures/slides/CU2_mF_sketch_217.{png,pdf}        (set point)
       docs/state_cleanup/figures/slides/CU2_mF_sketch_217_d13_f209.{...}  (other delta0 / f_AOM)

Only the P'=- branch of B J'=3 F1'=7/2 F'=3 and F'=4 is drawn, one bar per mF'. Within each manifold the bar
heights are the 217 V/cm transition frequencies on each beam's own detuning axis (main resonance = 0) from
scripts/plot_cu2_mF_scheme.py (cache tmp/cu2_mF_scheme_217/level_scheme.npz), compressed by U units per MHz; the
F'=3 - F'=4 gap is drawn much smaller than to scale (~200 MHz; 214.6 MHz between the two main lines).
Main pairs as in the CU2 simulation (scripts/cu2_laser_aom_scan.build): F'=3 beam X F=3 mF=2 -> F'=3 mF'=3,
F'=4 beam X F=3 mF=0 -> F'=4 mF'=1.
Blue (F'=3 beam): carrier Delta_laser = delta0 above its main line. Red (F'=4 beam, made with the AOM):
carrier Delta'_laser = delta0 + f_AOM - Delta_main above its main line. Both beams carry EOM sidebands
+-Delta_EOM = 17.5 MHz (beta = 1.5).
At the set point (14.5, 209) MHz: blue carrier on mF'=+-2, its +1 order on mF'=+-1; red Delta'_laser = 8.9 MHz,
carrier on mF'=0 (9.6 MHz), its -1 order near mF'=+-2 (-11.2 MHz).
With the extra argument `plus` (-> ..._Pplus.{png,pdf}) the P'=+ parent branch is drawn too (simulation keeps it,
retain_opposite_parity_levels=True). At 217 V/cm every |mF'|>=1 level of both branches is a ~40/60 parity mixture,
so both branches are driven; mF'=0 stays parity-pure and P'=+ mF'=0 is dark from J=2+ (strength 0.001, light grey).
At the set point the red +1 order lands on F'=4 P'=+ mF'=+-1 (24.3 MHz), ~29 % of the red-beam excitation.
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
PLUS = "plus" in sys.argv[1:]  # also draw the P'=+ (opposite-parity parent) branch
NUM = [a for a in sys.argv[1:] if a != "plus"]
D0, F_AOM = (float(NUM[0]), float(NUM[1])) if len(NUM) > 1 else (14.5, 209.0)
SFX = "" if (D0, F_AOM) == (14.5, 209.0) else f"_d{D0:g}_f{F_AOM:g}".replace(".", "p")
SFX += "_Pplus" if PLUS else ""
F_EOM = 17.5
U = 0.055                      # drawing units per MHz inside a manifold (compressed)
Y3, Y4 = (10.0, 18.6) if PLUS else (10.0, 16.0)   # F'=3 and F'=4 main resonances (gap NOT to scale)
W = 1.0                        # bar width = one mF' column

z = np.load(CACHE)
rows, dmain = z["rows"], float(z["delta_main"])
d1 = D0 + F_AOM - dmain
pos = {k: {int(r[1]): r[3] for r in rows if r[0] == k and r[2] < 0} for k in (0, 1)}   # P' = -
posp = {k: {int(r[1]): r[3] for r in rows if r[0] == k and r[2] > 0} for k in (0, 1)}  # P' = +
y3 = {m: Y3 + U * p for m, p in pos[0].items()}
y4 = {m: Y4 + U * p for m, p in pos[1].items()}
y3p = {m: Y3 + U * p for m, p in posp[0].items()}
y4p = {m: Y4 + U * p for m, p in posp[1].items()}
Y_HDR = (max(y4p.values()) + 0.55) if PLUS else (Y4 + U * (d1 + F_EOM) + 0.75)

fig, ax = plt.subplots(figsize=(13.0, 7.6) if PLUS else (13.0, 6.6))
ax.set_xlim(-8.6, 6.6)
ax.set_ylim(-0.9, Y_HDR + 0.9)
ax.axis("off")


def bar(m, y, lw=2.2, c="k"):
    ax.plot([m - W / 2, m + W / 2], [y, y], color=c, lw=lw, solid_capstyle="butt", zorder=3)


def arrow(p0, p1, c, lw=1.5, ms=12):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, lw=lw, color=c, zorder=5))


def bracket(x, ya, yb, label):
    ax.plot([x - 0.1, x, x, x - 0.1], [ya, ya, yb, yb], color=GREY, lw=1.1)
    ax.text(x + 0.12, (ya + yb) / 2, label, fontsize=11, color=GREY, va="center")


for m, y in y3.items():
    bar(m, y)
for m, y in y4.items():
    bar(m, y)
ax.text(-8.5, Y4 + 0.4, "$|B,\\ J'=3,\\ F_1'=7/2,\\ F'=4,\\ -\\rangle$", fontsize=13, va="center")
ax.text(-8.5, Y3 + 0.4, "$|B,\\ J'=3,\\ F_1'=7/2,\\ F'=3,\\ -\\rangle$", fontsize=13, va="center")
if PLUS:                                           # P'=+ parent branch (Stark-mixed Lambda doublet partner)
    for yy in (y3p, y4p):
        for m, y in yy.items():          # mF'=0 stays pure P=+ (no Stark mixing): dark from J=2+, light grey
            bar(m, y, lw=1.6 if m == 0 else 2.2, c="0.75" if m == 0 else "k")
    for F, yy in ((3, y3p), (4, y4p)):
        ax.text(-8.5, np.mean(list(yy.values())) + 0.15,
                f"$|B,\\ J'=3,\\ F_1'=7/2,\\ F'={F},\\ +\\rangle$", fontsize=13, va="center")
for m in range(-4, 5):
    ax.text(m, Y_HDR, f"{m:+d}" if m else "0", ha="center", va="bottom", fontsize=10.5, color=GREY)
ax.text(-4.9, Y_HDR, "$m_F'$ =", ha="right", va="bottom", fontsize=11, color=GREY)

# ~200 MHz brace between the two manifolds (gap drawn smaller than to scale)
xb = 5.6
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

# F'=3 beam (blue): carrier, EOM sidebands, detuning from its main line (mF'=+3 bar)
yc = Y3 + U * D0
arrow((-0.25, G1 + 0.05), (-0.25, yc - 0.03), BLUE, lw=1.8, ms=16)
ax.plot([-1.45, 1.45], [yc, yc], color=BLUE, lw=1.4, ls=(0, (4, 3)), zorder=4)
arrow((-0.5, yc), (-0.5, yc + U * F_EOM), BLUE)
arrow((-1.42, yc), (-1.42, yc - U * F_EOM), BLUE)
ax.text(-0.4, yc + 0.5 * U * F_EOM + 0.12, "$\\Delta_{EOM}$", fontsize=11, color=GREY, va="center")
ax.text(-1.32, yc - 0.5 * U * F_EOM - 0.12, "$\\Delta_{EOM}$", fontsize=11, color=GREY, va="center")
ax.plot([1.55, 2.6], [y3[3]] * 2, color=GREY, lw=0.6, ls=":")
bracket(2.75, yc, y3[3], "$\\Delta_{laser}$")

# F'=4 beam (red, made from the F'=3 beam with the AOM): carrier, EOM sidebands, detuning from mF'=+1 bar
y4c = Y4 + U * d1
arrow((0.25, G1 + 0.05), (0.25, y4c - 0.03), RED, lw=1.8, ms=16)
ax.plot([-0.6, 2.6], [y4c, y4c], color=RED, lw=1.4, ls=(0, (4, 3)), zorder=4)
arrow((0.62, y4c), (0.62, y4c + U * F_EOM), RED)
arrow((1.6, y4c), (1.6, y4c - U * F_EOM), RED)
ax.text(0.72, y4c + 0.5 * U * F_EOM + 0.12, "$\\Delta_{EOM}$", fontsize=11, color=GREY, va="center")
ax.text(1.7, y4c - 0.5 * U * F_EOM - 0.12, "$\\Delta_{EOM}$", fontsize=11, color=GREY, va="center")
ax.plot([1.55, 2.6], [y4[1]] * 2, color=GREY, lw=0.6, ls=":")
bracket(2.75, y4c, y4[1], "$\\Delta_{laser}'$")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"CU2_mF_sketch_217{SFX}.{ext}", dpi=170)
print(FIG / f"CU2_mF_sketch_217{SFX}.png", f"delta0 = {D0:g}, d1 = {d1:.2f} MHz")
print("F'=3 bars", {m: round(p, 1) for m, p in sorted(pos[0].items())})
print("F'=4 bars", {m: round(p, 1) for m, p in sorted(pos[1].items())})
