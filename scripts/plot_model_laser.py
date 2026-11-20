"""Laser model along the trajectory, for the Model section of docs/state_cleanup/main.tex.

    python scripts/plot_model_laser.py
    -> docs/state_cleanup/figures/model_laser_trajectory.{png,pdf}

Same definitions as scripts/cu2_laser_aom_scan.build: 13 passes 3.5 mm apart, first at 3.5 mm, two-prism optical
order (centrex_tlf.utils.multipass.multipass_prism_order), power fraction (1 - L)^(p - 1) for the p-th pass in
optical order, Gaussian passes sigma_s = 0.75 mm, sigma_y = height / 4, Rabi rate of the main pair
|Omega| = |m| d sqrt(2 I / (c eps0)) / hbar, d = 2.6675506e-30 C m. Doppler sign alternates from one spatial pass
to the next and flips half-way between passes (lindblad.alternating_sign). Polarization gate
g_Z(t) = [1 + sgn sin(2 pi f t + phi)] / 2 at f = 0.5 MHz, phi = 0; the F'=4 beam is in antiphase.
Phase modulation: sideband power J_n(beta)^2, beta = 1.5, f_EOM = 17.5 MHz.
Shown for CU2 at the set point conditions: 32 mW (F'=3, |m| = 0.095) + 43 mW (F'=4, |m| = 0.224), 3 cm beam,
5 % loss per pass, on axis, v_z = 184 m/s.
"""
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np
from scipy.special import jv

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from centrex_tlf import utils  # noqa: E402

FIG = ROOT / "docs" / "state_cleanup" / "figures"
INK, MUTED, GRID = "#262624", "#6b6b68", "#e4e3dc"
C0, C1 = "#2a78d6", "#eb6834"
mpl.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED,
                     "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK, "axes.labelcolor": INK,
                     "text.color": INK})

N_PASS, S1, DS, SS, HEIGHT, LOSS, VZ = 13, 3.5e-3, 3.5e-3, 0.75e-3, 3e-2, 0.05, 184.0
D, C, EPS0, HBAR = 2.6675506e-30, 2.99792458e8, 8.8541878128e-12, 1.054571817e-34
BEAMS = [("$F'=3$ beam, 32 mW", 0.032, 0.095, C0), ("$F'=4$ beam, 43 mW", 0.043, 0.224, C1)]
F_POL, F_EOM, BETA = 0.5e6, 17.5, 1.5

order = np.asarray(utils.multipass.multipass_prism_order(N_PASS))
a = (1 - LOSS) ** (order - 1)
centres = S1 + DS * np.arange(N_PASS)
sy = HEIGHT / 4


def rabi_MHz(s, P, m):
    I = P * np.sum(a[:, None] / (2 * np.pi * SS * sy) * np.exp(-(s[None] - centres[:, None]) ** 2 / (2 * SS**2)), 0)
    return m * D * np.sqrt(2 * I / (C * EPS0)) / HBAR / (2 * np.pi * 1e6)


def sign(s):
    return np.where(np.floor((s - (S1 - DS / 2)) / DS) % 2 == 0, 1, -1)


fig = plt.figure(figsize=(7.2, 4.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], width_ratios=[1.35, 1.0], hspace=0.62, wspace=0.32,
                      left=0.085, right=0.985, top=0.9, bottom=0.11)

# (a) Rabi rate of the main pair along the beam, Doppler sign bands
ax = fig.add_subplot(gs[0, :])
s = np.linspace(0, 0.05, 4001)
for k in range(N_PASS):
    lo, hi = S1 - DS / 2 + k * DS, S1 + DS / 2 + k * DS
    if sign(np.array([centres[k]]))[0] < 0:
        ax.axvspan(1e3 * lo, 1e3 * hi, color="#f1f0ea", lw=0, zorder=0)
for name, P, m, col in BEAMS:
    ax.plot(1e3 * s, rabi_MHz(s, P, m), color=col, lw=1.4, label=name, zorder=3)
top = rabi_MHz(centres, 0.043, 0.224)
for k in range(N_PASS):
    ax.text(1e3 * centres[k], top[k] + 0.04, f"{order[k]}", ha="center", va="bottom", fontsize=7.5, color=MUTED)

ax.set_xlim(0, 50)
ax.set_ylim(0, 1.22)
ax.set_xlabel("position along the molecular beam, $s$ (mm)")
ax.set_ylabel("$|\\Omega_k|/2\\pi$ (MHz)")
ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, loc="upper right", ncol=2, fontsize=8.5, borderaxespad=0.1)
sec = ax.secondary_xaxis("top", functions=(lambda x: x / VZ * 1e3, lambda t: t * VZ / 1e3))
sec.set_xlabel(f"time at {VZ:g} m/s ($\\mu$s)", fontsize=8.5, labelpad=2)
sec.tick_params(labelsize=8)
ax.text(-0.075, 1.02, "(a)", transform=ax.transAxes, fontweight="bold", va="bottom")
ax.text(1e3 * (centres[1]), 0.06, "$-k_Lv_x$", ha="center", fontsize=7.5, color=MUTED)
ax.text(1e3 * (centres[0]), 0.06, "$+k_Lv_x$", ha="center", fontsize=7.5, color=MUTED)

# (b) one pass: envelope and the antiphase polarization gates
ax = fig.add_subplot(gs[1, 0])
t = np.linspace(8e-6, 30e-6, 8001)
gz = 0.5 * (1 + np.where(np.sin(2 * np.pi * F_POL * t) >= 0, 1, -1))
fw = 2.3548 * SS / VZ / 2
ax.axvspan(1e6 * (S1 / VZ - fw), 1e6 * (S1 / VZ + fw), color="#f1f0ea", lw=0, zorder=0)
ax.text(1e6 * S1 / VZ, 2.62, "FWHM of pass 1", ha="center", va="bottom", fontsize=7.5, color=MUTED)
ax.plot(1e6 * t, 1.5 + gz, color=C0, lw=1.2)
ax.plot(1e6 * t, 1 - gz, color=C1, lw=1.2)
ax.text(30.4, 2.0, "$F'=3$\nbeam", color=C0, fontsize=7.5, va="center")
ax.text(30.4, 0.5, "$F'=4$\nbeam", color=C1, fontsize=7.5, va="center")
ax.set_xlim(8, 30)
ax.set_ylim(-0.15, 2.95)
ax.set_yticks([0, 1, 1.5, 2.5])
ax.set_yticklabels(["X", "Z", "X", "Z"])
ax.set_xlabel("time ($\\mu$s)")
ax.set_ylabel("polarization")
ax.spines[["top", "right"]].set_visible(False)
ax.text(-0.17, 1.02, "(b)", transform=ax.transAxes, fontweight="bold", va="bottom")

# (c) phase-modulation spectrum
ax = fig.add_subplot(gs[1, 1])
n = np.arange(-3, 4)
pw = jv(n, BETA) ** 2
ax.vlines(n * F_EOM, 0, pw, color=INK, lw=2.2)
for x, y in zip(n * F_EOM, pw):
    ax.text(x, y + 0.012, f"{y:.3f}", ha="center", va="bottom", fontsize=7)
ax.set_xticks(n * F_EOM, [f"{v:g}".replace("-", "−") for v in n * F_EOM], fontsize=7.5)
ax.set_xlabel("frequency from the carrier (MHz)")
ax.set_ylabel("power fraction $J_n(\\beta)^2$")
ax.set_ylim(0, 0.38)
ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.text(-0.25, 1.02, "(c)", transform=ax.transAxes, fontweight="bold", va="bottom")

FIG.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"model_laser_trajectory.{ext}", dpi=200)
print("wrote", FIG / "model_laser_trajectory.png")
print("pass powers", np.round(a, 3), "sum", a.sum().round(2), "peak Rabi", rabi_MHz(centres[:1], 0.032, 0.095),
      rabi_MHz(centres[:1], 0.043, 0.224))
