"""CU2 sensitivity to the AOM frequency, from the 217 V/cm set-point grids of scripts/cu2_power_check.py.

    python scripts/plot_cu2_aom_sensitivity.py
    -> docs/state_cleanup/figures/slides/CU2_aom_sensitivity_217.{pdf,png}

(a), (b): cleanup (six-trajectory proxy) vs f_AOM for each delta0 of the fine grid, SPB source levels and a uniform
start; (c): the same grid read along delta0 at f_AOM = 209 MHz, on the same scale of change, for comparison.
Grids: 32 + 43 mW, 5 %/pass loss, 3 cm beam, v_z 184 m/s, sigma_vx 0.65 m/s (the stored set-point runs).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import plot_cu2_setpoint_grid as g  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

FIG = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
TITLE = {"spb": "SPB source levels (lens, $m_J\\approx0$)", "uniform": "uniform $J=2$ start"}


def grid(name, init):
    d0s, aoms = g.grid_def(name)
    t = g.load(f"{name}_loss5_32_43_e217_{init}", len(d0s) * len(aoms))
    return np.array(d0s), np.array(aoms), t.reshape(len(d0s), len(aoms))


fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
cmap = plt.get_cmap("viridis")
for k, init in enumerate(["spb", "uniform"]):
    a = ax[k]
    d0s, aoms, t = grid("fine217", init)
    for i, d in enumerate(d0s):
        bold = d == 14.5
        a.plot(aoms, t[i], "-o", color=cmap(i / (len(d0s) - 1)), lw=2.4 if bold else 1.2, ms=5 if bold else 3.5,
               label=f"$\\delta_0$ = {d:g} MHz" + ("  (set point)" if bold else ""))
    a.axvline(209, color="k", lw=0.8)
    a.axvspan(205, 211, color="0.92", zorder=0)
    a.set(xlabel="$f_{AOM}$ (MHz)", ylabel="cleanup (%)", title=f"({'ab'[k]}) {TITLE[init]}: vs $f_{{AOM}}$",
          xticks=aoms)
    if init == "spb":
        a.axhspan(98.9 - 0.6, 98.9 + 0.6, color="#d1495b", alpha=0.12, zorder=0)
        a.axhline(98.9, color="#d1495b", lw=1, ls="--")
        a.text(207.2, 98.45, "CU2 test 98.9(6) % (group-solved; 20 + 30 mW, 156 m/s)", color="#a02c3c",
               fontsize=8, va="center")
        a.set_ylim(97.4, 100.0)
    else:
        a.set_ylim(95.0, 99.5)
    a.legend(fontsize=7, loc="lower left", ncol=2)

# (c) delta0 vs f_AOM: change from the set point along each axis
c = ax[2]
for init, col in [("spb", "#1f6fb4"), ("uniform", "#e08a00")]:
    d0s, aoms, t = grid("fine217", init)
    i0, j0 = list(d0s).index(14.5), list(aoms).index(209.0)
    c.plot(aoms - 209, t[i0] - t[i0, j0], "-o", color=col, label=f"{TITLE[init]}: move $f_{{AOM}}$")
    c.plot(d0s - 14.5, t[:, j0] - t[i0, j0], "--s", color=col, label=f"{TITLE[init]}: move $\\delta_0$")
c.axhline(0, color="k", lw=0.6)
c.axhspan(-0.6, 0.6, color="#d1495b", alpha=0.10, zorder=0)
c.text(-3.4, 0.42, "±0.6 %: error of the CU2 test", color="#a02c3c", fontsize=8)
c.set(xlabel="shift from the set point (14.5 MHz, 209 MHz) (MHz)", ylabel="change of cleanup (%)",
      title="(c) $f_{AOM}$ vs $\\delta_0$: same shift, different cost", ylim=(-3.8, 0.8))
c.legend(fontsize=7, loc="lower left")

fig.suptitle("CU2 at 217 V/cm, 32 + 43 mW, 5 %/pass, 3 cm beam, $v_z$ 184 m/s, $\\sigma_{v_x}$ 0.65 m/s "
             "(six-trajectory proxy). Grey band: 205–211 MHz. Test was at (13, 209 MHz), 20 + 30 mW, 156 m/s.",
             fontsize=9.5)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"CU2_aom_sensitivity_217.{ext}", dpi=160)
print(FIG / "CU2_aom_sensitivity_217.png")
