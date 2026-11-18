"""CU2 cleanup of the SPB start (lens levels) vs the laser-modulation parameters, from stored proxy runs (no new
simulation): f_AOM, EOM frequency, EOM modulation depth beta and the X/Z polarization-switching frequency.

    python scripts/plot_cu2_param_insensitive.py
    -> docs/state_cleanup/figures/slides/CU2_param_insensitive_spb.{png,pdf}

All panels use the six-trajectory proxy of scripts/cu2_power_check.py at v_z 184 m/s, sigma_vx 0.65 m/s, 5 %/pass,
3 cm beam, and the same y axis, to show that none of these needs a dedicated optimization. The runs differ:
  f_AOM:  217 V/cm, 32 + 43 mW, delta0 = 13 MHz (the CU2-test detuning)      [fine217/opt217/setpoints runs]
  f_EOM:  217 V/cm, 32 + 43 mW, (14.5, 209) MHz                              [stage eomfreq]
  beta:   200 V/cm, 35 + 25 mW, (12.5, 207.5) and (13, 206) MHz              [stage beta]
  pol:    200 V/cm, 35 + 25 mW, (12.5, 207.5) MHz, X/Z square-wave switching [stage polfreq]
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import plot_cu2_setpoint_grid as g  # noqa: E402
import plot_cu2_aom_range as A  # noqa: E402  (also writes its own figure; harmless)

import matplotlib.pyplot as plt  # noqa: E402

FIG = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
EOM_FREQS = [13.0, 14.0, 15.0, 15.5, 16.0, 17.0, 17.5, 18.5]
BETA_POINTS, BETAS = [(12.5, 207.5), (13.0, 206.0)], [1.0, 1.25, 1.5, 1.75, 2.0, 2.5]
POL_FREQS = [0.1, 0.25, 0.5, 1.0, 2.0]
RED, BLUE, GREEN = "#d1495b", "#1f6fb4", "#6a9a3a"

f_aom, t_aom = A.curves("spb")[13.0]
t_eom = g.load("eomfreq_loss5_d14p5_f209_32_43_e217_spb", len(EOM_FREQS))
EOM_RANGE = (15.0, 18.0)                                   # MHz the EOM can actually be tuned over
keep = (np.array(EOM_FREQS) >= EOM_RANGE[0]) & (np.array(EOM_FREQS) <= EOM_RANGE[1])
EOM_FREQS, t_eom = list(np.array(EOM_FREQS)[keep]), t_eom[keep]
t_beta = g.load("beta_loss5_35_25_spb", len(BETA_POINTS) * len(BETAS)).reshape(len(BETA_POINTS), len(BETAS))
t_pol = g.load("polfreq_loss5_35_25_spb", len(POL_FREQS))
for name, v in (("f_AOM", t_aom), ("f_EOM", t_eom), ("beta", t_beta), ("pol", t_pol)):
    print(name, np.round(v, 2))

fig, ax = plt.subplots(1, 4, figsize=(16, 4.6), sharey=True)
panels = [
    (ax[0], "(a) AOM frequency", "$f_{AOM}$ (MHz)", 209, "217 V/cm, 32 + 43 mW, $\\delta_0$ = 13 MHz",
     [(f_aom, t_aom, RED, "$\\delta_0$ = 13 MHz (CU2 test)")]),
    (ax[1], "(b) EOM frequency", "$f_{EOM}$ (MHz)", 17.5, "217 V/cm, 32 + 43 mW, (14.5, 209) MHz",
     [(np.array(EOM_FREQS), t_eom, BLUE, "(14.5, 209) MHz")]),
    (ax[2], "(c) EOM modulation depth", "$\\beta$ (rad)", 1.5, "200 V/cm, 35 + 25 mW",
     [(np.array(BETAS), t_beta[0], BLUE, "(12.5, 207.5) MHz"), (np.array(BETAS), t_beta[1], GREEN, "(13, 206) MHz")]),
    (ax[3], "(d) X/Z polarization switching", "switching frequency (MHz)", 0.5, "200 V/cm, 35 + 25 mW, (12.5, 207.5) MHz",
     [(np.array(POL_FREQS), t_pol, BLUE, "(12.5, 207.5) MHz")]),
]
for a, title, xl, used, cond, series in panels:
    a.axhspan(98.9 - 0.6, 98.9 + 0.6, color=RED, alpha=0.10, zorder=0)
    a.axvline(used, color="k", lw=0.8)
    for x, y, col, lab in series:
        a.plot(x, y, "-o", color=col, lw=1.8, ms=4.5,
               label=f"{lab}: spread {y.max() - y.min():.2f} %")
    a.set_title(f"{title}\n{cond}", fontsize=9.5)
    a.set_xlabel(xl)
    a.legend(fontsize=7.5, loc="lower left", bbox_to_anchor=(0, 0.06))
    a.text(used, 97.05, f"  used: {used:g}", fontsize=8, va="bottom")
ax[1].set_xlim(EOM_RANGE[0] - 0.2, EOM_RANGE[1] + 0.2)
ax[3].set_xscale("log")
ax[3].set_xticks(POL_FREQS, [f"{f:g}" for f in POL_FREQS])
ax[0].set_ylabel("cleanup of the SPB start (%)")
ax[0].set_ylim(97.0, 100.0)
ax[0].text(185, 98.38, "CU2 test 98.9(6) %", color="#a02c3c", fontsize=8.5)
fig.suptitle("CU2 cleanup of the SPB start (lens levels) vs the laser-modulation parameters: "
             "near the values we use, the changes stay below the ±0.6 % CU2-test error", fontsize=11.5)
fig.text(0.5, 0.005, "Stored six-trajectory proxy runs (no new simulation), $v_z$ 184 m/s, $\\sigma_{v_x}$ 0.65 m/s, "
         "5 %/pass, 3 cm beam. Black line: value used. Red band: CU2 test 98.9(6) %. Spread = max − min over the plotted range.", ha="center", va="bottom",
         fontsize=8.5, color="0.35")
fig.tight_layout(rect=(0, 0.04, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"CU2_param_insensitive_spb.{ext}", dpi=160)
print(FIG / "CU2_param_insensitive_spb.png")
