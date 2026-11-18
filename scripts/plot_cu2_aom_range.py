"""CU2 cleanup vs f_AOM over the full 185-215 MHz range (200 +- 15 MHz) at 217 V/cm, for three delta0.

    python scripts/plot_cu2_aom_range.py [delta0 ...]     # e.g. 13 -> only the CU2-test detuning, file suffix _d13
    -> docs/state_cleanup/figures/slides/CU2_aom_range_217_spb.{pdf,png}   (SPB start only; the earlier
       CU2_aom_range_217 with the uniform start as well is kept)

Merges three runs of scripts/cu2_power_check.py, all at 32 + 43 mW, 5 %/pass loss, 3 cm beam, v_z 184 m/s,
sigma_vx 0.65 m/s, six-trajectory proxy (plot_cu2_setpoint_grid.proxy):
  fine217 grid (205-213 MHz, 1-2 MHz steps), opt217 coarse grid (190-215 MHz, 5 MHz steps, delta0 = 12.5 only)
  and the "setpoints" fill-in run (185-200 and 215 MHz at delta0 = 13, 14.5; 185 MHz at 12.5).
delta0 = 14.5 MHz is the 217 V/cm set point, 13 MHz the CU2 test, 12.5 MHz the earlier set point.
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
FILL = ("d14p5_f185+d14p5_f190+d14p5_f195+d14p5_f200+d14p5_f215+d13_f185+d13_f190+d13_f195+d13_f200+d13_f205+"
        "d13_f209+d13_f215+d12p5_f185")
FILL_PTS = [(float(s.split("_f")[0][1:].replace("p", ".")), float(s.split("_f")[1])) for s in FILL.split("+")]
D0S = {14.5: ("#1f6fb4", "set point 217 V/cm"), 13.0: ("#d1495b", "CU2 test"), 12.5: ("#6a9a3a", "earlier set point")}


def curves(init):
    pts = {}
    for name in ("opt217", "fine217"):
        d0s, aoms = g.grid_def(name)
        t = g.load(f"{name}_loss5_32_43_e217_{init}", len(d0s) * len(aoms)).reshape(len(d0s), len(aoms))
        for i, d in enumerate(d0s):
            for j, f in enumerate(aoms):
                pts[(d, f)] = t[i, j]
    t = g.load(f"setpoints_loss5_32_43_e217_{init}_{FILL}", len(FILL_PTS))
    pts.update(dict(zip(FILL_PTS, t)))
    out = {}
    for d in D0S:
        fs = sorted(f for (dd, f) in pts if dd == d)
        out[d] = np.array(fs), np.array([pts[(d, f)] for f in fs])
    return out


# SPB start (lens levels) only; the uniform start is not drawn. Optional argv: delta0 values to draw (e.g. 13)
ONLY = [float(v) for v in sys.argv[1:]] or list(D0S)
SFX = "" if len(sys.argv) < 2 else "_d" + "_".join(f"{v:g}".replace(".", "p") for v in ONLY)
fig, a = plt.subplots(figsize=(10.5, 4.8))
for d, (f, t) in curves("spb").items():
    if d not in ONLY:
        continue
    col, what = D0S[d]
    a.plot(f, t, "-o", color=col, lw=2.2 if d == 14.5 else 1.4, ms=5 if d == 14.5 else 4,
           label=f"$\\delta_0$ = {d:g} MHz ({what}): {t.min():.2f}–{t.max():.2f} %, spread {t.max() - t.min():.2f} %")
    print("spb", d, dict(zip(f.tolist(), np.round(t, 3).tolist())))
a.axvspan(205, 211, color="0.92", zorder=0)
a.axvline(209, color="k", lw=0.8)
a.text(209.3, 98.62, " 209 MHz\n (test)", fontsize=8, va="bottom")
a.axhspan(98.9 - 0.6, 98.9 + 0.6, color="#d1495b", alpha=0.10, zorder=0)
a.text(184, 98.66, "CU2 test 98.9(6) % (group-solved; 20 + 30 mW, 156 m/s)", color="#a02c3c", fontsize=8.5)
a.set(xlabel="$f_{AOM}$ (MHz)", ylabel="cleanup of the SPB start (%)", xticks=np.arange(185, 216, 5), xlim=(183, 217),
      ylim=(98.6, 100.0))
a.set_title("CU2 cleanup of the SPB start (lens levels) vs $f_{AOM}$ = 200 ± 15 MHz, 217 V/cm", fontsize=11)
a.legend(fontsize=8.5 if len(ONLY) == 1 else 8, loc="upper left")
fig.text(0.5, 0.005, "32 + 43 mW, 5 %/pass, 3 cm beam, $v_z$ 184 m/s, $\\sigma_{v_x}$ 0.65 m/s (six-trajectory proxy). "
         "Grey band: 205–211 MHz. $f_{AOM}$ moves only the $F'=4$ beam.", ha="center", va="bottom", fontsize=8.5,
         color="0.35")
fig.tight_layout(rect=(0, 0.04, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"CU2_aom_range_217_spb{SFX}.{ext}", dpi=160)
print(FIG / f"CU2_aom_range_217_spb{SFX}.png")
