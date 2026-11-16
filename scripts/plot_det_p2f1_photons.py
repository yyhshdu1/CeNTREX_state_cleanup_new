"""P(2) F1'=1/2 F'=1 probe: OBE photons per molecule of each J=2 group vs DET power (top panel of
scripts/det_power_saturation.py only, landscape).

    python -I scripts/plot_det_p2f1_photons.py
    -> docs/state_cleanup/figures/slides/DET_power_P2F1o2F1_photons.{png,pdf}

Reads the cached probe OBE scans in tmp/detection/power_scan (DET plates off, E = 0, earth field 0.5 G, sigma+ about
k, flat top 1.5 x 1.5 cm, v_z 156 m/s, lens v_x distribution); same group definitions as det_power_saturation.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import det_power_saturation as S  # noqa: E402

import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LINE = "P1o2F1"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 10.5,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": "white"})
S.run_missing([LINE])
eta = np.array(S.group_yields([LINE])[LINE]["eta"])
P = np.array(S.POWERS) * 1e3
k24 = S.POWERS.index(24e-3)

fig, a = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
for j, (g, _) in enumerate(S.groups_of(LINE)):
    if eta[:, j].max() < 0.02:
        continue
    main = g == S.MAIN[LINE]
    a.plot(P, eta[:, j], "-" if main else "--", color=S.GCOL[g], lw=2.4 if main else 1.4, marker="o",
           ms=4 if main else 3, label=S.gname(g) + (" (main)" if main else "") + f":  {eta[k24, j]:.3f} at 24 mW")
a.axvline(24, color="#999", lw=0.9)
a.set_xscale("log")
tk = [0.3, 1, 3, 10, 24, 100]
a.set_xticks(tk, [f"{x:g}" for x in tk])
a.set(xlabel="DET power (mW)", ylabel="photons per molecule η", ylim=(0, None))
a.legend(frameon=False, fontsize=9, loc="upper left")
a.set_title("P(2) F1'=1/2 F'=1: photons per molecule vs DET power (plates off, σ+)", fontsize=11.5)
for ext in ("png", "pdf"):
    fig.savefig(S.OUT / f"DET_power_P2F1o2F1_photons.{ext}", dpi=200)
print(S.OUT / "DET_power_P2F1o2F1_photons.png")
print("eta at 24 mW", dict(zip([S.gname(g) for g, _ in S.groups_of(LINE)], np.round(eta[k24], 4))))
