"""Overview figure of the CU1 test-condition runs (scripts/cu1_test_conditions.py).

(a) remaining F=2 fraction vs detuning at 50 mW, (b) vs v_x at 15 and 19 MHz, (c) power dependence
against the 10-21 data. Writes tmp/cu1_test_conditions/cu1_test_conditions_overview.png.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import cu1_test_conditions as m  # noqa: E402

OUT = m.OUT
COL = {"old": "#1f6fb4", "test": "#d1495b"}
NAME = {"old": "150 V/cm, 145 m/s", "test": "156 V/cm, 156 m/s (test)"}


def load(cfg, name):
    z = np.load(OUT / cfg / name)
    labs = [tuple(x) for x in z["labels"].tolist()]
    return z, m.f2_left(z["pops"], labs)


fig, ax = plt.subplots(1, 3, figsize=(15, 6.2))

# (a) detuning
a = ax[0]
for cfg in m.CONFIGS:
    z, left = load(cfg, "scan_50mW.npz")
    on, avg = m.means(left)
    D = z["detunings"]
    a.plot(D, avg, "-o", color=COL[cfg], ms=4, label=f"{NAME[cfg]}, y-avg")
    a.plot(D, on, "--s", color=COL[cfg], ms=3, alpha=0.7, label=f"{NAME[cfg]}, on axis")
    a.plot(D, left[:, int(np.argmin(abs(m.VX))), 0], ":", color=COL[cfg], alpha=0.6,
           label=f"{NAME[cfg]}, $v_x=0$ only")
a.axvspan(13, 16, color="0.85", zorder=0, label="measured optimum (4–6 MHz below +19–20)")
a.axvline(19, color="k", lw=0.8)
a.text(19.15, 0.185, "+19 MHz\npredicted\nset point", fontsize=8, va="top")
a.axhline(0.032, color="g", lw=1, ls="-.", label="measured, 50 mW, at optimum (0.032)")
a.set(xlabel="CU1 detuning from main pair (MHz)", ylabel="remaining $F=2$ fraction",
      title="(a) 50 mW, $\\sigma_{v_x}=1.28$ m/s", ylim=(0, 0.19))
a.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)

# (b) v_x at 15 and 19 MHz, test config
b = ax[1]
z, left = load("test", "scan_50mW.npz")
D = list(z["detunings"])
wy = m.y_weights()
for Dk, ls in [(15.0, "-"), (19.0, "--")]:
    k = D.index(Dk)
    b.plot(m.VX, left[k] @ wy, ls, color="k", label=f"{Dk:g} MHz, y-avg")
    b.plot(m.VX, left[k, :, 0], ls, color="0.55", label=f"{Dk:g} MHz, on axis")
b2 = b.twinx()
g = np.exp(-0.5 * (m.VX / m.SIGMA_VX) ** 2)
b2.fill_between(m.VX, g, color="#d1495b", alpha=0.15)
b2.set(ylim=(0, 1.05), yticks=[])
b.text(-3.9, 0.52, "shaded: measured $v_x$ distribution\n($\\sigma=1.28$ m/s)", fontsize=8, color="#a02c3c")
b.set(xlabel="$v_x$ (m/s)", ylabel="remaining $F=2$ fraction", ylim=(0, 0.6),
      title="(b) test config, 50 mW: who is left")
b.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)

# (c) power
c = ax[2]
P = np.array(m.POWERS)
c.plot(P, [m.DATA_POWER[int(p)] for p in P], "ko", ms=7, label="data 10-21 (at measured optimum)", zorder=5)
for cfg in m.CONFIGS:
    for Dk, ls in [(15, "-"), (19, "--")]:
        _, left = load(cfg, f"power_D{Dk}.npz")
        _, avg = m.means(left)
        c.plot(P, avg, ls + "o", color=COL[cfg], ms=4, label=f"{NAME[cfg]}, {Dk} MHz")
c.set(xscale="log", yscale="log", xlabel="CU1 power (mW)", ylabel="remaining $F=2$ fraction (y-avg)",
      title="(c) power dependence", xticks=P, xticklabels=[f"{p:g}" for p in P])
c.xaxis.set_minor_formatter(plt.NullFormatter())
c.set_yticks([0.03, 0.05, 0.1, 0.2], ["0.03", "0.05", "0.1", "0.2"])
c.yaxis.set_minor_formatter(plt.NullFormatter())
c.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1)

fig.suptitle("CU1 model under test conditions: 2 cm beam, 17 passes, 3.5 %/pass loss, AP1 triplet start, "
             "$|v_x|\\leq4$ m/s; y-avg = Gaussian 10.95 mm cut at 8.5 mm (assumed)", fontsize=9.5)
fig.tight_layout()
f = OUT / "cu1_test_conditions_overview.png"
fig.savefig(f, dpi=160)
print(f)
