"""Remaining J=2 population of each J=2 level started alone, at two CU2 set points (the `perstate` stage of
scripts/cu2_power_check.py), mean in the detection acceptance on the reduced grid.

    python scripts/plot_cu2_perstate.py          -> docs/state_cleanup/figures/CU2_perstate_217_vs_155.{png,pdf}
    python scripts/plot_cu2_perstate.py 217      -> docs/state_cleanup/figures/CU2_perstate_217.{png,pdf} (217 V/cm only)
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys  # noqa: E402
sys.path.insert(0, str(ROOT))
from scripts import vz_paths  # noqa: E402
OUT = ROOT / "tmp" / "cu2_power"
FIG = ROOT / "docs" / "state_cleanup" / "figures"

import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK, MUTED, GRID = "#262624", "#6b6b68", "#e4e3dc"
RUNS = [("perstate_loss5_d14p5_f209_32_43_e217", "217 V/cm, (14.5, 209) MHz", "#2a78d6"),
        ("perstate_loss5_d12_f212_32_43_e155", "155 V/cm, (12, 212) MHz", "#eb6834")]
mpl.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED,
                     "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK, "axes.labelcolor": INK,
                     "text.color": INK})


def main(only217=False):
    data = []
    for stem, label, color in (RUNS[:1] if only217 else RUNS):
        z = np.load(vz_paths.find(OUT / f"{stem}_result.npz"), allow_pickle=True)
        keys = list(zip(z["F1"], z["F"].astype(int), z["mF"]))
        data.append((label, color, dict(zip(keys, 100 * z["acceptance_remaining"])), dict(zip(keys, z["p_mJ0"]))))
    keys = sorted(data[0][2], key=lambda k: (k[0], k[1], k[2]))
    x = np.arange(len(keys))
    fig, ax = plt.subplots(figsize=(10, 4.2))
    w = 0.6 if only217 else 0.38
    for k, (label, color, rem, _) in enumerate(data):
        ax.bar(x + (0 if only217 else (k - 0.5) * w), [rem[q] for q in keys], w, color=color, label=label, zorder=3)
    for b, q in enumerate(keys):
        if data[0][3][q] > 0.5:
            ax.text(b, -0.06, "m$_J$≈0", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7,
                    color=MUTED)
    ax.set_xticks(x, [f"{q[2]:+d}" for q in keys])
    groups = {}
    for b, q in enumerate(keys):
        groups.setdefault(q[:2], []).append(b)
    for (f1, f), bs in groups.items():
        ax.text(np.mean(bs), -0.15, f"F1={int(2 * f1)}/2, F={f}", transform=ax.get_xaxis_transform(), ha="center",
                va="top", fontsize=8.5)
        if bs[0] > 0:
            ax.axvline(bs[0] - 0.5, color=GRID, lw=1.0, zorder=1)
    ax.set_xlim(-0.6, len(keys) - 0.4)
    ax.set_ylabel("J=2 left after CU2 (%)\nmean in the detection acceptance")
    ax.text(-0.012, -0.02, "m$_F$", transform=ax.transAxes, ha="right", va="top")
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    if not only217:
        ax.legend(frameon=False, loc="upper center", ncol=2)
    ax.set_title("CU2 at 217 V/cm, (14.5, 209) MHz, each J=2 level started alone: 32 + 43 mW, 5 % loss, 3 cm beam"
                 if only217 else "CU2, each J=2 level started alone: 32 + 43 mW, 5 % loss per pass, 3 cm beam; "
                 "24-point acceptance grid", fontsize=9.5)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    stem = "CU2_perstate_217" if only217 else "CU2_perstate_217_vs_155"
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", dpi=200)
    print("wrote", FIG / f"{stem}.png")
    for q in keys:
        print(q, "  ".join(f"{d[2][q]:.2f}" for d in data))


if __name__ == "__main__":
    main(only217="217" in sys.argv[1:])
