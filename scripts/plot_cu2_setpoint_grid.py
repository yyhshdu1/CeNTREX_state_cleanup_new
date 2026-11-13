"""Proxy cleanup efficiency on a (delta0, f_AOM) set-point grid of scripts/cu2_power_check.py, one panel per start.

    python scripts/plot_cu2_setpoint_grid.py opt217 32 43 e217
    -> docs/state_cleanup/figures/CU2_setpoint_grid_opt217_32_43_e217.{pdf,png}

Proxy = 6 points (vx = -0.6, 0, 0.6 m/s with Gaussian weights, sigma 0.65 m/s; y0 = 0, 6 mm averaged), as in
cu2_power_check.proxy. Only complete grids are plotted (a start whose grid is unfinished is left out).
"""
import sys
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

PROXY_VX, PROXY_Y0 = [-0.6, 0.0, 0.6], [0.0, 6e-3]
INK, MUTED = "#262624", "#6b6b68"
LABEL = {"spb": "SPB start (4 levels, $m_J\\approx 0$)", "uniform": "uniform start (20 J=2 levels)"}
DELTA_MAIN_MHZ = {"155": 204.520, "200": 211.829, "217": 214.585}   # F'=4 minus F'=3 main line
mpl.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
                     "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "axes.labelcolor": INK, "text.color": INK})


def grid_def(name):
    """Read GRIDS from cu2_power_check.py without importing it (the import pulls in the OBE model code)."""
    src = (ROOT / "scripts" / "cu2_power_check.py").read_text(encoding="utf-8")
    line = next(l for l in src.splitlines() if l.strip().startswith(f'"{name}":'))
    return eval(line.split(":", 1)[1].split("#")[0].strip().rstrip(","))


def proxy(rem):
    w = np.exp(-0.5 * (np.asarray(PROXY_VX) / 0.65) ** 2); w /= w.sum()
    return np.einsum("v,svy->s", w, rem.reshape(-1, len(PROXY_VX), len(PROXY_Y0))) / len(PROXY_Y0)


def load(stem, n_set):
    parts = vz_paths.glob(OUT, f"{stem}_[0-9][0-9][0-9].npz")
    if not parts:
        return None
    rem = np.concatenate([np.load(p, allow_pickle=True)["remaining"] for p in parts])
    if rem.size < n_set * len(PROXY_VX) * len(PROXY_Y0):
        return None
    return 100 * (1 - proxy(rem))


def main(argv):
    grid, p0, p1, tag = argv[0], int(argv[1]), int(argv[2]), argv[3] if len(argv) > 3 else ""
    d0s, aoms = grid_def(grid)
    field = tag[1:] if tag else "200"
    sfx = f"_{tag}" if tag else ""
    data = {i: load(f"{grid}_loss5_{p0}_{p1}{sfx}_{i}", len(d0s) * len(aoms)) for i in ("spb", "uniform")}
    data = {i: v.reshape(len(d0s), len(aoms)) for i, v in data.items() if v is not None}
    if not data:
        sys.exit("no complete grid yet")
    fig, axes = plt.subplots(1, len(data), figsize=(5.2 * len(data), 3.9), squeeze=False)
    dm = DELTA_MAIN_MHZ.get(field)
    for ax, (init, tab) in zip(axes[0], data.items()):
        lo = np.floor(np.percentile(tab, 20))
        im = ax.imshow(tab, origin="lower", cmap="Blues", vmin=lo, vmax=tab.max(), aspect="auto")
        for (r, c), x in np.ndenumerate(tab):
            k = (r, c) == np.unravel_index(np.argmax(tab), tab.shape)
            ax.text(c, r, f"{x:.2f}", ha="center", va="center", fontsize=8, weight="bold" if k else "normal",
                    color="white" if x > lo + 0.6 * (tab.max() - lo) else INK)
        ax.set_xticks(range(len(aoms)), [f"{a:g}" for a in aoms])
        ax.set_yticks(range(len(d0s)), [f"{d:g}" for d in d0s])
        ax.set_xlabel("$f_{AOM}$ (MHz)" + (f";  $\\delta_1=\\delta_0+f_{{AOM}}-{dm:.2f}$" if dm else ""))
        ax.set_ylabel("$\\delta_0/2\\pi$ (MHz)")
        r, c = np.unravel_index(np.argmax(tab), tab.shape)
        ax.set_title(f"{LABEL[init]}: best {tab.max():.2f} % at ({d0s[r]:g}, {aoms[c]:g})", fontsize=9)
        fig.colorbar(im, ax=ax, label="proxy efficiency (%)", shrink=0.9)
    fig.suptitle(f"CU2 set-point grid at {field} V/cm, {p0} + {p1} mW\n5 % loss per pass, 3 cm beam, 6-point proxy",
                 fontsize=10)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    name = f"CU2_setpoint_grid_{grid}_{p0}_{p1}{sfx}"
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=200)
    print("wrote", FIG / f"{name}.png")
    for init, tab in data.items():
        print(init, "rows delta0", d0s, "cols f_AOM", aoms)
        for d, row in zip(d0s, tab):
            print(f"  {d:5.1f}: " + "  ".join(f"{x:6.2f}" for x in row))


if __name__ == "__main__":
    main(sys.argv[1:])
