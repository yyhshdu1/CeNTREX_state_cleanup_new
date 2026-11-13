"""CU2 correlated laser scan: AOM fixed, laser frequency scanned, so both carriers move together (a 45-degree line
in the (delta0, delta1) plane). On the beam axis, from the `corr` stage of scripts/cu2_power_check.py.

    python scripts/plot_cu2_corr_scan.py 200 32 43 e155 12.5   # f_AOM, powers, field tag, set point delta0 (MHz)
    -> docs/state_cleanup/figures/CU2_corr_scan_e155_f200.{pdf,png}
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
import matplotlib.pyplot as plt  # noqa: E402

INK, MUTED, GRID, SURFACE = "#262624", "#6b6b68", "#e4e3dc", "#ffffff"
COLOR = {"spb": "#2a78d6", "uniform": "#eb6834"}
LABEL = {"spb": "SPB start (4 levels, m_J ≈ 0)", "uniform": "uniform start (20 J=2 levels)"}
mpl.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "axes.labelcolor": INK, "text.color": INK, "savefig.facecolor": SURFACE,
})


DELTA_MAIN_MHZ = {"155": 204.520, "200": 211.829, "217": 214.585}   # F'=4 minus F'=3 main line, from the model log


def load(f_aom, p0, p1, tag, init, field):
    """Merge every saved chunk of the corr scans for this setting (so a scan still running can be plotted too)."""
    pat = f"corr_loss5_f{f_aom:g}_{p0}_{p1}{tag}_{init}_d*_[0-9][0-9][0-9].npz".replace(".", "p")
    pat = pat.replace("pnpz", ".npz")
    pts = {}
    for f in vz_paths.glob(OUT, pat):
        z = np.load(f, allow_pickle=True)
        for row, r in zip(z["rows"], z["remaining"]):
            pts[round(float(eval(str(row), {"np": np})[0][0]), 6)] = float(r)
    if not pts:
        return None
    d0 = np.array(sorted(pts))
    return d0, 100 * (1 - np.array([pts[d] for d in d0])), f_aom - DELTA_MAIN_MHZ[field]


def main(argv):
    f_aom, p0, p1 = float(argv[0]), int(argv[1]), int(argv[2])
    tag = f"_{argv[3]}" if len(argv) > 3 and argv[3].startswith("e") else ""
    d_set = float(argv[4]) if len(argv) > 4 else None
    field = tag[2:] if tag else "200"
    data = {i: load(f_aom, p0, p1, tag, i, field) for i in ("spb", "uniform")}
    data = {i: v for i, v in data.items() if v is not None}
    d1off = next(iter(data.values()))[2]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9), gridspec_kw={"width_ratios": [1.45, 1]})
    for ax, (lo, hi), ylim in zip(axes, [(-62, 142), (-12, 42)], [(0, 101), (90, 100.2)]):
        for init, (d0, eff, _) in data.items():
            m = (d0 >= lo) & (d0 <= hi)
            ax.plot(d0[m], eff[m], color=COLOR[init], lw=1.6, marker="o", ms=2.6, label=LABEL[init], zorder=3)
            k = int(np.argmax(eff))
            if lo <= d0[k] <= hi:
                ax.plot(d0[k], eff[k], marker="*", ms=10, color=COLOR[init], mec=SURFACE, mew=0.8, zorder=4)
        if d_set is not None:
            ax.axvline(d_set, color=MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
        ax.set_xlim(lo, hi); ax.set_ylim(*ylim)
        ax.set_xlabel("F'=3 carrier detuning δ0/2π (MHz)")
        ax.grid(True, color=GRID, lw=0.5, zorder=0); ax.set_axisbelow(True)
        top = ax.secondary_xaxis("top", functions=(lambda x: x + d1off, lambda x: x - d1off))
        top.set_xlabel("F'=4 carrier detuning δ1/2π (MHz)", color=INK)
        top.tick_params(colors=MUTED, labelcolor=INK)
    axes[0].set_ylabel("cleanup efficiency on axis (%)")
    axes[0].set_title(f"(a) laser scan, f_AOM = {f_aom:g} MHz fixed", loc="left")
    axes[1].set_title("(b) near the set point", loc="left")
    if d_set is not None:
        axes[1].text(d_set + 0.6, 90.3, f"set point\nδ0 = {d_set:g} MHz", fontsize=7, color=MUTED, va="bottom")
    axes[0].legend(loc="lower center", fontsize=7.5, frameon=False)
    fig.suptitle(f"CU2, {field} V/cm, {p0} + {p1} mW, 5 % loss/pass, 13 passes, vx = 0, y0 = 0; "
                 f"δ1 = δ0 {d1off:+.2f} MHz; ★ = maximum", fontsize=8.5, color=MUTED, y=0.995)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    stem = FIG / f"CU2_corr_scan{tag}_f{f_aom:g}".replace(".", "p")
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=200)
    print(stem)
    for init, (d0, eff, _) in data.items():
        k = int(np.argmax(eff))
        print(f"{init}: {len(d0)} points, max {eff[k]:.3f} % at delta0 {d0[k]:g} MHz"
              + (f", at set point {np.interp(d_set, d0, eff):.3f} %" if d_set is not None and d0[0] <= d_set <= d0[-1] else ""))


if __name__ == "__main__":
    main(sys.argv[1:])
