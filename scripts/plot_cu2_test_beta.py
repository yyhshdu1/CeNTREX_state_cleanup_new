"""CU2 test 2026-11-07 (B1 scan, triplet start) against the on-axis model at several EOM modulation depths beta.

    python -I scripts/plot_cu2_test_beta.py <CU2_test_processed_results_2026-11-07_08.h5>
    -> docs/state_cleanup/figures/CU2_test_2026-11-07_beta.{png,pdf}

Model: corrdet stage of scripts/cu2_power_check.py, 20 + 30 mW, f_AOM 209 MHz, 217 V/cm, 5 % loss per pass, vx = 0,
y0 = 0, triplet start (F1=5/2 F=3 mF=0,±1), F=3 left. beta = 1.5 is the default (no beta token). The R(2) F'=4 line
detects F=3 only, so the measured ratio is compared with F=3 left directly.
"""
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import vz_paths  # noqa: E402
OUT = ROOT / "tmp" / "cu2_power"
FIG = ROOT / "docs" / "state_cleanup" / "figures"
import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK, MUTED, GRID = "#262624", "#6b6b68", "#e4e3dc"
mpl.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED,
                     "ytick.color": MUTED, "xtick.labelcolor": INK, "ytick.labelcolor": INK, "axes.labelcolor": INK,
                     "text.color": INK})
STEM = "corrdet_loss5_f209_20_30_e217_trip_vx0"
RUNS = [(0.0, f"{STEM}_d-66_94_4_beta0", "#eb6834"), (0.75, f"{STEM}_d-66_94_4_beta0p75", "#c9a227"),
        (1.0, f"{STEM}_d-66_94_4_beta1", "#1baf7a"), (1.5, f"{STEM}_d-130_170_4", "#2a78d6")]
WINDOW = (-20.0, 20.0)          # MHz IR, the range every run covers; the other beam's lines lie outside


def curve(stem):
    z = np.load(vz_paths.find(OUT / f"{stem}_result.npz"), allow_pickle=True)
    return (z["delta0_MHz"] - 14.5) / 4, np.asarray(z["F3"]).reshape(-1)


def crossings(x, y, level=0.5):
    """Outermost level crossings around the deepest point (linear interpolation)."""
    i0 = int(np.argmin(y))
    lo = next((x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
               for i in range(i0 - 1, -1, -1) if y[i] >= level > y[i + 1]), np.nan)
    hi = next((x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
               for i in range(i0, len(x) - 1) if y[i] < level <= y[i + 1]), np.nan)
    return lo, hi


def main(h5path):
    with h5py.File(h5path, "r") as f:
        b1 = f["cu2_frequency_scans/B1_pm50MHz_on_R2F4"]
        xd, Rd, ed = (b1[k][()] for k in ("dIR_MHz", "R_template", "sem_bootstrap_R_template"))
    o = np.argsort(xd); xd, Rd, ed = xd[o], Rd[o], ed[o]
    m = (xd >= WINDOW[0]) & (xd <= WINDOW[1])
    lo_d, hi_d = crossings(xd, Rd)
    print(f"data: R=0.5 at {lo_d:+.1f} / {hi_d:+.1f} MHz IR, width {hi_d - lo_d:.1f}")
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(9.5, 6.6), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    ax.errorbar(xd, Rd, ed, fmt="o", ms=3.5, color=INK, ecolor=MUTED, elinewidth=0.8, capsize=0, zorder=5,
                label=f"measured 11-07 (R = 0.5 at {lo_d:+.1f} / {hi_d:+.1f})")
    rows = []
    for beta, stem, col in RUNS:
        x, y = curve(stem)
        lo, hi = crossings(x, y)
        yi = np.interp(xd[m], x, y)
        rms = float(np.sqrt(np.mean((Rd[m] - yi) ** 2)))
        rows.append((beta, lo, hi, rms))
        ax.plot(x, y, color=col, lw=1.6, zorder=3,
                label=f"model β = {beta:g} (R = 0.5 at {lo:+.1f} / {hi:+.1f}; rms {rms:.3f})")
        axr.plot(xd[m], Rd[m] - yi, "o-", ms=2.5, lw=0.9, color=col)
        print(f"beta {beta:4g}: R=0.5 at {lo:+.1f} / {hi:+.1f}, width {hi - lo:.1f}, rms(data-model) {rms:.3f} "
              f"over {WINDOW}")
    for a in (ax, axr):
        a.axvline(-0.36, color=MUTED, lw=0.7, ls=(0, (2, 2)))
        a.spines[["top", "right"]].set_visible(False)
        a.axhline(0, color=GRID, lw=0.8, zorder=0)
    ax.axhline(1, color=GRID, lw=0.8, zorder=0)
    ax.text(-0.06, 1.32, "11-08 final", fontsize=7, color=MUTED, ha="right")
    ax.set_xlim(*WINDOW); ax.set_ylim(-0.08, 1.45)
    ax.set_ylabel("F1=5/2 F=3 left\nS(CU2 on) / S(CU2 off)")
    ax.legend(frameon=False, loc="upper left", fontsize=7.8)
    ax.set_title("CU2 triplet scan vs EOM modulation depth β (sidebands ±1, ±2 at ±4.4, ±8.8 MHz IR); "
                 "on axis, 20 + 30 mW, 209 MHz, 217 V/cm", loc="left", fontsize=9)
    axr.set_ylabel("measured − model")
    axr.set_xlabel("CU2 laser − 275858.76318 GHz (MHz IR, Cs-corrected)  =  (δ0 − 14.5 MHz) / 4")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"CU2_test_2026-11-07_beta.{ext}", dpi=200)
    print("wrote", FIG / "CU2_test_2026-11-07_beta.png")


if __name__ == "__main__":
    main(sys.argv[1])
