"""Overlay the CU2 resonances (217 V/cm model) on the measured CU2 frequency-scan figure (2026-11-07, 20 + 30 mW).

    python scripts/annotate_cu2_scan_transitions.py <measured.png>
    -> docs/state_cleanup/figures/CU2_scan_2026-11-07_transitions.png
    python scripts/annotate_cu2_scan_transitions.py <measured.png> --sim <corrdet result.npz> [...]
    -> docs/state_cleanup/figures/CU2_scan_2026-11-07_sim.png   (simulated F=3 left on the lower panel)

The measured figure has no data file here, so the overlay is drawn in its pixel frame: the x axis is calibrated from
its grid lines (x = -40, -20, 0, 20, 40 MHz IR at columns 444.5, 782.5, 1121.5, 1459.5, 1797.5 of the 1950 px wide
image). Resonance ranges: X J=2 F1=5/2 F=3 (all mF) -> B J'=3 F1'=7/2 F'=3, 4 (both Lambda-doublet parities), X + Z
light, E = 217 V/cm, carrier only (no EOM sidebands, no Doppler); x_IR = (delta0 - 14.5 MHz) / 4, where the F'=4
beam is the laser + 209 MHz (UV). Computed in this session from generate_reduced_X/B_hamiltonian.
"""
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

X0_PX, PX_PER_MHZ = 1121.5, (1797.5 - 444.5) / 80.0
PANEL_BOTTOM = (723, 1378)            # rows of the lower panel frame
PANEL_TOP = (72, 669)
C3, C4, INK, MUTED = "#7a3db8", "#1baf7a", "#262624", "#6b6b68"     # F'=3 beam, F'=4 beam
# (beam, excited manifold, lo, hi, weighted centre) in MHz IR on the measured axis
LINES = [("F'=3 beam", "F'=3, P'=−1", 0.2, 7.3, 3.7), ("F'=3 beam", "F'=3, P'=+1", 11.8, 19.6, 14.4),
         ("F'=3 beam", "F'=4, P'=−1", 41.5, 52.4, 47.4), ("F'=3 beam", "F'=4, P'=+1", 56.0, 65.1, 60.9),
         ("F'=4 beam", "F'=3, P'=−1", -52.1, -45.0, -48.6), ("F'=4 beam", "F'=3, P'=+1", -40.4, -32.6, -37.8),
         ("F'=4 beam", "F'=4, P'=−1", -10.7, 0.1, -4.9), ("F'=4 beam", "F'=4, P'=+1", 3.8, 12.8, 8.6)]


def px(x):
    return X0_PX + PX_PER_MHZ * np.asarray(x)


Y0_PX, PX_PER_UNIT = 1337.0, 512.5        # lower panel: remaining fraction 0 and 1 rows from its grid lines


def overlay_sim(path, files):
    img = np.asarray(Image.open(path).convert("RGB"))
    H, W, _ = img.shape
    extra = 120
    fig = plt.figure(figsize=(W / 150, (H + extra) / 150), dpi=150)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, extent=(0, W, H, 0))
    ax.set_xlim(0, W); ax.set_ylim(H + extra, 0); ax.axis("off")
    styles = [("#262624", "-"), ("#7a3db8", "-"), ("#1baf7a", "--")]
    notes = []
    for (col, ls), f in zip(styles, files):
        z = np.load(f, allow_pickle=True)
        x = (z["delta0_MHz"] - 14.5) / 4
        vx = np.asarray(z["vx"]); f3 = np.asarray(z["F3"])
        w = np.exp(-0.5 * (vx / 0.65) ** 2); w /= w.sum()
        y = w @ f3
        lab = "on axis" if len(vx) == 1 else "vx " + ", ".join(f"{v:+g}" for v in vx) + " m/s (Gaussian weights)"
        if "y0" in z.files:
            lab += ", y0 " + ", ".join(f"{1e3 * v:g}" for v in z["y0"]) + " mm"
        ax.plot(px(x), Y0_PX - PX_PER_UNIT * np.clip(y, -0.08, 1.2), color=col, lw=2.0, ls=ls, zorder=5)
        p0, p1 = (int(v) for v in z["powers_mW"])
        notes.append((col, f"simulation, F=3 left ({lab}), {p0} + {p1} mW, triplet start"))
    for k, (col, t) in enumerate(notes):
        ax.plot([160, 230], [H + 35 + 32 * k] * 2, color=col, lw=2.0)
        ax.text(245, H + 35 + 32 * k, t, fontsize=9.5, va="center", color=col)
    ax.text(W - 20, H + 35, "217 V/cm, (14.5, 209) MHz set point; x = (δ0 − 14.5 MHz)/4; same-beam couplings only\n"
            "(the outer dips, one beam on the other beam's line, are not in the model)", fontsize=8.5,
            color="#6b6b68", va="top", ha="right")
    out = ROOT / "docs" / "state_cleanup" / "figures" / "CU2_scan_2026-11-07_sim.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)


def main(path):
    img = np.asarray(Image.open(path).convert("RGB"))
    H, W, _ = img.shape
    extra = 330
    fig = plt.figure(figsize=(W / 150, (H + extra) / 150), dpi=150)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, extent=(0, W, H, 0))
    ax.set_xlim(0, W); ax.set_ylim(H + extra, 0); ax.axis("off")
    xmin, xmax = (140 - X0_PX) / PX_PER_MHZ, (1915 - X0_PX) / PX_PER_MHZ
    for beam, man, lo, hi, c in LINES:
        col = C3 if beam.startswith("F'=3") else C4
        a = 0.16 if "−1" in man else 0.08
        lo_c, hi_c = max(lo, xmin), min(hi, xmax)
        if lo_c < hi_c:
            for r0, r1 in (PANEL_TOP, PANEL_BOTTOM):
                ax.fill_between([px(lo_c), px(hi_c)], r0 + 2, r1 - 2, color=col, alpha=a, lw=0, zorder=2)
            if xmin < c < xmax:
                ax.plot([px(c)] * 2, [PANEL_BOTTOM[0] + 2, PANEL_BOTTOM[1] - 2], color=col, lw=1.2, ls=(0, (5, 3)))
    # ruler below the figure: one row per beam
    rows = {"F'=3 beam": H + 80, "F'=4 beam": H + 185}
    for beam, y in rows.items():
        col = C3 if beam.startswith("F'=3") else C4
        ax.text(px(xmin), y + 16, beam + (" = laser" if beam.startswith("F'=3") else " = laser + 209 MHz (UV AOM)"),
                color=col, fontsize=9.5, va="top", ha="left", weight="bold")
        ax.plot([px(xmin), px(xmax)], [y, y], color=MUTED, lw=0.5)
    for beam, man, lo, hi, c in LINES:
        col = C3 if beam.startswith("F'=3") else C4
        y = rows[beam]
        lo_c, hi_c = max(lo, xmin), min(hi, xmax)
        if lo_c >= hi_c:
            ax.text(px(xmax), y + 16, f"(→ {man} at +{lo:.0f}…+{hi:.0f} MHz, off scale)", color=col, fontsize=8.5,
                    ha="right", va="top")
            continue
        solid = "−1" in man
        ax.plot([px(lo_c), px(hi_c)], [y, y], color=col, lw=9 if solid else 5, alpha=1.0 if solid else 0.55,
                solid_capstyle="butt")
        edge = hi_c >= xmax - 1e-9
        ax.text(px(hi_c) if edge else px(0.5 * (lo_c + hi_c)), y - 14, f"→ {man}", color=INK, fontsize=8.5,
                ha="right" if edge else "center", va="bottom")
    ax.text(20, H + 260,
            "Model resonances at 217 V/cm (carrier only; EOM sidebands add ±4.4 and ±8.8 MHz IR, no Doppler): "
            "X J=2 F1=5/2 F=3, all mF → B J'=3 F1'=7/2 F', X + Z light.\n"
            "Shaded: range of the components, dark = main parity P'=−1, light = Λ-doublet partner P'=+1 "
            "(Stark-mixed); dashed: strength-weighted centre.  x_IR = (δ0 − 14.5 MHz)/4.",
            fontsize=8.5, color=MUTED, va="top")
    out = ROOT / "docs" / "state_cleanup" / "figures" / "CU2_scan_2026-11-07_transitions.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    if "--sim" in sys.argv:
        overlay_sim(sys.argv[1], sys.argv[sys.argv.index("--sim") + 1:])
    else:
        main(sys.argv[1])
