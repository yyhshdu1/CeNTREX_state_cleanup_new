"""Slide schematic: the CU2 test 2026-11-08 at a glance -- J=2 populations with CU2 off (fitted), after CU2 (model)
and the cleanup efficiency each detection line measures, against the model.

    python -I scripts/cu2_test_result_schematic.py <CU2_test_processed_results_2026-11-07_08.h5> [power_mW] [ft_from_r]
    -> docs/state_cleanup/figures/slides/CU2_test_2026-11-08_result_schematic[_ft044].{png,pdf}
ft_from_r: start split f_t = 0.44 from S_B/S_A and the model numbers stored in cu2_residual_lines/numbers.json (the
per-group residual solve), instead of the four-signal start fit (41 : 59).

Numbers from det_plates_level_diagrams.residual_cu2test (start fitted to the CU2-off signals with the probe photon
yields, scripts/cu2_test_start_fit.py; CU2 transfer stage on the 24-point acceptance grid at the default v_z) and the
HDF5 attrs clean_eff / clean_eff_lower90.
"""
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import det_plates_level_diagrams as D  # noqa: E402

import matplotlib as mpl  # noqa: E402

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

OUT = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
INK, GRAY, LGRAY = "#111111", "#555555", "#c8c8c8"
COLORS = {(2.5, 3): "#2a78d6", (2.5, 2): "#e8673a", (1.5, 2): "#2e8b57", (1.5, 1): "#7a3db8"}
NAMES = {(2.5, 3): "F₁=5/2, F=3", (2.5, 2): "F₁=5/2, F=2", (1.5, 2): "F₁=3/2, F=2", (1.5, 1): "F₁=3/2, F=1"}
ROLE = {(2.5, 3): "triplet (lens)", (2.5, 2): "singlet (lens)", (1.5, 2): "", (1.5, 1): ""}
YLEV = {(2.5, 3): 3.0, (2.5, 2): 2.0, (1.5, 2): 1.0, (1.5, 1): 0.0}
LINES = [("R(2) F₁′=7/2, F′=4", "R2F4_triplet_nltl2_2068", "F4", (2.5, 3)),
         ("R(2) F₁′=7/2, F′=3", "R2F3_singlet_nltl2_2068", "R7o2F3", (2.5, 2)),
         ("R(2) F₁′=5/2, F′=3", "R2F5o2F3_F1eq3o2F2_nltl2_2068", "F1p5F3", (1.5, 2)),
         ("R(2) F₁′=5/2, F′=2", "R2F5o2F2_F1eq3o2F1_nltl2_2068", "F1p5F2", (1.5, 1))]
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
                     "savefig.facecolor": "white"})


NUMBERS = ROOT / "docs" / "state_cleanup" / "figures" / "cu2_residual_lines" / "numbers.json"
GROUPS = [(2.5, 3), (2.5, 2), (1.5, 2), (1.5, 1)]          # order of numbers.json (cu2_test_start_fit.GROUPS)


def from_numbers():
    """Start split f_t from r = S_B/S_A (the one the per-group residual uses), F1=3/2 empty; model populations after
    CU2, model on/off ratios and CU2-off signal shares as stored by cu2_test_line_combinations.py (no new run)."""
    import json
    d = json.loads(NUMBERS.read_text())
    ft = d["ft"]
    before = dict(zip(GROUPS, (ft, 1 - ft, 0.0, 0.0)))
    after = {g: d["model"][D_NAME[g]] for g in GROUPS}
    ratio = dict(zip([k for _, _, k, _ in LINES], d["model_R"]))
    shares = {k: dict(zip(GROUPS, s)) for (_, _, k, _), s in zip(LINES, d["share_off"])}
    return before, after, ratio, shares


D_NAME = {(2.5, 3): "F1=5/2 F=3", (2.5, 2): "F1=5/2 F=2", (1.5, 2): "F1=3/2 F=2", (1.5, 1): "F1=3/2 F=1"}


def main(h5path, power_mW=24.0, ft_from_r=False):
    _, _, t = D.residual_cu2test(h5path, power_mW)
    if ft_from_r:
        before, after, ratio, shares = from_numbers()
        t = dict(t, ratio=ratio, shares=shares)
    else:
        before = {g: sum(v for l, v in t["rho0"].items() if (l[0], l[1]) == g) for g in YLEV}
        after = {g: sum(v for l, v in t["after"].items() if (l[0], l[1]) == g) for g in YLEV}
    with h5py.File(h5path, "r") as f:
        meas = {k: (100 * float(f["cu2_onoff"][run].attrs["clean_eff"]),
                    100 * float(f["cu2_onoff"][run].attrs["clean_eff_lower90"])) for _, run, k, _ in LINES}
    model = {k: 100 * (1 - t["ratio"][k]) for _, _, k, _ in LINES}

    fig = plt.figure(figsize=(12.4, 5.2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 12.4); ax.set_ylim(-1.55, 4.75); ax.axis("off")

    def panel(x0, title, sub):
        ax.text(x0, 4.45, title, fontsize=11, weight="bold", color=INK, va="center")
        ax.text(x0, 4.08, sub, fontsize=8.3, color=GRAY, va="center")

    def levels(x0, w, vals, scale, fmt):
        for g, y in YLEV.items():
            ax.plot([x0, x0 + w], [y, y], color=INK, lw=1.2, solid_capstyle="butt")
            v = vals[g]
            L = min(w, w * v / scale)
            if L > 0:
                ax.add_patch(FancyBboxPatch((x0, y + 0.06), L, 0.36, boxstyle="square,pad=0", fc=COLORS[g], ec="none",
                                            alpha=0.85))
            ax.text(x0 + w + 0.1, y + 0.22, fmt(v), fontsize=10, color=COLORS[g], va="center", weight="bold")

    # (1) CU2 off
    x1, w1 = 1.55, 1.9
    panel(0.15, "① CU2 off", "F = 3 share from S$_B$/S$_A$ of the CU2-off signals, F₁ = 3/2 empty" if ft_from_r
          else "start populations fitted to the four CU2-off signals")
    for g, y in YLEV.items():
        ax.text(1.42, y + 0.22, NAMES[g], fontsize=9, ha="right", va="center", color=INK)
        if ROLE[g]:
            ax.text(1.42, y - 0.17, ROLE[g], fontsize=7.5, ha="right", va="center", color=GRAY)
    levels(x1, w1, before, 0.6, lambda v: f"{100 * v:.0f} %")
    ax.text(0.15, -0.6, "X ¹Σ⁺, J = 2 (bar length: share of all molecules)", fontsize=7.8, color=GRAY)

    # CU2 arrow
    ax.add_patch(FancyArrowPatch((4.15, 1.6), (5.15, 1.6), arrowstyle="simple", mutation_scale=34, color="#f2c14e"))
    ax.text(4.65, 2.55, "CU2", fontsize=11, weight="bold", ha="center", color=INK)
    ax.text(4.65, 0.75, "R(2) F₁′=7/2,\nF′ = 3 + 4\n217 V/cm\nδ₀ 13 MHz, f$_\\mathrm{AOM}$ 209\n20 + 30 mW",
            fontsize=7.6, ha="center", va="top", color=GRAY, linespacing=1.25)

    # (2) after CU2
    x2, w2 = 5.55, 1.9
    panel(5.45, "② after CU2 (model)", f"v$_z$ = {t['vz']:g} m/s, mean over the acceptance")
    levels(x2, w2, after, 0.006, lambda v: f"{100 * v:.2f} %")
    tot = sum(after.values())
    ax.text(5.45, -0.6, f"J = 2 left: {100 * tot:.2f} % (bar scale ×100 of ①);\n"
                        f"{100 * (after[(1.5, 2)] + after[(1.5, 1)]) / tot:.0f} % of it is F₁ = 3/2, pumped in by CU2",
            fontsize=7.8, color=GRAY, va="top", linespacing=1.3)

    # (3) detection lines
    x3 = 8.65
    panel(x3, "③ measured vs model", "cleanup efficiency seen by each line, DET plates off, "
                                    f"{power_mW:g} mW")
    ax.text(x3 + 2.05, 3.62, "measured", fontsize=8.5, ha="center", color=INK, weight="bold")
    ax.text(x3 + 3.15, 3.62, "model", fontsize=8.5, ha="center", color=INK, weight="bold")
    for (name, _, k, g), y in zip(LINES, (3.0, 2.0, 1.0, 0.0)):
        m, lo = meas[k]
        ok = lo <= model[k] <= 100
        ax.text(x3, y + 0.22, name, fontsize=8.8, color=COLORS[g], va="center", weight="bold")
        ax.text(x3 + 2.05, y + 0.32, f"{m:.1f} %", fontsize=10, ha="center", va="center", color=INK)
        ax.text(x3 + 2.05, y + 0.0, f"> {lo:.1f} % (90 % CL)", fontsize=7.2, ha="center", va="center", color=GRAY)
        ax.text(x3 + 3.15, y + 0.22, f"{model[k]:.1f} %", fontsize=10, ha="center", va="center", color=INK)
        ax.text(x3 + 3.6, y + 0.22, "✓" if ok else "✗", fontsize=12, va="center", color="#2e8b57" if ok else "#c0392b")
    sh = t["shares"]
    ax.text(x3, -0.6, "with CU2 off the two F₁ = 3/2 lines see mostly F₁ = 5/2 molecules\n"
                      f"(R(2) 5/2, F′=3: {100 * sh['F1p5F3'][(2.5, 3)]:.0f} % F=3; "
                      f"R(2) 5/2, F′=2: {100 * sh['F1p5F2'][(2.5, 2)]:.0f} % singlet); ✓ = model inside the 90 % CL",
            fontsize=7.8, color=GRAY, va="top", linespacing=1.3)
    for xs in (8.45,):
        ax.plot([xs, xs], [-0.4, 4.6], color=LGRAY, lw=0.8)
    OUT.mkdir(parents=True, exist_ok=True)
    name = "CU2_test_2026-11-08_result_schematic" + ("_ft044" if ft_from_r else "")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=220, bbox_inches="tight", pad_inches=0.08)
    print("wrote", OUT / f"{name}.png")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "ft_from_r"]
    main(args[0], float(args[1]) if len(args) > 1 else 24.0, ft_from_r="ft_from_r" in sys.argv)
