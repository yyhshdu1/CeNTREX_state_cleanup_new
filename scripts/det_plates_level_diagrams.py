"""Level diagrams of the three CU2-test detection lines with the DET plates on (E along k) or off (E = 0).

    python scripts/det_plates_level_diagrams.py [power_mW]       # plates on, 15 (default) or 3 mW
    python scripts/det_plates_level_diagrams.py 15 off           # plates off, zero field (15 mW only)
    python scripts/det_plates_level_diagrams.py 15 off Q3o2F1    # only the named lines (comma separated);
    python scripts/det_plates_level_diagrams.py 15 off R7o2F3,F1p5F3,F1p5F2 <CU2-test h5>
                                       # CU2 test 2026-11-08 overlays (217 V/cm, (13, 209) MHz, 20 + 30 mW):
                                       # lens start fitted to the CU2-off signals -> *_cu2test217.{png,pdf}
    (append "notext" for the same diagram without the text block underneath -> *_notext.{png,pdf})
                                                                 # Q(2) F1'=3/2 F'=1 and P(2) F1'=1/2 F'=0, 1 have
                                                                 # plates-off yields only

Plates off: E = 0, the quantization axis is the earth field (vertical, = the CU2 field axis, perpendicular to k),
so the sigma+ light about k is 1/2 pi + 1/4 sigma+ + 1/4 sigma- about the axis (arrows Delta mF = 0, +-1). The field
falls to zero adiabatically after CU2, so every level keeps its label and mF about the vertical. Photons per level
are those of the plates-off CU2-test table (cu2_detection_test.level_photons, POL "k": probe computed with the axis
along k, where the light is pure sigma+, and each vertical-axis level rotated onto k with |d^F(pi/2)|^2), vx = 0.

The diagrams are utils.plotting.plot_transition_level_diagram at 154 V/cm (levels segmented by zero-field
parent character: X by (F1, F), B by the two parities), with σ+ light about the field and these overlays on J=2:
    yellow       populated with CU2 off: triplet-only lens start (RC on P(2) F'=2), F1=5/2 F=3, mF=0, ±1
    blue number  photons per molecule starting in that level, probe on the strongest line of that level,
                 vx = 0, flat top 1.5 x 1.5 cm at the given power (field-on probe OBE)
    blue outline the level gives >= 0.1 photon, i.e. the line sees it
    black number population left in that level after CU2 (units 1e-4), triplet start, acceptance mean
The E field turns 90 deg from CU2 (vertical) to DET (along k) adiabatically and never passes zero, so each
CU2 level arrives in the DET level with the same label and m. Probe OBE results (vx = 0, per starting level,
detuning scan) are cached in tmp/detection/plates_on/ (scratch script spec_obe.py of 2026-11-04).
"""
import dataclasses
import json
import sys
import textwrap
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys  # noqa: E402
sys.path.insert(0, str(ROOT))
from scripts import vz_paths  # noqa: E402
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT))     # this repo's centrex_tlf (the conda env holds an older non-editable copy)
import cu2_detection_test as T  # noqa: E402

CACHE = ROOT / "tmp" / "detection" / "plates_on"
FIG_DIR = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
E_DET = 154.0
LINES = {"F4": ("R2_F1_7o2_F4", "R(2) $F_1'=7/2$, $F'=4$", "R2F4"),
         "Q5o2F3": ("Q2_F1_5o2_F3", "Q(2) $F_1'=5/2$, $F'=3$", "Q2F5o2F3"),
         "Q3o2F2": ("Q2_F1_3o2_F2", "Q(2) $F_1'=3/2$, $F'=2$", "Q2F3o2F2"),
         "Q3o2F1": ("Q2_F1_3o2_F1", "Q(2) $F_1'=3/2$, $F'=1$", "Q2F3o2F1"),
         "P1o2F1": ("P2_F1_1o2_F1", "P(2) $F_1'=1/2$, $F'=1$", "P2F1o2F1"),
         "P1o2F0": ("P2_F1_1o2_F0", "P(2) $F_1'=1/2$, $F'=0$", "P2F1o2F0"),
         "P3o2F1": ("P2_F1_3o2_F1", "P(2) $F_1'=3/2$, $F'=1$", "P2F3o2F1"),
         "P3o2F2": ("P2_F1_3o2_F2", "P(2) $F_1'=3/2$, $F'=2$", "P2F3o2F2"),
         "R7o2F3": ("R2_F1_7o2_F3", "R(2) $F_1'=7/2$, $F'=3$", "R2F7o2F3"),
         "F1p5F3": ("R2_F1_5o2_F3", "R(2) $F_1'=5/2$, $F'=3$", "R2F5o2F3"),
         "F1p5F2": ("R2_F1_5o2_F2", "R(2) $F_1'=5/2$, $F'=2$", "R2F5o2F2")}
DEFAULT_LINES = ("F4", "Q5o2F3", "Q3o2F2")
MIN_SEP = 0.02


def residual_triplet():
    """CU2-on population per J=2 level (F1, F, mF), triplet start, mean over the detection acceptance."""
    m = T.Model((12.5, 200.0), 32, 43, "_e155")
    trip = [s for s in m.fin if s[1] == 3]
    per = {lv: sum(m.mean(m.fin[s][..., k]) for s in trip) / len(trip) for k, lv in enumerate(m.levels)}
    return per, trip


ONOFF_LINES = {"R2F4_triplet_nltl2_2068": "F4", "R2F3_singlet_nltl2_2068": "R7o2F3",
               "R2F5o2F3_F1eq3o2F2_nltl2_2068": "F1p5F3", "R2F5o2F2_F1eq3o2F1_nltl2_2068": "F1p5F2"}
START_GROUPS = {(2.5, 3): (0, -1, 1), (2.5, 2): (0,), (1.5, 2): (-2, -1, 0, 1, 2), (1.5, 1): (-1, 0, 1)}


def residual_cu2test(h5path, power_mW=24):
    """CU2 test 2026-11-08: start populations fitted (non-negative least squares) to the four CU2-off signals with the
    plates-off probe yields (power_mW, sigma+ about k, Gaussian vx average, sigma 0.65 m/s), the lens levels at their
    m_J ~ 0 sublevels and F1=3/2 uniform over mF; then the CU2 transfer stage (217 V/cm, (13, 209) MHz, 20 + 30 mW,
    5 % loss per pass, 24-point acceptance mean) -> population per J=2 level after CU2."""
    import cu2_test_start_fit as SF                 # the start-population fit (run it alone to see every step)
    import plot_cu2_test_vs_sim as Pv
    assert list(SF.LINES) == list(ONOFF_LINES) and SF.GROUPS == START_GROUPS
    S_off = np.array([v["S_off"] for v in SF.measured(h5path).values()])
    z = np.load(Pv.SIM_TRANSFER_ACC, allow_pickle=True)
    lab = [(float(a), int(b), int(c)) for a, b, c in zip(z["F1"], z["F"], z["mF"])]
    st = [(float(a), int(b), int(c)) for a, b, c in zip(z["start_F1"], z["start_F"], z["start_mF"])]
    Y = SF.level_yields(power_mW, lab)
    A = SF.yield_matrix(Y, lab)
    _, wN, S_pred = SF.fit(A, S_off)
    N = wN * S_pred[0] / (A[0] @ wN)
    w = dict(zip(START_GROUPS, wN))
    fin = np.einsum("vy,svyl->sl", Pv.acceptance_weights(z["vx"], z["y0"], svx=0.65), z["final"])
    mirror = {l: lab.index((l[0], l[1], -l[2])) for l in lab}

    def evolve(l):
        f_ = fin[st.index((l[0], l[1], abs(l[2])))]
        return f_ if l[2] >= 0 else np.array([f_[mirror[q]] for q in lab])
    rho0 = np.zeros(len(lab))
    for g_, ms in START_GROUPS.items():
        for m in ms:
            rho0[lab.index((g_[0], g_[1], m))] = w[g_] / len(ms)
    after = sum(r * evolve(l) for l, r in zip(lab, rho0) if r > 0)
    start = [l for l, r in zip(lab, rho0) if r > 0]
    ratio = {line: float(Y[line] @ after / (Y[line] @ rho0)) for line in ONOFF_LINES.values()}
    shares = {line: {g_: float(sum(Y[line][lab.index((g_[0], g_[1], m))] * rho0[lab.index((g_[0], g_[1], m))]
                                   for m in ms) / (Y[line] @ rho0)) for g_, ms in START_GROUPS.items()}
              for line in ONOFF_LINES.values()}
    S_fit = dict(zip(ONOFF_LINES.values(), A @ N))
    S_meas = dict(zip(ONOFF_LINES.values(), S_off))
    vz = 184 if "vz184" in str(Pv.SIM_TRANSFER_ACC) else 156
    txt = dict(
        start=("Yellow: populated with CU2 off, fitted to the four measured CU2-off signals with these probe yields "
               "(lens levels at m_J ~ 0, F1=3/2 uniform in mF): " +
               ", ".join(f"F1={int(2 * a)}/2 F={b} {100 * v:.0f} %" for (a, b), v in w.items()) +
               ". With the saturated probe, the lens levels alone account for the CU2-off signals on the F1=3/2 lines."),
        resid=(f"Black number right of each J=2 level: population left in that level after CU2, units 1e-4 of the "
               f"starting population (CU2 test: 217 V/cm, (13, 209) MHz, 20 + 30 mW, v_z = {vz} m/s, mean over |vx| <= "
               f"1 m/s and |y0| <= 8.5 mm). J=2 total {100 * after.sum():.2f} %. Model on/off ratio of this line, "
               "same probe: {ratio}."),
        ratio=ratio, vz=vz, shares=shares, weights=w, S_fit=S_fit, S_meas=S_meas, power_mW=power_mW,
        Y=Y, rho0=dict(zip(lab, rho0)), after=dict(zip(lab, after)), labels=lab)
    return dict(zip(lab, after)), start, txt


def line_ratio(line, power_mW, txt):
    """Model CU2 on/off ratio of any line in the CU2-test overlay, as text (vx-averaged probe yields)."""
    if line in txt["ratio"]:
        return f"{txt['ratio'][line]:.4f}"
    import cu2_test_start_fit as SF
    lab = txt["labels"]
    Y = SF.level_yields(power_mW, lab)[line] if line in dict(SF.LINES.values()) else None
    if Y is None:
        T.POL = "k"
        w = np.exp(-0.5 * (SF.VXS / SF.SVX) ** 2)
        Y = (w / w.sum()) @ T.level_photons(line, lab, SF.VXS, power=1e-3 * power_mW)
    on = sum(Y[i] * txt["after"][l] for i, l in enumerate(lab))
    off = sum(Y[i] * txt["rho0"][l] for i, l in enumerate(lab))
    if off < 0.02 * on or off < 1e-4:
        return f"n/a (no CU2-off signal; CU2 on: {on:.4f} photons per J=2 molecule)"
    return f"{on / off:.4f}"


def photons(line, power_mW, plates="on"):
    if plates == "off":
        m = T.Model((12.5, 200.0), 32, 43, "_e155")
        T.POL = "k"
        return dict(zip(m.levels, T.level_photons(line, m.levels, [0.0], power=1e-3 * power_mW)[0]))
    sfx = "" if power_mW == 15 else f"_{power_mW:g}mW"
    d = json.loads(vz_paths.find(CACHE / f"specobe_{line}_154{sfx}.json", CACHE.parent).read_text())
    ph = np.asarray(d["photons"])
    return {tuple(lab): float(r.max()) for lab, r in zip(d["labels"], ph)}


def plot(power_mW=15, plates="on", lines=DEFAULT_LINES, h5path=None, notext=False):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, Rectangle
    from matplotlib.text import Annotation
    from centrex_tlf import transitions
    from centrex_tlf.utils import plotting as lp

    YELLOW, INK, BLUE = "#f6e3a1", "#111111", "#1537a6"
    width, pad = 0.30, 0.03
    res, trip, txt = residual_cu2test(h5path, power_mW) if h5path else (*residual_triplet(), None)
    on = plates == "on"
    E = E_DET if on else 0.0
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    def fmt(v):
        return f"{v:.0f}" if v >= 10 else f"{v:.1f}" if v >= 1 else f"{v:.2f}"

    for line in lines:
        tname, tlabel, short = LINES[line]
        trans = getattr(transitions, tname)
        ph = photons(line, power_mW, plates)
        structure = lp.calculate_transition_level_structure(trans, E=E, B=0.0)
        if not on:      # E = 0: parity is good; the Lambda-doublet partner with the ground parity does not couple
            structure = dataclasses.replace(
                structure, excited=tuple(e for e in structure.excited if e.P == trans.P_excited))
        diag = lp.plot_transition_level_diagram(
            structure=structure, figsize=(25.0, 12.5), min_level_separation=MIN_SEP, level_width=width,
            title=(f"Detection, DET plates on  |  {tlabel}  |  E = {E_DET:g} V/cm along k  |  σ+ {power_mW:g} mW" if on
                   else f"Detection, DET plates off  |  {tlabel}  |  E = 0, earth field vertical (⊥ k)  |  "
                        f"σ+ about k, {power_mW:g} mW"))
        ax, fig = diag.ax, diag.fig
        ground, excited = list(structure.ground), list(structure.excited)
        for t in list(ax.texts):             # drop the generic arrow of the diagram
            if (isinstance(t, Annotation) and t.get_text() == "") or (
                    "\n" in t.get_text() and t.get_transform() == ax.transData):
                t.remove()

        ys = lp._band_y(ground, *lp._Y_GROUND_BAND, MIN_SEP)
        ys_e = lp._band_y(excited, *lp._Y_EXCITED_BAND, MIN_SEP)
        missing = [(lv.F1, lv.F, lv.mF) for lv in ground if (lv.F1, lv.F, lv.mF) not in res]
        assert not missing, missing
        for lv, y in zip(ground, ys):
            key = (lv.F1, lv.F, lv.mF)
            x0 = lv.mF - width / 2 - pad
            if key in trip:
                ax.add_patch(Rectangle((x0, y - 0.009), width + 2 * pad, 0.018, facecolor=YELLOW,
                                       edgecolor="none", zorder=0))
            n = ph.get(key, 0.0)
            if n >= 0.1:
                ax.add_patch(Rectangle((x0, y - 0.009), width + 2 * pad, 0.018, facecolor="none",
                                       edgecolor=BLUE, lw=1.6, zorder=3))
            ax.text(x0 - 0.03, y, fmt(n), ha="right", va="center", fontsize=9, color=BLUE,
                    fontweight="bold" if n >= 0.1 else "normal")
            ax.text(lv.mF + width / 2 + pad + 0.03, y, fmt(1e4 * res[key]), ha="left", va="center",
                    fontsize=9, color=INK)
        # arrows from the CU2-off levels to the excited levels they drive: plates on sigma+ (mF + 1); plates off
        # the light is 1/2 pi + 1/4 sigma+ + 1/4 sigma- about the vertical axis (line width ~ weight)
        # if the line does not see the CU2-off levels at all (Q(2) F'=1: F=3 -> F'=1 is Delta F = 2), the arrows
        # start from the line's main ground level instead
        trip_seen = any(ph.get(k, 0.0) >= 0.1 for k in trip)
        main = T.LINES[line][1]
        if txt is not None:                  # CU2 test: arrows from the start levels this line sees
            trip_seen = True
        seen_start = {k for k in trip if ph.get(k, 0.0) >= 0.1}
        if txt is not None and not seen_start:      # the line sees none of the CU2-off levels
            trip_seen = False
        arrow_from = (seen_start if txt is not None and seen_start else set(trip) if trip_seen else
                      {(lv.F1, lv.F, lv.mF) for lv in ground if (lv.F1, lv.F) == (main["F1"], main["F"])})
        arrow_text = ("from the levels populated with CU2 off that this line sees" if txt is not None and seen_start else
                      "from the levels populated with CU2 off" if trip_seen else
                      f"from the main level F1={main['F1']:g} F={main['F']} (the line does not see the levels "
                      f"populated with CU2 off)")
        steps = {1: 1.4} if on else {0: 1.8, 1: 0.9, -1: 0.9}
        for lv, y in zip(ground, ys):
            if (lv.F1, lv.F, lv.mF) not in arrow_from:
                continue
            for e, ye in zip(excited, ys_e):
                if e.mF - lv.mF in steps:
                    ax.add_patch(FancyArrowPatch((lv.mF, y + 0.012), (e.mF, ye - 0.012), arrowstyle="-|>",
                                                 mutation_scale=13, lw=steps[e.mF - lv.mF], color=BLUE, alpha=0.75,
                                                 shrinkA=0, shrinkB=0, zorder=2.5))

        f15 = sum(v for k, v in res.items() if k[0] == 1.5)
        seen15 = [k for k in res if k[0] == 1.5 and ph.get(k, 0.0) >= 0.1]
        key_on = [
            ("σ+ about the field (E along the probe k): ground mF -> excited mF' = mF + 1. Arrows: " + arrow_text +
             ". The E field splits the excited mF' by several MHz (resolved lines); the 20 J=2 "
             f"levels lie within {np.ptp([lv.energy_MHz for lv in ground]):.1f} MHz, so levels of equal mF are not "
             "resolved, and the field mixes F1=5/2 and F1=3/2 (segmented bars)."),
            ("Yellow: populated with CU2 off, triplet-only lens start (RC on P(2) F'=2): F1=5/2 F=3, mF=0, ±1, 1/3 each."),
            (f"Blue number left of each J=2 level: photons per molecule starting in that level, probe on the strongest "
             f"line of that level, vx = 0, {power_mW:g} mW, flat top 1.5 x 1.5 cm, 81.5 us (field-on probe OBE). "
             f"Blue outline: >= 0.1 photon. F1=3/2 levels seen: {len(seen15)} of 8."),
            (f"Black number right of each J=2 level: population left in that level after CU2, units 1e-4 (155 V/cm, "
             f"(12.5, 200) MHz, 32 + 43 mW, triplet start, mean over |vx| <= 1 m/s and |y0| <= 8.5 mm). "
             f"J=2 total {100 * sum(res.values()):.2f} %, of it F1=3/2 {100 * f15 / sum(res.values()):.0f} %."),
            ("The field turns 90 deg from CU2 (vertical) to DET (along k) adiabatically without passing zero: every "
             "level keeps its label and mF. Energies to scale within each manifold (MHz axes); levels closer than 2 % of "
             "the band height are spread apart for legibility."),
        ]
        if txt is not None:
            key_on[1] = txt["start"]
            key_on[3] = txt["resid"].format(ratio=line_ratio(line, power_mW, txt))
        key_off = [
            ("E = 0: the quantization axis is the earth field (0.5 G, vertical = CU2 field axis, ⊥ k). σ+ about k is "
             "½ π + ¼ σ+ + ¼ σ− about this axis: arrows ΔmF = 0 (thick) and ±1 (thin) " + arrow_text +
             f". Parity is good at E = 0: only the P'={'+' if trans.P_excited > 0 else '−'} Λ-doublet partner couples "
             f"to J=2 (P=+); the other partner is not drawn. The excited mF' of each level are "
             "degenerate: one line per (F1', F', P'), no mF resolution. The 20 J=2 levels lie within "
             f"{np.ptp([lv.energy_MHz for lv in ground]):.2f} MHz and keep pure F1 (no Stark mixing)."),
            key_on[1],
            (f"Blue number left of each J=2 level: photons per molecule starting in that level, probe on the line, vx = 0, "
             f"{power_mW:g} mW, flat top 1.5 x 1.5 cm, {1e6 * T.T_PROBE:.1f} us at {T.VZ:g} m/s (probe OBE with the axis along k, where the light is pure "
             f"σ+; each level rotated onto k with |d^F(π/2)|², as in the plates-off detection table). "
             f"Blue outline: >= 0.1 photon. F1=3/2 levels seen: {len(seen15)} of 8."),
            key_on[3],
            ("The field falls from CU2 to zero adiabatically: every level keeps its label and mF about the vertical. "
             "Energies to scale within each manifold (MHz axes); levels closer than 2 % of the band height are spread "
             "apart for legibility."),
        ]
        key = key_on if on else key_off
        if notext:
            key = []
        wrapped = [w_ for par in key for w_ in textwrap.wrap(par, 200)]
        bottom = 0.05 + 0.02 * len(wrapped) if wrapped else 0.04
        fig.subplots_adjust(bottom=bottom)
        for k, ln in enumerate(wrapped):
            fig.text(0.02, bottom - 0.05 - 0.02 * k, ln, ha="left", va="center", fontsize=11, color=INK)
        sfx = "" if power_mW == 15 else f"_{power_mW:g}mW"
        name = f"DET_plates154_level_diagram_{short}{sfx}" if on else f"DET_platesoff_level_diagram_{short}{sfx}"
        if txt is not None:
            name += "_cu2test217"
        if notext:
            name += "_notext"
        for ext in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("wrote", name, "| F1=3/2 levels seen:", seen15)


if __name__ == "__main__":
    notext = "notext" in sys.argv
    args = [a for a in sys.argv[1:] if a != "notext"]
    plot(float(args[0]) if len(args) > 0 else 15, args[1] if len(args) > 1 else "on",
         tuple(args[2].split(",")) if len(args) > 2 else DEFAULT_LINES,
         args[3] if len(args) > 3 else None, notext)
