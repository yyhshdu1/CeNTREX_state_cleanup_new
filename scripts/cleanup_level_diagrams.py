"""Field-dressed level diagrams of CU1 and CU2 with the levels that most need cleaning marked.

    python scripts/cleanup_level_diagrams.py compute     # residual per starting sublevel, CU1 and CU2 (R1, R2)
    python scripts/cleanup_level_diagrams.py plot        # diagrams (slides style) from the cached residuals

The diagrams are utils.plotting.plot_transition_level_diagram (levels segmented by zero-field parent
character), with three overlays on the ground manifold:
    number      population left in the cleaned J after the stage, for a molecule starting in that level
    yellow      SPB-related level (CU1: J=1 |m_J|=1 levels SPB transfers into; CU2: J=2 m_J=0 levels it
                transfers from), identified from the m_J content in the uncoupled basis
    outline     residual above twice the average over the manifold
OBE eigenstates are matched to the diagram levels by energy rank inside each m_F block (same Hamiltonian,
same fields), so the diagram's adiabatically tracked labels are used throughout.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cleanup_spb_doppler import (CU1_DETUNING, CU1_E, CU1_SOLVER, CU2_POINTS, build_cu1,  # noqa: E402
                                         describe, solve_states)
from scripts.cu2_laser_aom_scan import B_FIELD, CONFIRM_SOLVER, log  # noqa: E402
from scripts.cu2_laser_aom_scan import build as build_cu2  # noqa: E402

CACHE = ROOT / "tmp" / "level_diagrams" / "residuals.json"
FIG_DIR = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
CU2_E = 200.0
MIN_SEP = 0.02          # cosmetic minimum vertical gap between levels of one m_F column (axes fraction)


def ground_table(obe, gidx):
    """m_F, energy (MHz, own origin) and m_J content of the ground levels of an OBE system."""
    import sympy as smp
    H = obe.H_symbolic
    rows = []
    for i in gidx:
        hii = H[i, i]
        e = complex(smp.N(hii.subs({s: 0 for s in hii.free_symbols}))).real / (2 * np.pi * 1e6)
        d = describe(obe.QN[i])
        rows.append(dict(idx=int(i), mF=int(round(d["mF"])), F1=float(d["F1"]), F=int(d["F"]), energy=e,
                         p_mj=d["p_mj"]))
    return rows


def rank_in_block(rows):
    """Rank of each row by energy inside its m_F block."""
    for mF in {r["mF"] for r in rows}:
        block = sorted((r for r in rows if r["mF"] == mF), key=lambda r: r["energy"])
        for k, r in enumerate(block):
            r["rank"] = k
    return rows


def compute():
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    out = {}
    from centrex_tlf import states
    X = states.ElectronicState.X

    # ---- CU1: all 12 J=1 sublevels -------------------------------------------------------------------
    prepared1, gidx1, obe1 = build_cu1(return_obe=True)
    rows1 = rank_in_block(ground_table(obe1, gidx1))
    t0 = time.perf_counter()
    pops = solve_states(prepared1, len(obe1.QN), [r["idx"] for r in rows1], ["laser_detuning"],
                        [np.array([2 * np.pi * 1e6 * CU1_DETUNING])] * len(rows1), CU1_SOLVER, threads=12)
    log(f"CU1: {len(rows1)} trajectories in {time.perf_counter() - t0:.0f} s")
    for r, p in zip(rows1, pops):
        r["residual"] = float(p[gidx1].sum())
        r["spb"] = r["p_mj"].get(1, 0) + r["p_mj"].get(-1, 0) > 0.5
    out["cu1"] = dict(E=CU1_E, detuning_MHz=CU1_DETUNING, levels=rows1)

    # ---- CU2: m_F >= 0 sublevels at R1 and R2, m_F < 0 from the +|m_F| partner ------------------------
    prepared2, _, gidx2, info, obe2 = build_cu2(return_obe=True)
    rows2 = rank_in_block(ground_table(obe2, gidx2))
    pos = [r for r in rows2 if r["mF"] >= 0]
    j1 = np.asarray(states.QuantumSelector(J=1, electronic=X).get_indices(obe2.QN), dtype=int).ravel()
    idx, par, tags = [], [], []
    for reg, (d0, faom) in CU2_POINTS.items():
        d1 = d0 + faom - info["delta_main_MHz"]
        for r in pos:
            idx.append(r["idx"]); par.append(2 * np.pi * 1e6 * np.array([d0, d1])); tags.append((reg, r["idx"]))
    t0 = time.perf_counter()
    pops = solve_states(prepared2, len(obe2.QN), idx, ["detuning0", "detuning1"], par, CONFIRM_SOLVER, threads=12)
    log(f"CU2: {len(idx)} trajectories in {time.perf_counter() - t0:.0f} s")
    res = {tag: (float(p[gidx2].sum()), float(p[j1].sum())) for tag, p in zip(tags, pops)}
    by_key = {(r["mF"], r["rank"]): r for r in rows2}
    for r in rows2:
        partner = by_key[(abs(r["mF"]), r["rank"])]
        r["residual"] = {reg: res[(reg, partner["idx"])][0] for reg in CU2_POINTS}
        r["to_J1"] = {reg: res[(reg, partner["idx"])][1] for reg in CU2_POINTS}
        r["spb"] = r["p_mj"].get(0, 0) > 0.5
    out["cu2"] = dict(E=CU2_E, points_MHz={k: list(v) for k, v in CU2_POINTS.items()},
                      delta_main_MHz=info["delta_main_MHz"], levels=rows2)
    CACHE.write_text(json.dumps(out, indent=1, default=float))
    for r in sorted(rows1, key=lambda r: -r["residual"]):
        log(f"CU1 F1={r['F1']} F={r['F']} mF={r['mF']:+d}: residual {r['residual']:.2e}{'  SPB' if r['spb'] else ''}")
    for r in sorted(rows2, key=lambda r: -r["residual"]["R2"]):
        log(f"CU2 F1={r['F1']} F={r['F']} mF={r['mF']:+d}: residual R1 {r['residual']['R1']:.2e} "
            f"R2 {r['residual']['R2']:.2e}; to J=1 R2 {r['to_J1']['R2']:.2e}{'  SPB' if r['spb'] else ''}")


def plot():
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, Rectangle
    from matplotlib.text import Annotation
    from scipy.special import jv
    from centrex_tlf import transitions
    from centrex_tlf.utils import plotting as lp

    data = json.loads(CACHE.read_text())
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    YELLOW, INK = "#f6e3a1", "#111111"
    BLUE, RED = "#1537a6", "#d0021b"          # laser colours: F'=3 beam and CU1 blue, F'=4 beam red
    F_EOM, BETA = 17.5, 1.5
    WINDOW = 4 * 1.6076                       # a comb line "targets" drivable sublevels within 4 Gamma/2pi
    width, pad = 0.48, 0.04                   # bar width and box margin (m_F units)

    def as_frac(v):
        return f"{int(round(2 * v))}/2" if abs(v - round(v)) > 0.25 else f"{v:.0f}"

    def fmt(v):
        return f"{v:.0f}" if v >= 10 else f"{v:.1f}" if v >= 1 else f"{v:.2f}"

    def mlabel(m):
        return "m'=0" if m == 0 else f"m'=±{abs(m)}"

    def figure(stage, trans, E, residual_of, point_text, fname, title, comb, extra=()):
        """comb: main_g (F1,F,mF), main_e (F1',F',mF',P'), parity (addressed P'), lines [(n, power, delta)],
        color, name."""
        levels = data[stage]["levels"]
        structure = lp.calculate_transition_level_structure(trans, E=E, B=float(B_FIELD[2]))
        diag = lp.plot_transition_level_diagram(structure=structure, title=title, figsize=(19.5, 11.2),
                                                min_level_separation=MIN_SEP, level_width=width)
        ax, fig = diag.ax, diag.fig
        ground, excited = list(structure.ground), list(structure.excited)
        J = trans.J_ground
        c = comb["color"]

        # the generic R(J) arrow of the diagram is replaced by the laser comb below
        for t in list(ax.texts):
            if (isinstance(t, Annotation) and t.get_text() == "") or (
                    "\n" in t.get_text() and t.get_transform() == ax.transData):
                t.remove()

        e_main = next(k for k, lv in enumerate(excited) if (lv.F1, lv.F, lv.mF, lv.P) == comb["main_e"])
        g_main = next(lv for lv in ground if (lv.F1, lv.F, lv.mF) == comb["main_g"])
        lines = [(n, pw, d, excited[e_main].energy_MHz + d) for n, pw, d in comb["lines"]]
        labels = [(f"laser   δ = {d:+.1f} MHz" if len(lines) == 1 else
                   f"{'carrier' if n == 0 else f'n = {n:+d}'}   {100 * pw:.0f} %   δ = {d:+.1f} MHz")
                  for n, pw, d, _ in lines]

        # room on the right for the comb labels (excited band) and the box note (ground band): move the
        # F / F1 key. A text of axes-width fraction w starting `start` units right of mf_max needs a shift
        # G = (start + margin + w * range) / (1 - w), since the shift itself widens the data range.
        mf_max = max(abs(lv.mF) for lv in ground + excited)
        mfg = max(abs(lv.mF) for lv in ground)
        energies = np.array([lv.energy_MHz for lv in ground])
        box_note = (f"dashed box: all {len(ground)} J={J}\nsublevels lie within {np.ptp(energies):.1f} MHz\n"
                    f"(Γ/2π = 1.6 MHz), not\nresolved by the laser")
        renderer = fig.canvas.get_renderer()
        ax_w = ax.get_window_extent(renderer).width

        def frac(texts, fs):
            probe = [ax.text(0, 0, s, fontsize=fs) for s in texts]
            out = max(t.get_window_extent(renderer).width for t in probe) / ax_w
            for t in probe:
                t.remove()
            return out

        x0, x1 = ax.get_xlim()
        need = [(width / 2 + 0.33, frac(labels, 10)),
                (mfg - mf_max + width / 2 + 0.7, frac(box_note.split("\n"), 10.5))]
        gutter = max((start + 0.25 + w * (x1 - x0)) / (1 - w) for start, w in need)
        cut = mf_max + 0.3
        for t in ax.texts:
            if t.get_transform() == ax.transData and t.get_position()[0] > cut:
                t.set_x(t.get_position()[0] + gutter)
        for ln in ax.lines:
            x = np.asarray(ln.get_xdata(), dtype=float)
            if x.size and x.min() > cut:
                ln.set_xdata(x + gutter)
        ax.set_xlim(x0, x1 + gutter)

        # ---- ground: residual numbers, SPB band, bottleneck outline ----------------------------------
        ys = lp._band_y(ground, *lp._Y_GROUND_BAND, MIN_SEP)
        order = {}      # diagram level <-> OBE eigenstate by energy rank inside the m_F block
        for mF in {lv.mF for lv in ground}:
            for k, j in enumerate(sorted((j for j, lv in enumerate(ground) if lv.mF == mF),
                                         key=lambda j: ground[j].energy_MHz)):
                order[j] = (mF, k)
        obe_by = {(r["mF"], r["rank"]): r for r in levels}
        vals = np.array([residual_of(obe_by[order[j]]) for j in range(len(ground))])
        mean = vals.mean()
        for j, (lv, y) in enumerate(zip(ground, ys)):
            r = obe_by[order[j]]
            if (r["F1"], r["F"]) != (lv.F1, lv.F):
                log(f"{stage}: label differs for mF={lv.mF:+d} rank {order[j][1]}: OBE largest "
                    f"({r['F1']}, {r['F']}) vs tracked parent ({lv.F1}, {lv.F})")
            if r["spb"]:
                ax.add_patch(Rectangle((lv.mF - width / 2 - pad, y - 0.009), width + 2 * pad, 0.018,
                                       facecolor=YELLOW, edgecolor="none", zorder=0))
            if vals[j] > 2 * mean:
                ax.add_patch(Rectangle((lv.mF - width / 2 - pad, y - 0.009), width + 2 * pad, 0.018,
                                       facecolor="none", edgecolor=INK, lw=1.8, zorder=3))
            ax.text(lv.mF + width / 2 + pad + 0.03, y, fmt(1e4 * vals[j]), ha="left", va="center",
                    fontsize=10, color=INK, fontweight="bold" if vals[j] > 2 * mean else "normal")

        # dashed box: the laser does not resolve the ground sublevels
        bx0, bx1 = -mfg - width / 2 - 0.15, mfg + width / 2 + 0.55
        by0, by1 = min(ys) - 0.024, max(ys) + 0.024
        ax.add_patch(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fill=False, ls=(0, (5, 3)), lw=1.5,
                               edgecolor="0.3", zorder=4))
        ax.text(bx1 + 0.15, 0.5 * (by0 + by1), box_note, ha="left", va="center", fontsize=10.5, color="0.25",
                linespacing=1.4)

        # ---- excited: comb lines, targets, arrows --------------------------------------------------------
        yl, yh = lp._Y_EXCITED_BAND
        lo, hi = min(e.energy_MHz for e in excited), max(e.energy_MHz for e in excited)

        def y_of(en):           # linear part of the diagram's energy -> height mapping
            return yl + (en - lo) / (hi - lo) * (yh - yl)

        ys_e = lp._band_y(excited, yl, yh, MIN_SEP)
        drivable = [k for k, e in enumerate(excited) if e.character.get(comb["parity"], 0.0) > 0.05]
        xl, xr = -mf_max - width / 2 - 0.1, mf_max + width / 2 + 0.15
        ax.text(excited[e_main].mF + width / 2 + pad + 0.03, ys_e[e_main], "δ = 0", ha="left", va="center",
                fontsize=9.5, color=c, zorder=6, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
        targets = []
        for (n, pw, d, en), lab in zip(lines, labels):
            y = y_of(en)
            strong = pw > 0.2
            ax.plot([xl, xr], [y, y], ls=(0, (6, 3)), lw=1.8 if strong else 1.0, color=c, zorder=1)
            ax.text(xr + 0.08, y, lab, ha="left", va="center", fontsize=10, color=c,
                    fontweight="bold" if strong else "normal")
            hits = {}       # one entry per (|m'|, P') pair inside the window
            for k in drivable:
                tgt = excited[k]
                if abs(tgt.energy_MHz - en) <= WINDOW:
                    key_ = (abs(tgt.mF), tgt.P)
                    if key_ not in hits or abs(tgt.energy_MHz - en) < abs(hits[key_].energy_MHz - en):
                        hits[key_] = tgt
            hits = sorted(hits.values(), key=lambda t: abs(t.energy_MHz - en))
            targets.append((n, [(t, en - t.energy_MHz) for t in hits]))
            for tgt in hits:
                for kk, e in enumerate(excited):
                    if abs(e.mF) == abs(tgt.mF) and e.P == tgt.P:
                        ax.add_patch(Rectangle((e.mF - width / 2 - pad, ys_e[kk] - 0.010), width + 2 * pad,
                                               0.020, facecolor="none", edgecolor=c, lw=2.2 if strong else 1.2,
                                               zorder=3))
                ax.add_patch(FancyArrowPatch((0.0, by1), (abs(tgt.mF), y), arrowstyle="-|>", mutation_scale=16,
                                             lw=2.0 if strong else 1.0, color=c, alpha=0.9 if strong else 0.6,
                                             shrinkA=0, shrinkB=0, zorder=2.5))

        # ---- key below the figure ----------------------------------------------------------------------
        dg = energies - g_main.energy_MHz
        g_lab = f"|F1={as_frac(g_main.F1)}, F={g_main.F}, mF={g_main.mF:+d}>"
        e_lv = excited[e_main]
        e_lab = f"|F1'={as_frac(e_lv.F1)}, F'={e_lv.F}, mF'={e_lv.mF:+d}, P'={'+' if e_lv.P > 0 else '-'}>"
        comb_text = ("single frequency, no phase modulation" if len(lines) == 1 else
                     f"carrier and sidebands n of the {F_EOM} MHz EOM (beta = {BETA}), % = power fraction J_n(beta)^2")
        def tg_one(n, hits):
            name = "laser" if len(lines) == 1 else "carrier" if n == 0 else f"n={n:+d}"
            if not hits:
                return f"{name}: none"
            return f"{name} -> " + " and ".join(f"{mlabel(t.mF)} (P'{'+' if t.P > 0 else '-'}) {off:+.1f} MHz"
                                                for t, off in hits)
        tg = "; ".join(tg_one(n, h) for n, h in targets)
        key = [
            (f"Dashed lines and arrows: {comb['name']}, {comb_text}. Drawn for the main ground level {g_lab}; "
             f"delta = detuning from its resonance with {e_lab} (marked delta = 0)."),
            (f"For the other ground levels the lines move by {dg.min():+.1f} to {dg.max():+.1f} MHz.   "
             f"Coloured outline and arrow: sublevels a line reaches, i.e. within 4 Gamma/2pi = {WINDOW:.1f} MHz and "
             f"with addressed-parity character > 5 % (line minus level): {tg}."),
            (f"Number right of each J={J} level: population left in J={J} after {stage.upper()} for a molecule "
             f"starting in that level, in units of 1e-4 ({point_text}). Average over the {len(ground)} levels: "
             f"{1e4 * mean:.1f}e-4."),
            ("Bold black outline: residual more than twice that average.   Yellow band: SPB-related level "
             + ("(J=1, |m_J|=1 content > 0.5; SPB transfers into these)." if stage == "cu1" else
                "(J=2, m_J=0 content > 0.5; SPB transfers out of these, what it leaves CU2 must remove).")),
            "Single trajectory on the beam axis, no transverse velocity; 13 passes, 10 % loss per pass, X/Z "
            "switching at 0.5 MHz.",
        ] + list(extra)
        import textwrap
        wrapped = [w_ for par in key for w_ in textwrap.wrap(par, 200)]
        bottom = 0.05 + 0.02 * len(wrapped)
        fig.subplots_adjust(bottom=bottom)
        for k, line in enumerate(wrapped):
            fig.text(0.02, bottom - 0.05 - 0.02 * k, line, ha="left", va="center", fontsize=11, color=INK)
        for ext in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"{fname}.{ext}", dpi=150, bbox_inches="tight")
        plt.close(fig)
        log(f"wrote {fname}; targets: {tg}")

    def comb_lines(delta):
        return [(n, float(jv(n, BETA) ** 2), delta + n * F_EOM) for n in (-2, -1, 0, 1, 2)]

    cu1 = data["cu1"]
    figure("cu1", transitions.R1_F1_5o2_F2, cu1["E"], lambda r: r["residual"],
           f"80 mW, delta = +{cu1['detuning_MHz']:.1f} MHz, 4.8 cm beam",
           "CU1_level_diagram_R1_F2_150Vcm_cleanup",
           "CU1  |  R(1) $F_1'=5/2$, $F'=2$  |  E = 150 V/cm  |  B = 1 mG (model)",
           comb=dict(main_g=(1.5, 2, 1), main_e=(2.5, 2, 2, 1), parity=+1, color=BLUE, name="CU1 laser",
                     lines=[(0, 1.0, cu1["detuning_MHz"])]))
    # CU1 at the set point chosen for the transverse distribution (no EOM, +19.0 MHz), residuals per
    # starting level from the on-axis trajectories of scripts/cu1_trajectories.py
    traj = ROOT / "tmp" / "cu1_trajectories" / "cu1_trajectories.npz"
    if traj.exists():
        z = np.load(traj)
        k = list(z["detunings"]).index(19.0)
        gidx = z["ground_idx"]
        res19 = {int(i): float(z["pops"][k, j, -1, gidx].sum()) for j, i in enumerate(gidx)}
        figure("cu1", transitions.R1_F1_5o2_F2, cu1["E"], lambda r: res19[r["idx"]],
               "80 mW, delta = +19.0 MHz, 4.8 cm beam, vx = 0, y0 = 0",
               "CU1_level_diagram_R1_F2_150Vcm_cleanup_19MHz",
               "CU1  |  R(1) $F_1'=5/2$, $F'=2$  |  E = 150 V/cm  |  laser at +19.0 MHz",
               comb=dict(main_g=(1.5, 2, 1), main_e=(2.5, 2, 2, 1), parity=+1, color=BLUE, name="CU1 laser",
                         lines=[(0, 1.0, 19.0)]),
               extra=("+19.0 MHz is the best detuning without EOM for a Gaussian vx distribution of 1.8 m/s "
                      "FWHM (truncated at 1 m/s): averaged over vx and |y0| <= 15 mm it gives 95.2 %, against "
                      "93.1 % at +20.2 MHz, although on the axis at vx = 0 it is worse (96.36 % against 99.965 %).",))
    cu2 = data["cu2"]
    d0, faom = cu2["points_MHz"]["R2"]
    d1 = d0 + faom - cu2["delta_main_MHz"]
    pt = f"setting R2: delta0 = {d0:.1f} MHz, f_AOM = {faom:.1f} MHz; 80 + 24 mW, 2.4 cm beam"
    lv2 = cu2["levels"]
    worst = sorted({(r["F1"], r["F"], abs(r["mF"])): r for r in lv2}.values(),
                   key=lambda r: -r["residual"]["R2"])[:3]
    d0r1, faomr1 = cu2["points_MHz"]["R1"]
    d1r1 = d0r1 + faomr1 - cu2["delta_main_MHz"]
    r1_line = (f"At setting R1 (delta0 = {d0r1:.1f} MHz, f_AOM = {faomr1:.1f} MHz) the same three worst levels give "
               + ", ".join(f"{fmt(1e4 * r['residual']['R1'])} (F1={as_frac(r['F1'])}, F={r['F']}, mF=±{abs(r['mF'])})"
                           if r["mF"] else
                           f"{fmt(1e4 * r['residual']['R1'])} (F1={as_frac(r['F1'])}, F={r['F']}, mF=0)"
                           for r in worst)
               + f"; average {1e4 * np.mean([r['residual']['R1'] for r in lv2]):.1f}e-4.")
    other = (f"Beam 1 (F'=4, delta1 = {d1:+.2f} MHz) lies about 200 MHz higher and is drawn in the F'=4 diagram.",
             f"Beam 0 (F'=3, delta0 = {d0:+.1f} MHz) lies about 200 MHz lower and is drawn in the F'=3 diagram.   "
             f"At R1 (delta1 = {d1r1:+.2f} MHz) the beam-1 lines sit at the same frequencies, each order shifted "
             f"by one.")
    figure("cu2", transitions.R2_F1_7o2_F3, cu2["E"], lambda r: r["residual"]["R2"], pt,
           "CU2_level_diagram_R2_F3_200Vcm_cleanup",
           "CU2 beam 0  |  R(2) $F_1'=7/2$, $F'=3$  |  E = 200 V/cm  |  B = 1 mG (model)",
           comb=dict(main_g=(2.5, 3, 2), main_e=(3.5, 3, 3, -1), parity=-1, color=BLUE,
                     name="CU2 beam 0 at setting R2", lines=comb_lines(d0)),
           extra=(other[0], r1_line))
    figure("cu2", transitions.R2_F1_7o2_F4, cu2["E"], lambda r: r["residual"]["R2"], pt,
           "CU2_level_diagram_R2_F4_200Vcm_cleanup",
           "CU2 beam 1  |  R(2) $F_1'=7/2$, $F'=4$  |  E = 200 V/cm  |  B = 1 mG (model)",
           comb=dict(main_g=(2.5, 3, 0), main_e=(3.5, 4, 1, -1), parity=-1, color=RED,
                     name="CU2 beam 1 at setting R2", lines=comb_lines(d1)),
           extra=(other[1], r1_line))


if __name__ == "__main__":
    {"compute": compute, "plot": plot}[sys.argv[1]]()
