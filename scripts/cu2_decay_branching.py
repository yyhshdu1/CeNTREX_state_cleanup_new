"""CU2 excited-state decay branching to the X-state J levels, drawn like the CU1 slide.

CU2 drives R(2) F1'=7/2 F'=3 and F'=4 (J=2+ -> J'=3-). Branching ratios to each X J level (v = 0)
come from centrex-tlf B-state eigenstates (J' = 1-5) at E = 0 and at the CU2 field, 155 V/cm:

    a) E = 0. J'=3- decays only to J=2+ and J=4+. The same holds at 155 V/cm for mF'=0, which has
       no linear Stark mixing with its Lambda-doublet partner.
    b) 155 V/cm, |mF'| >= 1. The Lambda doublet is Stark-mixed, so both mixed eigenstates of each
       |mF'| are laser-driven; they also decay to J=3- and J=5-. The ranges cover both mixed
       eigenstates, every |mF'| >= 1, and both F'=3 and F'=4.

J=1- is drawn in both panels. J'=3 cannot decay to J=1 directly (Delta J = -2); it gets there only
through the small hyperfine J'=2 admixture of the P'=+ partner. At E = 0 the driven P'=- state has
exactly zero branching to J=1 (parity); at 155 V/cm it is <= 2.6e-5 per decay (F'=3; <= 1e-6 for
F'=4). Small components are kept down to 1e-7 in amplitude so that these numbers are not truncated.

As a check, the same code reproduces the CU1 numbers (R(1) F1'=5/2 F'=2, 150 V/cm): 0.525 / 0.475
at E = 0, and 0.17-0.35, 0.31-0.64, 0.16-0.32, 0.02-0.04 to J = 1-4 at 150 V/cm.

    python scripts/cu2_decay_branching.py [E]        # prints the table, writes the slide figure (default 217 V/cm)
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # use the repo library, not the older copy in the conda env

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from centrex_tlf import couplings, hamiltonian, states  # noqa: E402

OUT = ROOT / "docs" / "state_cleanup" / "figures" / "slides"
E_CU2 = float(sys.argv[1]) if len(sys.argv) > 1 else 217.0   # V/cm; measured CU2 field (2026-11-05)
NAME = f"CU2_decay_branching_R2_F3_F4_{E_CU2:g}Vcm"
B_FIELD = np.array([0.0, 0.0, 1e-3])
JMAX_X = 6
TOL = 1e-7                      # keep small state components: the J=1 branching is ~1e-5
X, B = states.ElectronicState.X, states.ElectronicState.B


def branching(Jp, F1p, Fp, Pp, E):
    """{(P label, mF'): {J: BR}} for both Lambda-doublet partners of (J', F1', F') at field E (V/cm)."""
    warnings.simplefilter("ignore")
    ground = [1 * g for g in states.generate_coupled_states_X(
        states.QuantumSelector(J=list(range(JMAX_X + 1)), electronic=X))]
    gJ = np.array([g.largest.J for g in ground])
    approx = []
    for P in (Pp, -Pp):
        approx += list(states.generate_coupled_states_B(states.QuantumSelector(
            J=Jp, F1=F1p, F=Fp, P=P, Ω=1, electronic=B)))
    red = hamiltonian.generate_reduced_B_hamiltonian(approx, E=np.array([0.0, 0.0, E]), B=B_FIELD,
                                                     Jmin=1, Jmax=Jp + 2, stol=TOL)
    out = {}
    for a, q in zip(approx, red.QN_basis):
        L = q.largest
        assert (L.J, L.F1, L.F, L.mF, L.P) == (a.J, a.F1, a.F, a.mF, a.P), (a, L)
        br = couplings.calculate_br(q, ground, tol=TOL)
        out[(a.P, a.mF)] = {J: float(br[gJ == J].sum()) for J in range(JMAX_X + 1)}
    return out


def ranges(tables, P_laser, which):
    """Min/max BR per J over the selected eigenstates. which = "mF0" (driven partner, mF'=0) or
    "mixed" (both partners, |mF'| >= 1)."""
    rows = [br for tab in tables for (P, mF), br in tab.items()
            if (which == "mF0" and mF == 0 and P == P_laser) or (which == "mixed" and mF != 0)]
    return {J: (min(r[J] for r in rows), max(r[J] for r in rows)) for J in range(JMAX_X + 1)}


# ---- drawing -----------------------------------------------------------------
INK, GREY, LIGHT = "black", "0.55", "0.62"


def wavy(ax, p0, p1, amp=0.006, waves=None, lw=0.9):
    """Wavy decay arrow from p0 to p1 in axes coordinates (aspect corrected)."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    asp = ax.figure.get_figheight() * ax.get_position().height / (ax.figure.get_figwidth() * ax.get_position().width)
    scale = np.array([1.0, asp])
    d = (p1 - p0) / scale                     # isotropic coordinates
    length = np.hypot(*d)
    u, n = d / length, np.array([-d[1], d[0]]) / length
    head = 0.018
    s = np.linspace(0.0, length - head, 600)
    waves = waves or max(8, int(length / 0.03))
    w = amp * np.sin(2 * np.pi * waves * s / (length - head))
    pts = (p0 / scale + np.outer(s, u) + np.outer(w, n)) * scale
    ax.plot(pts[:, 0], pts[:, 1], color=INK, lw=lw, transform=ax.transAxes, clip_on=False)
    ax.annotate("", xy=p1, xytext=tuple(pts[-1]), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.18", color=INK, lw=lw,
                                shrinkA=0, shrinkB=0))


def fmt_range(lo, hi):
    a, b = f"{lo:.2f}", f"{hi:.2f}"
    return a if a == b else f"{a}–{b}"


def panel(ax, letter, values, reached, footer):
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(-0.02, 0.98, letter, transform=ax.transAxes, fontsize=17, fontweight="bold", va="top")
    top = (0.45, 0.90)
    ax.text(*top, r"$|\widetilde{J}'=3^-,\,\widetilde{F}_1'=7/2,\,F'=3,4\rangle$", fontsize=16,
            ha="center", va="center", transform=ax.transAxes)
    ax.text(0.90, 0.90, r"$B\,{}^3\Pi_1$", fontsize=17, ha="center", va="center", transform=ax.transAxes)
    ax.text(0.93, 0.21, r"$X\,{}^1\Sigma^+$", fontsize=17, ha="center", va="center", transform=ax.transAxes)
    # ground J levels, higher J sits higher (rotational energy)
    ground = {1: (0.07, 0.10, "-"), 2: (0.27, 0.17, "+"), 3: (0.47, 0.24, "-"), 4: (0.66, 0.32, "+"),
              5: (0.84, 0.43, "-")}
    for J, (x, y, p) in ground.items():
        ax.text(x, y, rf"$|J={J}^{p}\rangle$", fontsize=16, ha="center", va="center",
                color=INK if J in reached else LIGHT, transform=ax.transAxes)
    # laser: thick grey arrow J=2+ -> excited
    x0, y0, x1, y1 = 0.255, 0.215, 0.315, 0.81
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>,head_length=1.4,head_width=0.65", color=GREY, lw=5.5,
                                shrinkA=0, shrinkB=0))
    ang = np.degrees(np.arctan2((y1 - y0) * ax.figure.get_figheight() * ax.get_position().height,
                                (x1 - x0) * ax.figure.get_figwidth() * ax.get_position().width))
    ax.text(0.262, 0.52, "R2 F3, F4 (CU2)", rotation=ang, fontsize=13, family="serif",
            ha="center", va="center", transform=ax.transAxes)
    # decays: start points spread along the upper label, ends just above each ground label
    starts = {1: (0.24, 0.82), 2: (0.37, 0.82), 3: (0.47, 0.82), 4: (0.57, 0.82), 5: (0.67, 0.82)}
    frac = {1: 0.45, 2: 0.55, 3: 0.62, 4: 0.55, 5: 0.40}      # where along the arrow the number sits
    for J in reached:
        x, y, _ = ground[J]
        end = (x + 0.005, y + 0.045)
        wavy(ax, starts[J], end, lw=0.6 if J == 1 else 0.9)
        f = frac[J]
        # J=5 arrow is shallow: number above it; J=1 arrow runs to the left: number on its left
        dx, dy = {5: (0.0, 0.045), 1: (-0.018, 0.0)}.get(J, (0.018, 0.0))
        lx = starts[J][0] + f * (end[0] - starts[J][0]) + dx
        ly = starts[J][1] + f * (end[1] - starts[J][1]) + dy
        ax.text(lx, ly, values[J], fontsize=15, family="serif", ha="right" if J == 1 else "left", va="center",
                transform=ax.transAxes)
    ax.text(0.0, 0.02, footer, fontsize=13, family="serif", ha="left", va="bottom", transform=ax.transAxes)


def main():
    F_EXC = (3, 4)
    tab0 = [branching(3, 3.5, F, -1, 0.0) for F in F_EXC]
    tabE = [branching(3, 3.5, F, -1, E_CU2) for F in F_EXC]
    r0 = ranges(tab0, -1, "mF0")
    rE0 = ranges(tabE, -1, "mF0")          # mF'=0 at the CU2 field: same as E = 0
    rE = ranges(tabE, -1, "mixed")

    print(f"J   E=0 (F'=3,4)      {E_CU2:g} V/cm mF'=0    {E_CU2:g} V/cm |mF'|>=1 (both partners, F'=3,4)")
    for J in range(JMAX_X + 1):
        if max(r0[J][1], rE[J][1]) > 1e-3:
            print(f"{J}   {r0[J][0]:.3f}-{r0[J][1]:.3f}      {rE0[J][0]:.3f}-{rE0[J][1]:.3f}      "
                  f"{rE[J][0]:.3f}-{rE[J][1]:.3f}")
    for F, t in zip(F_EXC, tabE):
        print(f"F'={F} at {E_CU2:g} V/cm, |mF'|>=1:",
              {J: fmt_range(*ranges([t], -1, "mixed")[J]) for J in (2, 3, 4, 5)},
              f"J=1 max {ranges([t], -1, 'mixed')[1][1]:.2e}")
    print(f"J=1: E=0 driven state {r0[1][1]:.1e}, {E_CU2:g} V/cm mF'=0 {rE0[1][1]:.1e}, |mF'|>=1 {rE[1][1]:.1e}")

    plt.rcParams.update({"text.usetex": False, "mathtext.fontset": "cm", "font.family": "serif",
                         "font.serif": ["Times New Roman", "DejaVu Serif"]})
    fig, (a, b) = plt.subplots(1, 2, figsize=(17.2, 7.4))
    fig.subplots_adjust(left=0.02, right=0.99, top=0.97, bottom=0.03, wspace=0.08)
    vals0 = {J: f"{0.5 * (r0[J][0] + r0[J][1]):.2f}" for J in (2, 4)}
    panel(a, "a)", vals0, (2, 4), rf"$E=0$  (also $m_F'=0$ at {E_CU2:g} V/cm)")
    valsE = {J: fmt_range(*rE[J]) for J in (2, 3, 4, 5)}
    m, e = f"{rE[1][1]:.0e}".split("e")
    valsE[1] = rf"$\leq {int(m)}\times10^{{{int(e)}}}$"
    panel(b, "b)", valsE, (1, 2, 3, 4, 5), rf"$E={E_CU2:g}$ V/cm,  $|m_F'|\geq 1$")
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{NAME}.{ext}", dpi=200)
    print("wrote", OUT / f"{NAME}.png")
    summary = {"E0_mF0": {J: r0[J] for J in (1, 2, 4)}, f"E{E_CU2:g}_mixed": {J: rE[J] for J in (1, 2, 3, 4, 5)},
               f"E{E_CU2:g}_mF0": {J: rE0[J] for J in (1, 2, 4)}}
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
