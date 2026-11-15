"""CU2 test: residual per group with F1=3/2 empty before CU2 (report) vs 3 % in F1=3/2 (P(2) F'=1 detection).

    python scripts/plot_cu2_f32_start.py [x_percent]     # default 3
    -> docs/state_cleanup/figures/cu2_residual_lines/f32_start_0_vs_3.{png,pdf}

Same solver as scripts/cu2_test_line_combinations.py, which uses only the shutter on/off ratios R_l (flux drifts
cancel shot by shot) plus S_B/S_A of two adjacent runs for f_t. The F1=3/2 share before CU2 is split over its
8 sublevels uniformly (F=2 : F=1 = 5 : 3); f_t is re-derived analytically from S_B/S_A with that share included.
(a) J=2 left per group after CU2 (solved from the four lines, +-1 sigma MC) and the CU2 OBE prediction, both starts;
(b) on/off ratio of each line, measured vs the CU2 OBE prediction for both starts (no solving).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import cu2_test_line_combinations as L  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

X = float(sys.argv[1]) / 100 if len(sys.argv) > 1 else 0.03
W32 = (5 / 8, 3 / 8)
M = L.Model(24.0)
EOFF = M.eta_off()
PARTS = M.parts("res", "obe")
G = L.G


def f_of_factory(x):
    a, b = EOFF[0], EOFF[1]
    a32, b32 = W32[0] * a[2] + W32[1] * a[3], W32[0] * b[2] + W32[1] * b[3]

    def f_of(r):
        ft = ((1 - x) * (b[1] - r * a[1]) + x * (b32 - r * a32)) / ((1 - x) * (r * (a[0] - a[1]) - (b[0] - b[1])))
        return np.array([ft * (1 - x), (1 - ft) * (1 - x), x * W32[0], x * W32[1]])
    return f_of


def case(x):
    f_of = f_of_factory(x)
    f = f_of(L.RBA)
    n = L.solve(L.eta_on_from(PARTS, f), EOFF, f, L.R)
    S = L.run_mc(M, PARTS, EOFF, f_of, nmc=10000)
    resid = sum(w * M.pv(g, "res") for g, w in zip(G, f) if w > 0)            # [vx, level]
    nm = np.array([resid.sum(0)[M.gi[g]].sum() for g in G])
    Rm = np.einsum("lvi,vi->l", M.Y, resid) / (EOFF @ f)
    return dict(f=f, n=n, err=S.std(0), tot=n.sum(), tot_err=S.sum(1).std(), model=nm, model_tot=nm.sum(), Rm=Rm)


C0, C3 = case(0.0), case(X)
lab = ["(5/2, 3)\ntriplet", "(5/2, 2)\nsinglet", "(3/2, 2)", "(3/2, 1)", "$J=2$\ntotal"]
col0, col3 = "#1f6fb4", "#d1495b"
xs = np.arange(5)
fig, ax = plt.subplots(1, 2, figsize=(13, 4.9), gridspec_kw=dict(width_ratios=[1.25, 1]))

a = ax[0]
for C, dx, col, name in [(C0, -0.13, col0, "F1=3/2 empty before CU2 (report)"),
                         (C3, +0.13, col3, f"{100 * X:g} % in F1=3/2 before CU2 (5:3 over F=2, F=1)")]:
    y = 100 * np.append(C["n"], C["tot"]); e = 100 * np.append(C["err"], C["tot_err"])
    m = 100 * np.append(C["model"], C["model_tot"])
    a.errorbar(xs + dx, y, e, fmt="o", color=col, capsize=3, ms=6, label=f"solved from 4 lines: {name}")
    a.plot(xs + dx, m, "D", color=col, mfc="white", ms=7, mew=1.6, label=f"CU2 OBE, {name.split(' before')[0]}")
    for xi, yi in zip(xs + dx, y):
        a.text(xi, -0.95, f"{yi:.2f}", ha="center", fontsize=7.5, color=col)
a.axhline(0, color="k", lw=0.6)
a.axvline(3.5, color="0.7", lw=0.8)
a.set_xticks(xs, lab)
a.set(ylabel="left after CU2 (% of all $J=2$ before CU2)", ylim=(-1.1, 2.0),
      title="(a) residual per group: numbers below = solved central values")
a.legend(fontsize=7.5, loc="upper left")

b = ax[1]
lines = ["A  R(2) 7/2 F'=4", "B  R(2) 7/2 F'=3", "C  R(2) 5/2 F'=3", "D  R(2) 5/2 F'=2"]
xl = np.arange(4)
b.errorbar(xl, 100 * L.R, 100 * L.RERR, fmt="ks", capsize=3, ms=6, label="measured (shutter on/off)")
b.plot(xl - 0.12, 100 * C0["Rm"], "D", color=col0, mfc="white", ms=7, mew=1.6, label="CU2 OBE, F1=3/2 empty")
b.plot(xl + 0.12, 100 * C3["Rm"], "D", color=col3, mfc="white", ms=7, mew=1.6, label=f"CU2 OBE, {100 * X:g} % F1=3/2")
for xi, r0, r3 in zip(xl, C0["Rm"], C3["Rm"]):
    b.text(xi, 4.35, f"{100 * r0:.2f} → {100 * r3:.2f}", ha="center", fontsize=8)
b.set_xticks(xl, lines, fontsize=8)
b.set(ylabel="on/off ratio $R_l$ (%)", ylim=(-1, 4.8), title=f"(b) line on/off ratios, no solving (model: 0 → {100 * X:g} %)")
b.axhline(0, color="k", lw=0.6)
b.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(0, 0.9))

fig.suptitle("CU2 test 11-08 (217 V/cm, (13, 209) MHz, 20 + 30 mW, 156 m/s, DET 24 mW plates off): "
             f"effect of {100 * X:g} % of $J=2$ in $F_1=3/2$ before CU2.\n"
             "Only shutter on/off ratios enter (flux-drift free); $f_t$ from $S_B/S_A$ of adjacent runs.", fontsize=9)
fig.tight_layout()
out = ROOT / "docs" / "state_cleanup" / "figures" / "cu2_residual_lines" / f"f32_start_0_vs_{100 * X:g}"
for ext in ("png", "pdf"):
    fig.savefig(f"{out}.{ext}", dpi=160)
print(out.with_suffix(".png"))
for k, C in (("0", C0), ("x", C3)):
    print(k, "f_t", C["f"][0] / (C["f"][0] + C["f"][1]), "n", np.round(100 * C["n"], 3), "tot", round(100 * C["tot"], 3),
          "+-", round(100 * C["tot_err"], 3), "model", np.round(100 * C["model"], 3), round(100 * C["model_tot"], 3),
          "Rm", np.round(100 * C["Rm"], 3))
