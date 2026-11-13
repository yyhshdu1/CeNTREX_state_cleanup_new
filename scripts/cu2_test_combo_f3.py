"""F1=5/2 F=3 left after CU2, solved from each line combination: the values, and the (tiny) differences from line A.

    python -I scripts/cu2_test_combo_f3.py
    -> docs/state_cleanup/figures/cu2_residual_lines/fig4_combo_f3.png, combo_f3.json

Every combination contains line A, and A's equation is
    eta_on[A,F=3] n_F=3 + eta_on[A,F=2] n_F=2 + eta_on[A,3/2 F=2] n_3/2F=2 + ... = R_A (eta_off @ f)[A],
so the combinations differ only by which other n_F are allowed to be non-zero:
    n_F=3(combo) - n_F=3(A) = - sum_{F != F=3} eta_on[A,F] n_F(combo) / eta_on[A,F=3].
The differences are evaluated with the same Monte Carlo draws for every combination (paired), so their spread is
the error of the difference, not of the values.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cu2_test_line_combinations as L  # noqa: E402

COMBOS = {"A": (0,), "A+B": (0, 1), "A+C": (0, 2), "A+B+C": (0, 1, 2), "A+B+C+D": (0, 1, 2, 3)}
COLORS = {"A": "#9ecae1", "A+B": "#6baed6", "A+C": "#f4a261", "A+B+C": "#3182bd", "A+B+C+D": "#08306b"}
NMC = 20000


def main():
    M = L.Model(24.0)
    eoff = M.eta_off()
    ft0 = L.ft_of(eoff, L.RBA)
    f0 = np.array([ft0, 1 - ft0, 0, 0])
    parts = M.parts("res", "obe")
    eon = L.eta_on_from(parts, f0)
    model_f3 = sum(w * M.pv(g, "res").sum(0) for g, w in zip(L.G, f0) if w > 0)[M.gi[L.G[0]]].sum()
    n = {c: L.solve(eon, eoff, f0, L.R, ls) for c, ls in COMBOS.items()}
    S = {c: L.run_mc(M, parts, eoff, L.lens_f(eoff), ls, nmc=NMC, seed=1) for c, ls in COMBOS.items()}  # same draws
    out = {"eta_on_A": eon[0].tolist(), "model_F3": float(model_f3), "combos": {}}
    for c in COMBOS:
        d = n[c][0] - n["A"][0]
        terms = {L.GN[g]: float(-eon[0, j] * n[c][j] / eon[0, 0]) for j, g in enumerate(L.G) if j > 0 and n[c][j] != 0}
        out["combos"][c] = {"n_F3": float(n[c][0]), "n": n[c].tolist(), "err": float(S[c][:, 0].std()),
                            "diff_from_A": float(d), "diff_err": float((S[c][:, 0] - S["A"][:, 0]).std()),
                            "terms": terms, "corr_with_A": float(np.corrcoef(S[c][:, 0], S["A"][:, 0])[0, 1])}
    (L.OUT / "combo_f3.json").write_text(json.dumps(out, indent=1))
    figure(out)
    print(json.dumps(out, indent=1))


def figure(out):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9.5,
                         "savefig.facecolor": "white"})
    names = list(COMBOS)
    x = np.arange(len(names))
    C = out["combos"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.1), gridspec_kw={"width_ratios": [1, 1.35]},
                                 constrained_layout=True)
    # (a) values with their full error
    for k, c in enumerate(names):
        a1.errorbar(k, 100 * C[c]["n_F3"], 100 * C[c]["err"], fmt="o", color=COLORS[c], ms=6, capsize=0, lw=1.6)
    a1.axhline(100 * out["model_F3"], color="#d62728", ls="--", lw=1,
               label=f"CU2 OBE model, vz 156 m/s: {100 * out['model_F3']:.3f} %")
    a1.axhline(0, color="#999", lw=0.6)
    a1.set_xticks(x, names)
    a1.set_xlim(-0.6, len(names) - 0.4)
    a1.set_ylim(-0.3, 0.5)
    a1.set_ylabel("F1=5/2 F=3 left after CU2 (% of J=2 start)")
    lo, hi = min(C[c]["n_F3"] for c in names), max(C[c]["n_F3"] for c in names)
    a1.set_title(f"(a) value ± total error: {100 * lo:.4f}–{100 * hi:.4f} %, error ± {100 * C['A']['err']:.2f} %",
                 fontsize=9.5, loc="left")
    a1.legend(frameon=False, fontsize=8.5, loc="upper left")
    # (b) difference from line A alone, paired MC error
    u = 1e3 * 100                                     # 10^-3 % of the J=2 start
    for k, c in enumerate(names):
        a2.errorbar(k, u * C[c]["diff_from_A"], u * C[c]["diff_err"], fmt="o", color=COLORS[c], ms=6, capsize=3,
                    lw=1.6)
        a2.annotate(f"{u * C[c]['diff_from_A']:+.2f}", (k, u * C[c]["diff_from_A"]), xytext=(9, 0),
                    textcoords="offset points", va="center", fontsize=8.5, color="#262624")
    a2.axhline(0, color="#999", lw=0.6)
    a2.set_xticks(x, names)
    a2.set_xlim(-0.6, len(names) - 0.4)
    a2.set_ylim(-6.5, 1.5)
    a2.set_ylabel("n_F=3(combination) − n_F=3(A)   (10⁻³ % of J=2)")
    a2.set_title(f"(b) shift from line A alone (error of (a) = ±{u * C['A']['err']:.0f} in these units)",
                 fontsize=9.5, loc="left")
    fig.suptitle("F1=5/2 F=3 (triplet) left after CU2, NLTL2 2068: every combination is solved from line A; the "
                 "other lines only change it through line A's small crosstalk", fontsize=10)
    fig.savefig(L.OUT / "fig4_combo_f3.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
