"""LaTeX (matplotlib mathtext) renderings of the equations of the CU2-test line-combination method, for the slides.

    python -I scripts/cu2_test_equations.py
    -> docs/state_cleanup/figures/cu2_residual_lines/eq_*.png  (transparent, 300 dpi, 20 pt)
"""
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "docs" / "state_cleanup" / "figures" / "cu2_residual_lines"
mpl.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif"})
EQ = {
    "eq_signal": r"$S_l = \sum_F \eta_{lF}\, N_F$",
    "eq_ratio": r"$\sum_F \eta^{\mathrm{on}}_{lF}\, n_F \;=\; R_l \sum_F \eta^{\mathrm{off}}_{lF}\, f_F ,\qquad "
                r"R_l = \frac{S^{\mathrm{on}}_l}{S^{\mathrm{off}}_l},\qquad n_F = \frac{N^{\mathrm{on}}_F}"
                r"{N^{\mathrm{off}}_{J=2}},\qquad f_F = \frac{N^{\mathrm{off}}_F}{N^{\mathrm{off}}_{J=2}}$",
    "eq_ft": r"$f_F = (f_t,\ 1-f_t,\ 0,\ 0),\qquad f_t = \frac{\eta^{\mathrm{off}}_{B,s} - r\,\eta^{\mathrm{off}}_{A,s}}"
             r"{r\,(\eta^{\mathrm{off}}_{A,t}-\eta^{\mathrm{off}}_{A,s}) - (\eta^{\mathrm{off}}_{B,t}-\eta^{\mathrm{off}}_{B,s})}"
             r",\qquad r = \frac{S^{\mathrm{off}}_B}{S^{\mathrm{off}}_A}$",
    "eq_etaoff": r"$\eta^{\mathrm{off}}_{lF} = \frac{1}{\mathrm{\#}L_F}\sum_{m_F\in L_F}\ \sum_{v_x} G(v_x)\, \eta_l(F, m_F; v_x)$",
    "eq_etaon": r"$\eta^{\mathrm{on}}_{lF} = \frac{\sum_{i\in L} w_i \sum_{v_x, y_0} W(v_x, y_0) \sum_{m_F\in F} "
                r"P_{i;\,F,m_F}(v_x, y_0)\ \eta_l(F, m_F; v_x)}{\sum_{i\in L} w_i \sum_{v_x, y_0} W(v_x, y_0) \sum_{m_F\in F} "
                r"P_{i;\,F,m_F}(v_x, y_0)}$",
    # v2 (one definition per symbol): x = off, on; f_F defined once; f_t through r, without the t/s labels
    "eq_signal_v2": r"$S^x_l = \sum_F \eta^x_{lF}\, N^x_F ,\qquad x = \mathrm{off},\ \mathrm{on}$",
    "eq_ft_v2": r"$f_F \equiv \frac{N^{\mathrm{off}}_F}{N^{\mathrm{off}}_{J=2}} = (f_t,\ 1-f_t,\ 0,\ 0),\qquad "
                r"r \equiv \frac{S^{\mathrm{off}}_B}{S^{\mathrm{off}}_A} = \frac{\sum_F \eta^{\mathrm{off}}_{BF}\, f_F}"
                r"{\sum_F \eta^{\mathrm{off}}_{AF}\, f_F}\ \Rightarrow\ f_t = 0.44$",
    "eq_ratio_v2": r"$R_l \equiv \frac{S^{\mathrm{on}}_l}{S^{\mathrm{off}}_l} = \frac{\sum_F \eta^{\mathrm{on}}_{lF}\, n_F}"
                   r"{\sum_F \eta^{\mathrm{off}}_{lF}\, f_F}\ \ \mathrm{(measured)},\qquad n_F \equiv \frac{N^{\mathrm{on}}_F}"
                   r"{N^{\mathrm{off}}_{J=2}}\ \ \mathrm{(want\ to\ know)}$",
    # v3 of the yields: one population-weighted mean for x = off, on; only the population p^x differs
    "eq_eta_v3": r"$\eta^x_{lF} = \frac{\sum_{m_F\in F}\,\sum_{v_x} p^x_{F,m_F}(v_x)\; \eta_l(F, m_F; v_x)}"
                 r"{\sum_{m_F\in F}\,\sum_{v_x} p^x_{F,m_F}(v_x)},\qquad x = \mathrm{off},\ \mathrm{on}$",
    "eq_poff_v3": r"$p^{\mathrm{off}}_{F,m_F}(v_x) = \frac{f_F}{\mathrm{\#}L_F}\, G(v_x)\quad "
                  r"(m_F\in L_F;\ 0\ \mathrm{otherwise})$",
    "eq_pon_v3": r"$p^{\mathrm{on}}_{F,m_F}(v_x) = \sum_{i\in L} w_i \sum_{y_0} W(v_x, y_0)\; P_{i;\,F,m_F}(v_x, y_0)$",
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, tex in EQ.items():
        fig = plt.figure(figsize=(0.1, 0.1))
        fig.text(0, 0, tex, fontsize=20, color="#262624")
        fig.savefig(OUT / f"{name}.png", dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.04)
        plt.close(fig)
        print("wrote", name)


if __name__ == "__main__":
    main()
