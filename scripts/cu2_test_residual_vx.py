"""Section 7.4 of cu2_residual_from_line_combinations_2026-11-09: eta_on with the vx distribution of the CU2 residual,
from the CU2 OBE at the CU2-test setting, then the line-combination solve again.

    python -I scripts/cu2_test_residual_vx.py [power_mW]

CU2 OBE: transfer stage (scripts/cu2_power_check.py transfer 13,209 20 30 e217 acc), i.e. NLTL2 2068 = delta0 13 MHz,
f_AOM 209 MHz, 20 + 30 mW, 217 V/cm, 5 % loss per pass, v_z 156 m/s; final[start, vx, y0, level] on vx = +-0.2,
+-0.6, +-1.0 m/s and y0 = 0, 3, 6, 9 mm. Weights W[vx, y0] = acceptance_weights (Gaussian vx, sigma 0.65 m/s,
|vx| <= 1 m/s; y0 Gaussian sigma 10.35 mm, |y0| <= 8.5 mm) -- the same G(vx) x y0 acceptance as the document.

  eta_on[l, g] = sum_s w_s sum_{vx,y0} W  sum_{m in g} P_fin[s; g, m](vx, y0) eta_l(g, m; vx)
                 / sum_s w_s sum_{vx,y0} W  sum_{m in g} P_fin[s; g, m](vx, y0)
  start weights w_s: f_t over F=3 mF = 0, +-1, 1 - f_t on F=2 mF = 0 (mF < 0 starts mirror the mF > 0 ones).
eta_off: lens m_F, same vx weights. Solve sum_g eta_on n_g = R_l sum_g eta_off f_g, f = (f_t, 1 - f_t, 0, 0), f_t from
S_B/S_A; Monte Carlo over R_l and S_B/S_A (20000).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import cu2_test_start_fit as SF  # noqa: E402
import plot_cu2_test_vs_sim as P  # noqa: E402

T = SF.T
G = list(SF.GROUPS)
NAME = SF.GROUP_NAME
LINES = ["F4", "R7o2F3", "F1p5F3", "F1p5F2"]                    # A, B, C, D
R = np.array([0.0019, 0.0057, 0.026, 0.013]); RERR = np.array([0.0059, 0.0037, 0.014, 0.019])
RBA, RBA_ERR = 1.56, 0.17
TRANSFER = ROOT / "tmp" / "cu2_power" / "transfer_loss5_d13_f209_20_30_e217_acc_result.npz"


def main(power_mW=24.0):
    z = np.load(TRANSFER, allow_pickle=True)
    lab = [(float(a), int(b), int(c)) for a, b, c in zip(z["F1"], z["F"], z["mF"])]
    st = [(float(a), int(b), int(c)) for a, b, c in zip(z["start_F1"], z["start_F"], z["start_mF"])]
    vx, y0 = np.asarray(z["vx"], float), np.asarray(z["y0"], float)
    W = P.acceptance_weights(vx, y0, svx=0.65)                       # [vx, y0], sums to 1
    wv = W.sum(1)                                                    # G(vx) on the vx nodes
    fin = np.asarray(z["final"])                                     # [start, vx, y0, level]
    mirror = [lab.index((l[0], l[1], -l[2])) for l in lab]

    def final_of(start):                                             # [vx, y0, level] for any start level
        if start[2] >= 0:
            return fin[st.index(start)]
        return fin[st.index((start[0], start[1], -start[2]))][..., mirror]
    T.POL = "k"
    Y = {k: T.level_photons(k, lab, vx, power=1e-3 * power_mW) for k in LINES}    # [vx, level]
    gi = {g: [i for i, l in enumerate(lab) if (l[0], l[1]) == g] for g in G}

    def eta_off():
        return np.array([[np.mean([wv @ Y[k][:, lab.index((g[0], g[1], m))] for m in SF.GROUPS[g]]) for g in G]
                         for k in LINES])

    def eta_on(ft, residual_vx=True):
        num = np.zeros((len(LINES), len(G))); den = np.zeros(len(G))
        for g0, wg in (((2.5, 3), ft), ((2.5, 2), 1 - ft)):
            for m in SF.GROUPS[g0]:
                Pf = final_of((g0[0], g0[1], m)) * wg / len(SF.GROUPS[g0])            # [vx, y0, level]
                if residual_vx:
                    pv = np.einsum("vy,vyl->vl", W, Pf)                              # [vx, level], residual-weighted
                else:                                                                # lens vx: same total, G(vx)
                    pv = np.outer(wv, np.einsum("vy,vyl->l", W, Pf))
                for j, g in enumerate(G):
                    den[j] += pv[:, gi[g]].sum()
                    for i, k in enumerate(LINES):
                        num[i, j] += (pv[:, gi[g]] * Y[k][:, gi[g]]).sum()
        return num / den, den
    eoff = eta_off()
    a0, a1, b0, b1 = eoff[0, 0], eoff[0, 1], eoff[1, 0], eoff[1, 1]
    ft_of = lambda r: (b1 - r * a1) / (r * (a0 - a1) - (b0 - b1))                  # noqa: E731
    ft0 = ft_of(RBA)
    eon_res, n_mod = eta_on(ft0, True)
    eon_lens, _ = eta_on(ft0, False)
    print(f"CU2 OBE 217 V/cm, (13, 209) MHz, 20 + 30 mW, v_z {T.VZ:g} m/s; probe {power_mW:g} mW; f_t = {ft0:.3f}")
    print("\neta_on(residual vx) / eta_on(lens vx):")
    for k, tag, row_r, row_l in zip(LINES, "ABCD", eon_res, eon_lens):
        print(f"  {tag} {k:7s} " + "  ".join(f"{NAME[g]} {a / b:.2f}" if b > 0.05 else f"{NAME[g]}   - "
                                             for g, a, b in zip(G, row_r, row_l)))
    print("  (document 7.3, old r(vx): A F=3 0.94, B F=2 0.93, C 3/2F=2 0.94, D 3/2F=1 0.97; "
          "B F=3 0.79, C F=3 0.62, D 5/2F=2 0.63, D 3/2F=2 0.79)")
    print("\nmodel residual (% of J=2):", "  ".join(f"{NAME[g]} {100 * n:.3f}" for g, n in zip(G, n_mod)),
          f"  total {100 * n_mod.sum():.3f}")
    rng = np.random.default_rng(1)
    for label, res_vx in (("lens vx (no residual vx)", False), ("residual vx from the CU2 OBE", True)):
        def solve(Rv, r):
            ft = ft_of(r)
            eon, _ = eta_on(ft, res_vx)
            return np.linalg.solve(eon, Rv * (eoff @ np.array([ft, 1 - ft, 0, 0])))
        n0 = solve(R, RBA)
        S = np.array([solve(R + RERR * rng.standard_normal(4), RBA + RBA_ERR * rng.standard_normal())
                      for _ in range(4000)])
        e = S.std(0)
        tot = S.sum(1).std()
        print(f"\n=== solve with eta_on from {label}")
        for g, n, s in zip(G, n0, e):
            print(f"  {NAME[g]:11s} {100 * n:+.3f} +- {100 * s:.2f} %")
        print(f"  J=2 total   {100 * n0.sum():.3f} +- {100 * tot:.2f} %  -> cleanup {100 * (1 - n0.sum()):.2f} %; "
              f"F=3 {100 * (1 - n0[0] / ft0):.2f} %, F1=5/2 F=2 {100 * (1 - n0[1] / (1 - ft0)):.2f} %")


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 24.0)
