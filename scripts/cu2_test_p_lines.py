"""Can the P(2) detection lines measure the CU2 cleanup, and what do they add to the four R(2) lines?

    python -I scripts/cu2_test_p_lines.py <CU2_test_processed_results_2026-11-07_08.h5> [power_mW]

DET plates off (E = 0), sigma+ about k, power_mW (24 mW measured), v_z 156 m/s; probe photon yields from
cu2_detection_test.level_photons (vx-averaged, sigma 0.65 m/s). CU2: transfer stage on the acceptance grid at 156 m/s,
start = lens levels with the split of cu2_residual_from_line_combinations_2026-11-09 (F=3 43.7 %, singlet 56.3 %).
For each line:
  - photons per molecule of each J=2 group (CU2 off: lens m_F; CU2 on: m_F of the model residual);
  - CU2-off signal predicted on the scale of the measured R lines (K = S_off / (eta . f) from lines A and B), and the
    share of each group in it, CU2 off and on; the model on/off ratio;
  - expected error of R from the measured noise floor: R_err * S_off ~ const (500-shot runs; lines A-D);
  - error on each solved group when the line is added to the A-D solve (linearised weighted least squares, f fixed).
A line that sees only F1=3/2 has no CU2-off signal: it measures the CU2-on population directly, and the noise floor
gives the 1-sigma sensitivity in % of J=2.
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
F_T = 0.437
R_LINES = {"A": ("F4", "R(2) F1'=7/2 F'=4", 0.0019, 0.0059, 7.09), "B": ("R7o2F3", "R(2) F1'=7/2 F'=3", 0.0057, 0.0037, 12.05),
           "C": ("F1p5F3", "R(2) F1'=5/2 F'=3", 0.026, 0.014, 2.52), "D": ("F1p5F2", "R(2) F1'=5/2 F'=2", 0.013, 0.019, 2.99)}
P_LINES = {"P0": ("P1o2F0", "P(2) F1'=1/2 F'=0"), "P1": ("P1o2F1", "P(2) F1'=1/2 F'=1"),
           "P2": ("P3o2F1", "P(2) F1'=3/2 F'=1"), "P3": ("P3o2F2", "P(2) F1'=3/2 F'=2")}
TRANSFER = ROOT / "tmp" / "cu2_power" / "transfer_loss5_d13_f209_20_30_e217_acc_result.npz"


def model():
    z = np.load(TRANSFER, allow_pickle=True)
    lab = [(float(a), int(b), int(c)) for a, b, c in zip(z["F1"], z["F"], z["mF"])]
    st = [(float(a), int(b), int(c)) for a, b, c in zip(z["start_F1"], z["start_F"], z["start_mF"])]
    fin = np.einsum("vy,svyl->sl", P.acceptance_weights(z["vx"], z["y0"], svx=0.65), z["final"])
    mirror = {l: lab.index((l[0], l[1], -l[2])) for l in lab}
    rho0 = np.zeros(len(lab))
    after = np.zeros(len(lab))
    for g, w in (((2.5, 3), F_T), ((2.5, 2), 1 - F_T)):
        for m in SF.GROUPS[g]:
            i = lab.index((g[0], g[1], m))
            rho0[i] += w / len(SF.GROUPS[g])
            f_ = fin[st.index((g[0], g[1], abs(m)))]
            after += w / len(SF.GROUPS[g]) * (f_ if m >= 0 else np.array([f_[mirror[q]] for q in lab]))
    return lab, rho0, after


def main(h5path, power_mW=24.0):
    lab, rho0, after = model()
    keys = [v[0] for v in R_LINES.values()] + [v[0] for v in P_LINES.values()]
    T.POL = "k"
    w = np.exp(-0.5 * (SF.VXS / SF.SVX) ** 2); w /= w.sum()
    Y = {k: w @ T.level_photons(k, lab, SF.VXS, power=1e-3 * power_mW) for k in keys}
    gi = {g: [i for i, l in enumerate(lab) if (l[0], l[1]) == g] for g in G}
    f = np.array([F_T, 1 - F_T, 0, 0])
    eta_off = {k: np.array([np.mean([Y[k][lab.index((g[0], g[1], m))] for m in SF.GROUPS[g]]) for g in G]) for k in keys}
    eta_on = {k: np.array([Y[k][gi[g]] @ after[gi[g]] / after[gi[g]].sum() for g in G]) for k in keys}
    n_mod = np.array([after[gi[g]].sum() for g in G])
    K = np.mean([R_LINES[x][4] / (eta_off[R_LINES[x][0]] @ f) for x in ("A", "B")])     # mV ms per photon/molecule
    c = np.mean([v[3] * v[4] for v in R_LINES.values()])                                # noise floor, mV ms
    print(f"probe {power_mW:g} mW, v_z {T.VZ:g} m/s; scale K = {K:.2f} mV ms per (photon per J=2 molecule); "
          f"noise floor R_err*S_off = {c:.3f} mV ms (500 shots)\n")
    print(f"{'line':22s} | photons per molecule F1=5/2 F=3, F=2, F1=3/2 F=2, F=1 (CU2 on m_F) | S_off pred | "
          "CU2-off share | CU2-on share | R model | R_err est")
    rows = {}
    for tag, v in list(R_LINES.items()) + list(P_LINES.items()):
        k, name = v[0], v[1]
        S_off = K * (eta_off[k] @ f)
        S_on = K * (eta_on[k] @ n_mod)
        sh_off = eta_off[k] * f / max(eta_off[k] @ f, 1e-12)
        sh_on = eta_on[k] * n_mod / (eta_on[k] @ n_mod)
        R = S_on / S_off if S_off > 1e-6 else np.inf
        Rerr = c / S_off if S_off > 1e-6 else np.inf
        rows[tag] = (k, name, S_off, S_on, R, Rerr)
        fmt = lambda v: " ".join(f"{100 * x:3.0f}" for x in v)                        # noqa: E731
        print(f"{tag:3s} {name:18s} | {' '.join(f'{x:5.2f}' for x in eta_on[k])}  (off: {' '.join(f'{x:5.2f}' for x in eta_off[k])}) "
              f"| {S_off:7.3f} | {fmt(sh_off)} % | {fmt(sh_on)} % | "
              + (f"{100 * R:6.2f} % | {100 * Rerr:5.2f} %" if np.isfinite(R) else
                 f"no CU2-off signal; CU2-on S {S_on:.4f} mV ms = {S_on / c:.2f} x noise floor"))
    # error propagation: rows y_l = eta_on[l] . n, sigma_l = R_err_l * (eta_off[l] . f)
    def cov(lines):
        A, s = [], []
        for tag in lines:
            k = rows[tag][0]
            A.append(eta_on[k])
            if tag in R_LINES:
                s.append(R_LINES[tag][3] * (eta_off[k] @ f))
            elif np.isfinite(rows[tag][5]):
                s.append(rows[tag][5] * (eta_off[k] @ f))
            else:                                    # background-free line: S_on/K known to c/K, i.e. eta_on . n to c/K
                s.append(c / K)
        A, s = np.array(A), np.array(s)
        return np.linalg.inv(A.T @ np.diag(1 / s ** 2) @ A)
    print("\n1-sigma error on each group after CU2 (% of J=2), linearised, start split fixed:")
    print(f"{'lines':22s} " + "  ".join(f"{NAME[g]:>11s}" for g in G) + "   F1=5/2 sum  F1=3/2 sum   total")
    sets = [("A B C D", ["A", "B", "C", "D"])] + [(f"A B C D + {t}", ["A", "B", "C", "D", t]) for t in P_LINES] + \
           [("A B C D + all P", ["A", "B", "C", "D"] + list(P_LINES)), ("A B + all P", ["A", "B"] + list(P_LINES))]
    for name, ls in sets:
        try:
            C = cov(ls)
        except np.linalg.LinAlgError:
            print(f"{name:22s} singular"); continue
        e = np.sqrt(np.diag(C))
        u = lambda v: np.sqrt(v @ C @ v)                                                  # noqa: E731
        print(f"{name:22s} " + "  ".join(f"{100 * x:11.2f}" for x in e) +
              f"   {100 * u(np.array([1, 1, 0, 0])):10.2f}  {100 * u(np.array([0, 0, 1, 1])):10.2f}  "
              f"{100 * u(np.ones(4)):7.2f}")
    print("\nmodel after CU2 (% of J=2):", "  ".join(f"{NAME[g]} {100 * n:.3f}" for g, n in zip(G, n_mod)))


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 24.0)
