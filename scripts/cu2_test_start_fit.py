"""CU2 test 2026-11-08: populations of the four X J=2 groups with CU2 off, fitted to the four measured CU2-off signals.

    python -I scripts/cu2_test_start_fit.py <CU2_test_processed_results_2026-11-07_08.h5> [power_mW] [--weak s ...]

Model:  S_off(line) = sum_g  A[line, g] * N_g

  S_off(line)  measured CU2-off (shutter closed) fluorescence per shot, V ms: HDF5 cu2_onoff/<run>.attrs
               ["fl_closed_per_shot_Vms"], the four lines at NLTL2 2068 (11-08).
  A[line, g]   photons per molecule that starts in group g, from the DET probe OBE (scripts/cu2_detection_test.py,
               level_photons: plates off, E = 0, sigma+ about k, power_mW (24 mW measured), flat top 1.5 x 1.5 cm,
               v_z 156 m/s), averaged over vx (Gaussian, sigma 0.65 m/s, grid -1...1 m/s) and over the m_F the group is
               assumed to occupy:
                  F1=5/2 F=3   m_F = 0, +-1    (lens levels, m_J ~ 0)
                  F1=5/2 F=2   m_F = 0         (lens level, m_J ~ 0)
                  F1=3/2 F=2   all m_F, uniform
                  F1=3/2 F=1   all m_F, uniform
  N_g          molecules in group g (common, unknown scale: same beam flux and light collection on every line);
               printed normalised to sum 1.

Solved two ways: exactly (np.linalg.solve; can go negative) and with N >= 0 (scipy.optimize.nnls, used by the
figures). --weak s multiplies the off-diagonal entries of A (a group seen on another group's line; at E = 0 the 20
J=2 levels lie within 0.36 MHz, so these are the saturated weak couplings) by s, to show how much the fit depends on
them. The HDF5 has no uncertainty on S_off, so the fit is unweighted and gives no errors.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import cu2_detection_test as T  # noqa: E402

LINES = {"R2F4_triplet_nltl2_2068": ("F4", "R(2) F1'=7/2 F'=4"),             # HDF5 run -> (probe key, name)
         "R2F3_singlet_nltl2_2068": ("R7o2F3", "R(2) F1'=7/2 F'=3"),
         "R2F5o2F3_F1eq3o2F2_nltl2_2068": ("F1p5F3", "R(2) F1'=5/2 F'=3"),
         "R2F5o2F2_F1eq3o2F1_nltl2_2068": ("F1p5F2", "R(2) F1'=5/2 F'=2")}
GROUPS = {(2.5, 3): (0, -1, 1), (2.5, 2): (0,), (1.5, 2): (-2, -1, 0, 1, 2), (1.5, 1): (-1, 0, 1)}  # (F1, F): m_F
GROUP_NAME = {(2.5, 3): "F1=5/2 F=3", (2.5, 2): "F1=5/2 F=2", (1.5, 2): "F1=3/2 F=2", (1.5, 1): "F1=3/2 F=1"}
LEVELS = [(F1, F, m) for F1, F in ((1.5, 1), (1.5, 2), (2.5, 2), (2.5, 3)) for m in range(-F, F + 1)]
VXS = np.round(np.arange(-1.0, 1.0001, 0.2), 6)
SVX = 0.65


def measured(h5path):
    """CU2-off and CU2-on signals and on/off ratios of the four lines -> dict run -> dict."""
    import h5py
    out = {}
    with h5py.File(h5path, "r") as f:
        for run in LINES:
            a = f["cu2_onoff"][run].attrs
            out[run] = dict(S_off=float(a["fl_closed_per_shot_Vms"]), S_on=float(a["fl_open_per_shot_Vms"]),
                            R=float(a["R_template"]), R_err=float(a["R_template_err"]))
    return out


def level_yields(power_mW, levels=LEVELS, svx=SVX):
    """Photons per molecule for every J=2 level on every line, vx-averaged -> dict probe key -> array[level]."""
    T.POL = "k"
    w = np.exp(-0.5 * (VXS / svx) ** 2)
    w /= w.sum()
    return {key: w @ T.level_photons(key, list(levels), VXS, power=1e-3 * power_mW) for key, _ in LINES.values()}


def yield_matrix(Y, levels=LEVELS, weak=1.0):
    """A[line, group]: mean photons per molecule over the group's m_F; off-diagonal entries scaled by weak."""
    lv = list(levels)
    A = np.array([[np.mean([Y[key][lv.index((g[0], g[1], m))] for m in ms]) for g, ms in GROUPS.items()]
                  for key, _ in LINES.values()])
    off = ~np.eye(len(GROUPS), dtype=bool)
    A[off] *= weak
    return A


def fit(A, S):
    """-> (exact solution, NNLS solution), both normalised to sum 1, and the NNLS prediction A @ N (unnormalised)."""
    from scipy.optimize import nnls
    exact = np.linalg.solve(A, S)
    N, _ = nnls(A, S)
    return exact / exact.sum(), N / N.sum(), A @ N


def main(h5path, power_mW=24.0, weaks=(1.0,)):
    m = measured(h5path)
    S = np.array([m[r]["S_off"] for r in LINES])
    Y = level_yields(power_mW)
    names = [n for _, n in LINES.values()]
    print(f"CU2 test 2026-11-08, NLTL2 2068; probe {power_mW:g} mW, sigma+ about k, E = 0, v_z {T.VZ:g} m/s, "
          f"vx Gaussian sigma {SVX} m/s\n")
    print("measured CU2-off signal S_off (V ms per shot)")
    for r, n in zip(LINES, names):
        print(f"  {n:18s} {m[r]['S_off']:.5f}")
    for weak in weaks:
        A = yield_matrix(Y, weak=weak)
        exact, nn, pred = fit(A, S)
        print(f"\n=== off-diagonal (weak-coupling) yields x {weak:g}")
        print("A = photons per molecule  (rows: line; columns: " + ", ".join(GROUP_NAME.values()) + ")")
        for n, row in zip(names, A):
            print(f"  {n:18s} " + "  ".join(f"{x:6.3f}" for x in row))
        print("fitted CU2-off populations (normalised)")
        print(f"  {'group':12s} {'exact':>8s} {'N >= 0':>8s}")
        for g, a, b in zip(GROUPS, exact, nn):
            print(f"  {GROUP_NAME[g]:12s} {100 * a:7.1f}% {100 * b:7.1f}%")
        print("CU2-off signal predicted by the N >= 0 fit / measured - 1")
        for n, p, s in zip(names, pred, S):
            print(f"  {n:18s} {p:.5f} vs {s:.5f}  ({100 * (p / s - 1):+.0f} %)")
        Nn = nn * pred[0] / (A[0] @ nn)           # any common scale; shares do not depend on it
        print("share of each line's CU2-off signal per group (N >= 0 fit)")
        for n, row in zip(names, A):
            sh = row * Nn / (row @ Nn)
            print(f"  {n:18s} " + ", ".join(f"{GROUP_NAME[g]} {100 * x:.0f} %" for g, x in zip(GROUPS, sh) if x > 0.005))


if __name__ == "__main__":
    args = sys.argv[1:]
    weaks = [float(args[i + 1]) for i, a in enumerate(args) if a == "--weak"] or [1.0]
    rest = [a for i, a in enumerate(args) if a != "--weak" and (i == 0 or args[i - 1] != "--weak")]
    main(rest[0], float(rest[1]) if len(rest) > 1 else 24.0, weaks)
