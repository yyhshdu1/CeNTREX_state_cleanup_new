"""Absolute laser frequencies for the CU2 test (RC, CU2 at 217 V/cm, DET with plates off or on), transition-finder style.

    python scripts/cu2_test_line_frequencies.py            # -> docs/state_cleanup/CU2_test_line_frequencies.md
    python scripts/cu2_test_line_frequencies.py --cs 351730.547   # with the cesium reading of the day (GHz, IR)

Calibration as in CeNTREX-transition-finder (github.com/ograsdijk/CeNTREX-transition-finder, calibration.py and
transition_finder.calibration_offset_from_grid): each named line (branch, J, F1', F') is the E1-strength-weighted
centroid of its components at E = 0, divided by 4 for the IR, and shifted so that R(2) F1'=7/2 F'=3 sits at
275858770.677 MHz IR. A cesium reading different from the calibration one (351730547.000 MHz IR) shifts every line
by (reading - 351730547.000).

Fields: CU2 lines at the measured 217 V/cm use the main pairs of the CU2 model (F'=3: X F=3 mF=2 -> B F'=3 mF'=3;
F'=4: X F=3 mF=0 -> B F'=4 mF'=1), whose detunings delta0, delta1 are the set-point coordinates. The DET lines are
given at E = 0 (plates off) and, with the plates on (154 V/cm along k), as the strongest sigma+ components about k.
Frequencies the lasers must be set to:
  CU2: laser (F'=3 beam) = F'=3 main line + delta0 (UV); F'=4 beam = laser + f_AOM (single-pass UV AOM, +).
  DET: the DET AOM shifts the light by -80 MHz in the UV, so the laser is set 80 MHz UV (20 MHz IR) above the line.
"""
import argparse
import warnings
from pathlib import Path

import numpy as np

warnings.simplefilter("ignore")
from centrex_tlf import hamiltonian, states  # noqa: E402
from centrex_tlf.hamiltonian.matrix_elements_electric_dipole import generate_ED_ME_mixed_state  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
X, B = states.ElectronicState.X, states.ElectronicState.B
CAL = ("R", 2, 3.5, 3)                    # calibration line of the transition finder
CAL_IR_MHZ = 275858770.677
CAL_CS_MHZ = 351730547.000
TWO_PI_MHZ = 2 * np.pi * 1e6
DJ = {"P": -1, "Q": 0, "R": 1}
B_SMALL = np.array([0.0, 0.0, 1e-3])
E_CU2, E_DET = 217.0, 154.0
SET_POINT = (14.5, 209.0)                  # (delta0, f_AOM) in MHz, 217 V/cm, 32 + 43 mW
DET_AOM_UV = -80.0
SP = np.array([1, -1j, 0]) / 2 ** 0.5     # sigma+ about z (centrex: x - iy drives Delta mF = +1)

# measured in the CU1 test (docs/state_cleanup/CU1_handoff_2026-10-22_v2.md), IR wavemeter
MEASURED = {("P", 2, 1.5, 2): (275842.0401e3, "RC (no AOM)"),
            ("R", 1, 2.5, 3): (275855.4647e3, "Det2 laser (DET AOM -80 MHz UV)"),
            ("Q", 1, 0.5, 1): (275845.27034e3, "Det2 laser (DET AOM -80 MHz UV)")}

# CU2-off and CU2-on peaks of the plates-on CU2-test spectra (centrex probe OBE, sigma+ 15 mW, Doppler averaged;
# CU2 at 217 V/cm, (14.5, 209) MHz, triplet-only lens start; figures CU2_test_detplates154_*_cu2_217Vcm.png).
# Offsets in UV MHz from the line's zero-field position; (offset, photons per J=2 molecule, main contributors).
OBE_PEAKS = {
    ("R", 2, 3.5, 4): {"off": [(-3.31, 0.666, "F=3 m=0, -1"), (+11.69, 0.453, "F=3 m=0"), (+18.19, 0.443, "F=3 m=+1")],
                       "on": [(-2.81, 2.42e-4, "F1=3/2"), (-9.31, 2.27e-4, "F1=3/2 F=1 m=0"),
                              (-16.81, 2.09e-4, "F1=3/2")]},
    ("Q", 2, 2.5, 3): {"off": [(+2.13, 0.897, "F=3 m=-1, 0"), (+18.13, 0.623, "F=3 m=+1"), (-34.87, 0.442, "F=3 m=+1"),
                               (-22.37, 0.437, "F=3 m=0")],
                       "on": [(+5.13, 4.96e-4, "F1=3/2"), (-1.87, 3.43e-4, "F1=3/2"), (+11.63, 3.21e-4, "F1=3/2")]},
    ("Q", 2, 1.5, 2): {"off": [(+1.31, 0.281, "F=3 m=-1"), (-83.69, 0.196, "F=3 m=+1"), (+96.81, 0.182, "F=3 m=+1"),
                               (-37.69, 0.036, "F=3 m=0"), (+53.31, 0.012, "F=3 m=0")],
                       "on": [(-41.19, 5.07e-4, "F1=3/2 m=0"), (-36.69, 4.96e-4, "F1=3/2 m=0"),
                              (+49.81, 4.23e-4, "F1=3/2 m=0"), (+55.31, 4.01e-4, "F1=3/2 m=0"),
                              (+0.81, 3.95e-4, "F1=3/2 m=-1"), (-84.19, 1.70e-4, "F1=3/2 F=2 m=+1")]}}


def name(line):
    br, j, f1, f = line
    return f"{br}({j}) F1'={int(2 * f1)}/2 F'={f}"


def cbs(J, F1, F, mF, P, el):
    return states.CoupledBasisState(J=J, F1=F1, F=F, mF=mF, I1=1 / 2, I2=1 / 2, Ω=0 if el == X else 1, P=P,
                                    electronic_state=el)


def zero_field_levels():
    """X (J 0-5) and B (J' 1-5, both parities) levels at E = 0, keyed (J, F1, F, P) -> energy (MHz, UV)."""
    xs = list(states.generate_coupled_states_X(states.QuantumSelector(J=list(range(0, 5)), electronic=X)))
    bs = list(states.generate_coupled_states_B(states.QuantumSelector(J=list(range(1, 5)), electronic=B,
                                                                      P=[-1, 1], Ω=1)))
    hx = hamiltonian.generate_reduced_X_hamiltonian(xs, E=np.zeros(3), B=B_SMALL, Jmin=0, Jmax=6)
    hb = hamiltonian.generate_reduced_B_hamiltonian(bs, E=np.zeros(3), B=B_SMALL, Jmin=1, Jmax=6)
    out = {}
    for h, el in ((hx, X), (hb, B)):
        for i, q in enumerate(h.QN_basis):
            s = q.largest
            out.setdefault((el, int(s.J), float(s.F1), int(s.F), int(s.P)), []).append(np.real(h.H[i, i]) / TWO_PI_MHZ)
    return {k: float(np.mean(v)) for k, v in out.items()}


def strength(jg, f1g, fg, pg, je, f1e, fe):
    """Sum over mF, mF' and x, y, z of |<e|d|g>|^2 at zero field (pure basis states)."""
    tot = 0.0
    for mg in range(-fg, fg + 1):
        g = 1 * cbs(jg, f1g, fg, mg, pg, X)
        for me in range(-fe, fe + 1):
            if abs(me - mg) > 1:
                continue
            e = 1 * cbs(je, f1e, fe, me, -pg, B)
            tot += sum(abs(generate_ED_ME_mixed_state(e, g, pol_vec=p)) ** 2 for p in np.eye(3))
    return tot


def centroid_uv(line, lev):
    """E1-weighted centroid (MHz, UV, model origin) of a named line over its ground hyperfine levels."""
    br, j, f1e, fe = line
    je, pg = j + DJ[br], (-1) ** j
    f1gs = [abs(j - 0.5), j + 0.5] if j > 0 else [0.5]
    comps = []
    for f1g in f1gs:
        for fg in (int(f1g - 0.5), int(f1g + 0.5)):
            s = strength(j, f1g, fg, pg, je, f1e, fe)
            if s > 1e-9:
                comps.append((lev[(B, je, f1e, fe, -pg)] - lev[(X, j, f1g, fg, pg)], s, f1g, fg))
    f = np.array([c[0] for c in comps]); w = np.array([c[1] for c in comps])
    return float(np.average(f, weights=w)), comps


def main_pair_uv(g, e, E):
    hx = hamiltonian.generate_reduced_X_hamiltonian([g], E=np.array([0, 0, E]), B=B_SMALL, Jmin=0, Jmax=5)
    hb = hamiltonian.generate_reduced_B_hamiltonian([e], E=np.array([0, 0, E]), B=B_SMALL, Jmin=1, Jmax=5)
    return (np.real(hb.H[0, 0]) - np.real(hx.H[0, 0])) / TWO_PI_MHZ


def plates_on_components(line, E, top=6):
    """Strongest sigma+ components (about the field along k) of a J=2 line at E, offsets from the E=0 line (UV MHz)."""
    br, j, f1e, fe = line
    je, pg = j + DJ[br], (-1) ** j
    xs = list(states.generate_coupled_states_X(states.QuantumSelector(J=j, electronic=X)))
    bs = list(states.generate_coupled_states_B(states.QuantumSelector(J=[je], electronic=B, P=[-1, 1], Ω=1)))
    res = {}
    for field in (0.0, E):
        hx = hamiltonian.generate_reduced_X_hamiltonian(xs, E=np.array([0, 0, field]), B=B_SMALL)
        hb = hamiltonian.generate_reduced_B_hamiltonian(bs, E=np.array([0, 0, field]), B=B_SMALL)
        res[field] = (hx, hb)
    hx0, hb0 = res[0.0]
    lab0x = [q.largest for q in hx0.QN_basis]; lab0b = [q.largest for q in hb0.QN_basis]
    ref = np.mean([np.real(hb0.H[b, b]) - np.real(hx0.H[a, a])
                   for a, sx in enumerate(lab0x) for b, sb in enumerate(lab0b)
                   if sb.F1 == f1e and sb.F == fe and sb.P == -pg and sx.mF == sb.mF]) / TWO_PI_MHZ
    hx, hb = res[E]
    comps = []
    for a, qx in enumerate(hx.QN_basis):
        for b, qb in enumerate(hb.QN_basis):
            sb = qb.largest
            if not (sb.F1 == f1e and sb.F == fe):          # both Lambda-doublet partners: parity mixes in the field
                continue
            s = abs(generate_ED_ME_mixed_state(qb, qx, pol_vec=SP)) ** 2
            if s > 1e-4:
                f = (np.real(hb.H[b, b]) - np.real(hx.H[a, a])) / TWO_PI_MHZ - ref
                sx = qx.largest
                comps.append((f, s, f"F1={sx.F1:g} F={sx.F} m={sx.mF:+d} -> m'={sb.mF:+d}, P'={sb.P:+d}"))
    smax = max(c[1] for c in comps)
    return sorted(comps, key=lambda c: -c[1])[:top], smax


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cs", type=float, default=CAL_CS_MHZ / 1e3, help="cesium reading of the day, GHz (IR)")
    args = ap.parse_args(argv)
    cs_shift = args.cs * 1e3 - CAL_CS_MHZ
    lev = zero_field_levels()
    cal_uv, _ = centroid_uv(CAL, lev)

    def ir(line):                                          # zero-field line, IR MHz, calibrated
        return (centroid_uv(line, lev)[0] - cal_uv) / 4 + CAL_IR_MHZ + cs_shift

    rows = []
    out = ["# CU2 test: laser frequencies (transition-finder calibration)", "",
           f"Generated by `scripts/cu2_test_line_frequencies.py`. Calibration: R(2) F1'=7/2 F'=3 = "
           f"{CAL_IR_MHZ / 1e3:.6f} GHz IR at cesium {CAL_CS_MHZ / 1e3:.6f} GHz"
           + (f"; cesium reading used {args.cs:.6f} GHz (shift {cs_shift:+.3f} MHz IR)." if cs_shift else
              " (no cesium correction applied)."),
           "All IR frequencies in GHz, wavemeter units (UV = 4 x IR). Offsets in MHz.", ""]

    # 1. check against the CU1-test measurements
    out += ["## 1. Check against lines measured in the CU1 test", "",
            "| line | model (E = 0), IR GHz | measured, IR GHz | measured - model (IR MHz) | what was measured |",
            "|---|---|---|---|---|"]
    for line, (meas, what) in MEASURED.items():
        f = ir(line)
        out.append(f"| {name(line)} | {f / 1e3:.5f} | {meas / 1e3:.5f} | {meas - f:+.1f} | {what} |")
    out += ["", "Det2 lines: the -80 MHz UV DET AOM alone puts the laser +20.0 MHz IR above the line, so the "
                "residual is (measured - model) - 20 MHz: RC and Q(1) agree to about 1 MHz IR; R(1) F'=3 is "
                "about 11 MHz IR off, unexplained (that day had a cavity-1 mode hop). The transition finder "
                "also applies the cesium reading of the day (--cs).", ""]

    # 2. RC
    out += ["## 2. RC (E = 0, no AOM)", "", "| line | line, IR GHz | laser, IR GHz |", "|---|---|---|"]
    for line in (("P", 2, 1.5, 2), ("P", 2, 1.5, 1)):
        f = ir(line)
        out.append(f"| {name(line)} | {f / 1e3:.5f} | {f / 1e3:.5f} |")
    out += ["", "The measured RC P(2) F'=2 line (275842.0401 GHz) is the better number to use.", ""]

    # 3. CU2 at 217 V/cm
    d0, faom = SET_POINT
    g3, e3 = cbs(2, 2.5, 3, 2, 1, X), cbs(3, 3.5, 3, 3, -1, B)
    g4, e4 = cbs(2, 2.5, 3, 0, 1, X), cbs(3, 3.5, 4, 1, -1, B)
    sh3 = main_pair_uv(g3, e3, E_CU2) - main_pair_uv(g3, e3, 0.0)
    sh4 = main_pair_uv(g4, e4, E_CU2) - main_pair_uv(g4, e4, 0.0)
    z3, z4 = ir(("R", 2, 3.5, 3)), ir(("R", 2, 3.5, 4))
    m3, m4 = z3 + sh3 / 4, z4 + sh4 / 4                     # main lines at 217 V/cm, IR
    laser = m3 + d0 / 4
    beam4_uv = 4 * laser + faom
    d1 = beam4_uv - 4 * m4
    out += [f"## 3. CU2 at {E_CU2:g} V/cm, set point delta0 = {d0:g} MHz, f_AOM = {faom:g} MHz", "",
            "| quantity | IR GHz | UV GHz | note |", "|---|---|---|---|",
            f"| R(2) F1'=7/2 F'=3, E = 0 | {z3 / 1e3:.5f} | {4 * z3 / 1e3:.4f} | calibration line |",
            f"| R(2) F1'=7/2 F'=4, E = 0 | {z4 / 1e3:.5f} | {4 * z4 / 1e3:.4f} | F'=4 - F'=3 = {4 * (z4 - z3):.3f} MHz UV |",
            f"| F'=3 main pair at {E_CU2:g} V/cm | {m3 / 1e3:.5f} | {4 * m3 / 1e3:.4f} | Stark shift {sh3:+.3f} MHz UV |",
            f"| F'=4 main pair at {E_CU2:g} V/cm | {m4 / 1e3:.5f} | {4 * m4 / 1e3:.4f} | Stark shift {sh4:+.3f} MHz UV; "
            f"Delta_main = {4 * (m4 - m3):.3f} MHz UV |",
            f"| **CU2 laser (F'=3 beam)** | **{laser / 1e3:.5f}** | {4 * laser / 1e3:.4f} | delta0 = +{d0:g} MHz UV "
            f"= +{d0 / 4:.3f} MHz IR above the F'=3 main pair |",
            f"| F'=4 beam (after the AOM) | -- | {beam4_uv / 1e3:.4f} | laser + {faom:g} MHz UV; delta1 = {d1:+.2f} MHz |",
            "", f"Laser relative to the zero-field R(2) F'=3 line: {(laser - z3):+.3f} MHz IR "
                f"({4 * (laser - z3):+.2f} MHz UV).", ""]

    # 4. DET
    det = [("R", 2, 3.5, 4), ("Q", 2, 2.5, 3), ("Q", 2, 1.5, 2), ("Q", 2, 1.5, 1), ("Q", 2, 2.5, 2),
           ("R", 2, 3.5, 3), ("R", 2, 2.5, 3), ("R", 2, 2.5, 2), ("R", 1, 2.5, 3), ("R", 0, 1.5, 2)]
    out += [f"## 4. DET (AOM {DET_AOM_UV:+g} MHz UV: laser = line + {-DET_AOM_UV:g} MHz UV = line + "
            f"{-DET_AOM_UV / 4:g} MHz IR)", "", "### Plates off (E = 0)", "",
            "| line | line, IR GHz | **DET laser, IR GHz** |", "|---|---|---|"]
    for line in det:
        f = ir(line)
        out.append(f"| {name(line)} | {f / 1e3:.5f} | **{(f - DET_AOM_UV / 4) / 1e3:.5f}** |")
    out += ["", f"### Plates on ({E_DET:g} V/cm along k), strongest sigma+ components about k", "",
            "Offsets are from the line's zero-field position, UV MHz (molecule frame); the DET laser is at "
            "(zero-field DET laser above) + offset/4 in the IR. Strengths relative to the strongest component.", ""]
    for line in det[:3]:
        comps, smax = plates_on_components(line, E_DET, top=8)
        f = ir(line)
        out += [f"**{name(line)}** (E = 0 laser {(f - DET_AOM_UV / 4) / 1e3:.5f} GHz IR)", "",
                "| component (ground about k -> excited) | offset, UV MHz | DET laser, IR GHz | strength |",
                "|---|---|---|---|"]
        for fo, s, lab in comps:
            out.append(f"| {lab} | {fo:+.2f} | {(f - DET_AOM_UV / 4 + fo / 4) / 1e3:.5f} | {s / smax:.2f} |")
        out.append("")
        if line in OBE_PEAKS:
            out += ["Expected CU2-test signal peaks (OBE, Doppler and power broadened):", "",
                    "| CU2 | offset, UV MHz | **DET laser, IR GHz** | photons per J=2 molecule | from |",
                    "|---|---|---|---|---|"]
            for k in ("off", "on"):
                for fo, ph, frm in OBE_PEAKS[line][k]:
                    out.append(f"| {k} | {fo:+.2f} | **{(f - DET_AOM_UV / 4 + fo / 4) / 1e3:.5f}** | {ph:.3g} | {frm} |")
            out.append("")
    path = ROOT / "docs" / "state_cleanup" / "CU2_test_line_frequencies.md"
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    print("\n".join(out))
    print("wrote", path)


if __name__ == "__main__":
    main()
