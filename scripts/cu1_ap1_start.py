"""CU1 starting from the J=1 states that RC + SPA AP1 actually populate, with every J=1 sublevel tracked.

    python scripts/cu1_ap1_start.py compute --power 80 --height 3    # 2 starts x 21 x 11 (vx, y0) at +19.0 MHz
    python scripts/cu1_ap1_start.py analyze --power 80 --height 3    # efficiencies, composition, detection readings

Initial state (zero-field labels, adiabatic through SPA -> zero field -> CU1 field; these are the four mainly-m_J=0
levels of the upper group at 150 V/cm):
    triplet  J=0, F=1 -> J=1, F1=3/2, F=2, mF = -1, 0, +1 (equal weights)
    singlet  J=0, F=0 -> J=1, F1=3/2, F=1, mF = 0
RC with P(2) F1'=3/2 F'=1: the excited level decays to J=0, F=1 with 0.147 and to J=0, F=0 with 0.337, so the
default weights are triplet : singlet = 0.147 : 0.337 (AP1 transfers all four J=0 sublevels alike).
The two starts are run separately (the OBE is linear in rho0), so any weights can be combined in analyze.
Model as scripts/cu1_conditions_scan.py: build_cu1, +19.0 MHz, no EOM, Doppler +-k vx alternating pass by pass;
final populations of the 12 J=1 sublevels are kept. Means: acceptance weights of cu1_3cm_scan.py.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.cleanup_spb_doppler import CU1_SOLVER, build_cu1, describe  # noqa: E402
from scripts.cu1_3cm_scan import SY_BEAM, VX_F, Y_F  # noqa: E402
from scripts.cu1_conditions_scan import cond_tag  # noqa: E402
from scripts.cu2_laser_aom_scan import T_END, log  # noqa: E402
from scripts.plot_cleanup_maps_current import acceptance_mean  # noqa: E402

HDF = ROOT / "results" / "SPB_J1_cleanup_scans.h5"
DETUNING = 19.0
STARTS = {"triplet": [(1.5, 2, -1), (1.5, 2, 0), (1.5, 2, 1)], "singlet": [(1.5, 1, 0)]}
W_TRIPLET, W_SINGLET = 0.147, 0.337
CHUNK, THREADS = 60, 20
SLOTS = ["laser_detuning", "vx", "y0", "laser_power"]


def labels_of(obe, gidx):
    out = []
    for i in gidx:
        d = describe(obe.QN[i])
        out.append((float(d["F1"]), int(d["F"]), int(round(d["mF"]))))
    return out


def run_name(start, ct):
    return f"doppler_cu1_ap1_{start}_{ct}_det19MHz"


def compute(power_mw, height_cm):
    from centrex_tlf import lindblad
    ct = cond_tag(power_mw, height_cm)
    out = ROOT / "tmp" / "cu1_ap1" / ct
    out.mkdir(parents=True, exist_ok=True)
    prepared, gidx, obe = build_cu1(return_obe=True, doppler=True, sigma_y_m=height_cm / 100 / 4)
    labs = labels_of(obe, gidx)
    pts = np.array([(v, y) for v in VX_F for y in Y_F])
    for start, keys in STARTS.items():
        idx = [gidx[labs.index(k)] for k in keys]
        for i in idx:                                   # the AP1 levels are the mainly-m_J=0 upper group
            assert describe(obe.QN[i])["p_mj"].get(0, 0) > 0.99
        rho0 = np.zeros((len(obe.QN),) * 2, dtype=np.complex128)
        rho0[idx, idx] = 1.0 / len(idx)
        f = out / f"{start}.npz"
        if f.exists():
            log(f"{start}: already computed ({f})")
            continue
        pops = []
        for s in range(0, len(pts), CHUNK):
            chunk = pts[s:s + CHUNK]
            batch = np.column_stack([np.full(len(chunk), 2 * np.pi * 1e6 * DETUNING), chunk[:, 0], chunk[:, 1],
                                     np.full(len(chunk), power_mw * 1e-3)]).astype(np.complex128)
            t0 = time.perf_counter()
            res = lindblad.parameter_scan(prepared, rho0, (0.0, T_END), parameter_slots=SLOTS, parameter_batch=batch,
                                          threads=THREADS, **CU1_SOLVER)
            pops.append(np.asarray(res.values).real[:, gidx])
            log(f"{start} {ct}: {s + len(chunk)}/{len(pts)} points, chunk {time.perf_counter() - t0:.0f} s")
        pops = np.concatenate(pops).reshape(len(VX_F), len(Y_F), len(gidx))
        np.savez(f, pops=pops, labels=np.array(labs), vx=VX_F, y=Y_F)
        save_hdf(start, ct, power_mw, height_cm, pops, labs)


def save_hdf(start, ct, power_mw, height_cm, pops, labs):
    from scripts.scan_hdf import save_scan_hdf
    name = run_name(start, ct)
    try:
        save_scan_hdf(HDF, name, pops,
                      axes={"vx_m_per_s": (VX_F, "m/s"), "y0_mm": (Y_F * 1e3, "mm"),
                            "level_index": (np.arange(len(labs)), "")},
                      quantity="final_population_of_each_X_J1_sublevel",
                      fixed_parameters={"laser_detuning_MHz": DETUNING, "laser_power_mW": power_mw,
                                        "beam_height_4sigma_cm": height_cm, "eom": "none",
                                        "initial_state": f"AP1 {start}: " + json.dumps(STARTS[start]),
                                        "level_labels_F1_F_mF": json.dumps(labs),
                                        "doppler": "k*vx, sign alternating pass by pass",
                                        "y_symmetry": "results for y0 >= 0; symmetric in y0"},
                      model={"script": "scripts/cu1_ap1_start.py", "build": "cleanup_spb_doppler.build_cu1"},
                      solver=dict(CU1_SOLVER), notes="CU1 from the AP1-populated J=1 levels; all J=1 sublevels kept.",
                      overwrite=False)
        log(f"saved {name}")
    except Exception as exc:
        log(f"could not save {name}: {exc}")


def load(ct):
    out = ROOT / "tmp" / "cu1_ap1" / ct
    d = {s: np.load(out / f"{s}.npz") for s in STARTS}
    labs = [tuple(x) for x in d["triplet"]["labels"].tolist()]
    labs = [(float(a), int(b), int(c)) for a, b, c in labs]
    return {s: d[s]["pops"] for s in STARTS}, labs


def probe_yields():
    """Saturated photons per molecule (50 mW/cm2, centre, vx = 0) from scripts/detection_level_diagrams.py."""
    data = json.loads((ROOT / "tmp" / "detection" / "probe_yields.json").read_text())

    def at(run):
        r = data["runs"][run]
        k = [tuple(c) for c in r["combos"]].index((500.0, 0.0, 0.0))
        return {tuple((float(a), int(b), int(c))): v for (a, b, c), v in zip(r["labels"], r["photons"][k])}
    return at("R1F3_par"), at("Q1F1_par"), at("Q1F1_perp")


def analyze(power_mw, height_cm, w_t=W_TRIPLET, w_s=W_SINGLET):
    ct = cond_tag(power_mw, height_cm)
    pops, labs = load(ct)
    w = {"triplet": w_t / (w_t + w_s), "singlet": w_s / (w_t + w_s)}
    i0 = int(np.argmin(abs(VX_F)))

    def acc(grid):                                       # mean in acceptance of a (vx, y0) grid
        return acceptance_mean(np.maximum(grid, 1e-30), VX_F, Y_F, SY_BEAM)    # log-interpolated in y0

    # acceptance-weighted mean population of each sublevel: mean of each level map (linear, so this is exact)
    per = {s: np.array([acc(pops[s][:, :, k]) for k in range(len(labs))]) for s in STARTS}
    axis = {s: pops[s][i0, 0] for s in STARTS}
    mix = w["triplet"] * per["triplet"] + w["singlet"] * per["singlet"]
    mix_axis = w["triplet"] * axis["triplet"] + w["singlet"] * axis["singlet"]
    print(f"CU1 {power_mw:g} mW, {height_cm:g} cm, +19.0 MHz; weights triplet {w['triplet']:.3f} "
          f"(each sublevel {w['triplet'] / 3:.3f}), singlet {w['singlet']:.3f}")
    for s in STARTS:
        print(f"  start {s:8s}: J=1 left on axis {100 * axis[s].sum():.3f} %, mean in acceptance {100 * per[s].sum():.3f} %"
              f"  -> efficiency {100 * (1 - per[s].sum()):.2f} %")
    print(f"  mixture      : J=1 left on axis {100 * mix_axis.sum():.3f} %, mean in acceptance {100 * mix.sum():.3f} %"
          f"  -> efficiency {100 * (1 - mix.sum()):.2f} %")
    print("  where the remaining J=1 sits (mean in acceptance, % of all molecules):")
    order = sorted(range(len(labs)), key=lambda k: (-labs[k][0], -labs[k][1], labs[k][2]))
    for k in order:
        f1, f, mf = labs[k]
        print(f"     F1={f1:g} F={f} mF={mf:+d}: triplet start {100 * per['triplet'][k]:.4f}  "
              f"singlet start {100 * per['singlet'][k]:.4f}  mixture {100 * mix[k]:.4f}")
    yR, yP, yX = probe_yields()
    init = np.array([w["triplet"] / 3 if (l[0], l[1]) == (1.5, 2) and abs(l[2]) <= 1 else
                     (w["singlet"] if l == (1.5, 1, 0) else 0.0) for l in labs])
    groups = {}
    for k, l in enumerate(labs):
        groups.setdefault(l[:2], []).append(k)
    yI = {l: np.mean([yP[labs[j]] for j in groups[l[:2]]]) for l in labs}
    print("  detection (saturated probe; ratio S_on / S_off of the AP1-difference signal):")
    for name, y in (("R(1) F'=3", yR), ("Q(1) F'=1, pol par", yP), ("Q(1) F'=1, pol perp", yX),
                    ("Q(1) F'=1, scrambled", yI)):
        yy = np.array([y[l] for l in labs])
        print(f"     {name:22s}: {100 * (mix @ yy) / (init @ yy):.3f} %")
    nF2_on = sum(mix[k] for k, l in enumerate(labs) if l[:2] == (1.5, 2))
    nF2_off = w["triplet"]
    sQ_par = mix @ np.array([yP[l] for l in labs])
    bound = (nF2_on + sQ_par / yP[(0.5, 0, 0)]) / nF2_off
    trip_f2 = w["triplet"] * sum(per["triplet"][k] for k, l in enumerate(labs) if l[:2] == (1.5, 2))
    sing_f2 = w["singlet"] * sum(per["singlet"][k] for k, l in enumerate(labs) if l[:2] == (1.5, 2))
    print(f"     R(1) F'=3 numerator: from triplet {100 * trip_f2 / nF2_off:.3f} %, from singlet pumped into F=2 "
          f"{100 * sing_f2 / nF2_off:.3f} % (of the CU1-off F=2 signal)")
    print(f"     bound (S_R,on/2.1 + S_Q,on/4.4)/(S_R,off/2.1), pol par: {100 * bound:.3f} %  "
          f"(true J=1 left / CU1-off J=1: {100 * mix.sum():.3f} %; the bound normalizes to F=2 only)")
    return dict(per=per, mix=mix, labs=labs)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["compute", "analyze"])
    ap.add_argument("--power", type=float, default=80.0)
    ap.add_argument("--height", type=float, default=3.0)
    ap.add_argument("--triplet", type=float, default=W_TRIPLET, help="weight of the triplet (3 sublevels together)")
    ap.add_argument("--singlet", type=float, default=W_SINGLET, help="weight of the singlet")
    a = ap.parse_args()
    if a.stage == "compute":
        compute(a.power, a.height)
    else:
        analyze(a.power, a.height, a.triplet, a.singlet)
