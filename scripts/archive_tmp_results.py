"""Write the results that so far lived only in tmp/ into the two scan files (schema centrex-cleanup-scan-v1).

    python scripts/archive_tmp_results.py

New runs in SPB_J1_cleanup_scans.h5 (CU1):
    eom_beta_detuning_vx_cu1_spb_4p8cm          beta x detuning x vx at y0 = 0 (NaN = not computed)
    doppler_cu1_spb_4p8cm_coarse_settings       vx x y0 coarse grids for the 9 compared (beta, detuning) settings
    detuning_vx_y_cu1_spb_3cm_coarse            detuning x vx x y0, 3 cm beam
    trajectories_cu1_vx0_19p0_20p2MHz_4p8cm     level populations vs time for each starting J=1 sublevel
    per_level_cu1_20p2MHz_4p8cm                 remaining population for each starting J=1 sublevel
    per_state_cu1_spb_20p2MHz_4p8cm             SPB-related levels: remaining and destinations
New runs in SPB_J2_cleanup_scans.h5 (CU2, 200 V/cm):
    trajectories_cu2_R1_R2_2p4cm                level populations vs time at the two optima
    per_level_cu2_R1_R2_2p4cm                   remaining and leakage to J=1 for each starting J=2 sublevel
    per_state_cu2_spb_R1_R2_2p4cm               SPB source levels: remaining and destinations
    level_scheme_cu2_200V_sublevels             excited sublevels: position, strength, spread (per beam)
    level_scheme_cu2_200V_pairs                 coupled (excited, ground) pairs: position and strength
Existing runs are never overwritten. The checkpoints in tmp/ that are already stored as runs are not repeated.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import cu1_3cm_scan as c3  # noqa: E402
from scripts import cu1_eom_doppler_scan as e1  # noqa: E402
from scripts.cleanup_spb_doppler import CU1_SOLVER, build_cu1  # noqa: E402
from scripts.cu2_laser_aom_scan import CONFIRM_SOLVER, build  # noqa: E402
from scripts.scan_hdf import save_scan_hdf  # noqa: E402

J1 = ROOT / "results" / "SPB_J1_cleanup_scans.h5"
J2 = ROOT / "results" / "SPB_J2_cleanup_scans.h5"
TMP = ROOT / "tmp"
CU1_COMMON = {"electric_field_V_per_cm": [0.0, 0.0, 150.0], "magnetic_field_G": [0.0, 0.0, 1e-3],
              "transition": "R(1) F1'=5/2 F'=2", "laser_power_mW": 80.0, "passes": 13, "loss_per_pass": 0.10,
              "polarization_switching_MHz": 0.5, "molecular_vz_m_per_s": 184.0}
CU2_COMMON = {"electric_field_V_per_cm": [0.0, 0.0, 200.0], "magnetic_field_G": [0.0, 0.0, 1e-3],
              "transitions": ["R(2) F1'=7/2 F'=3", "R(2) F1'=7/2 F'=4"], "powers_mW": [80.0, 24.0],
              "eom_frequency_MHz": 17.5, "eom_beta": 1.5, "aom": "single pass", "passes": 13, "loss_per_pass": 0.10,
              "polarization_switching_MHz": 0.5, "polarization_scheme": "orthogonal beams, antiphase X/Z",
              "molecular_vz_m_per_s": 184.0, "beam_height_4sigma_cm": 2.4, "transverse_velocity": 0.0, "y0_mm": 0.0}
DOPPLER = "k*vx, sign alternating pass by pass (two-prism multipass)"
saved = []


def save(fname, name, results, **kw):
    try:
        save_scan_hdf(fname, name, np.asarray(results), overwrite=False, **kw)
        saved.append(f"{fname.name}: {name} {np.asarray(results).shape}")
    except Exception as exc:                       # e.g. already stored
        saved.append(f"{fname.name}: {name} NOT written ({exc})")


def labels_of(obe):
    out = []
    for q in obe.QN:
        s = q.largest
        el = s.electronic_state.name
        if s.mF is None or s.F is None:                # compacted spectator level
            out.append(f"{el} J={s.J} (compacted)")
            continue
        P = "" if el == "X" else (f", P={'+' if s.P > 0 else '-'}" if s.P is not None else "")
        out.append(f"{el} J={s.J}, F1={s.F1}, F={s.F}, mF={s.mF:+.0f}{P}")
    return out


def cu1_eom_grid():
    d = e1.load("coarse")
    betas = sorted({k[1] for k in d}); dets = sorted({k[0] for k in d}); vx = e1.VX_COARSE
    g = np.full((len(betas), len(dets), len(vx)), np.nan)
    for i, b in enumerate(betas):
        for j, de in enumerate(dets):
            for k, v in enumerate(vx):
                g[i, j, k] = d.get(e1.key((de, b, v, 0.0)), np.nan)
    save(J1, "eom_beta_detuning_vx_cu1_spb_4p8cm", g,
         axes={"eom_beta": (np.array(betas), "rad"), "laser_detuning_MHz": (np.array(dets), "MHz"),
               "vx_m_per_s": (vx, "m/s")},
         quantity="final_population_remaining_in_X_J1",
         fixed_parameters={**CU1_COMMON, "beam_height_4sigma_cm": 4.8, "eom_frequency_MHz": 1.6, "y0_mm": 0.0,
                           "initial_state": "spb (8 levels with |m_J|=1 content > 0.5, equal weights)",
                           "doppler": DOPPLER, "missing": "NaN where the (beta, detuning) pair was not computed"},
         model={"script": "scripts/cu1_eom_doppler_scan.py (stages coarse, extend)", "levels": 26,
                "source": "tmp/cu1_eom_doppler/coarse_*.npz"},
         solver=dict(CU1_SOLVER), notes="CU1 with a 1.6 MHz phase EOM: modulation depth x detuning x vx (report Sec. 5.8).")


def cu1_settings():
    d = e1.load("yavg")
    pairs = sorted({(k[1], k[0]) for k in d})
    vx, y = e1.VX_COARSE, e1.Y_COARSE
    pairs = [p for p in pairs if all(e1.key((p[1], p[0], v, yy)) in d for v in vx for yy in y)]
    g = np.array([[[d[e1.key((de, b, v, yy))] for yy in y] for v in vx] for b, de in pairs])
    save(J1, "doppler_cu1_spb_4p8cm_coarse_settings", g,
         axes={"setting_index": (np.arange(len(pairs)), "index"), "vx_m_per_s": (vx, "m/s"),
               "y0_mm": (y * 1e3, "mm")},
         quantity="final_population_remaining_in_X_J1",
         fixed_parameters={**CU1_COMMON, "beam_height_4sigma_cm": 4.8, "eom_frequency_MHz": 1.6,
                           "settings": [{"eom_beta": b, "laser_detuning_MHz": de} for b, de in pairs],
                           "initial_state": "spb (8 levels, equal weights)", "doppler": DOPPLER,
                           "y_symmetry": "results for y0 >= 0; symmetric in y0"},
         model={"script": "scripts/cu1_eom_doppler_scan.py (stage yavg)", "levels": 26,
                "source": "tmp/cu1_eom_doppler/yavg_*.npz"},
         solver=dict(CU1_SOLVER), notes="Coarse vx x y0 grids of the CU1 settings compared in report Table 5.4.")


def cu1_3cm_detuning():
    d = c3.load("detuning")
    g = np.array([[[d[c3.key((de, v, yy))] for yy in c3.Y_C] for v in c3.VX_C] for de in c3.DETUNINGS])
    save(J1, "detuning_vx_y_cu1_spb_3cm_coarse", g,
         axes={"laser_detuning_MHz": (c3.DETUNINGS, "MHz"), "vx_m_per_s": (c3.VX_C, "m/s"),
               "y0_mm": (c3.Y_C * 1e3, "mm")},
         quantity="final_population_remaining_in_X_J1",
         fixed_parameters={**CU1_COMMON, "beam_height_4sigma_cm": 3.0, "eom": "none",
                           "initial_state": "spb (8 levels, equal weights)", "doppler": DOPPLER,
                           "y_symmetry": "results for y0 >= 0; symmetric in y0"},
         model={"script": "scripts/cu1_3cm_scan.py (stage detuning)", "levels": 26,
                "source": "tmp/cu1_3cm/detuning_*.npz"},
         solver=dict(CU1_SOLVER), notes="CU1 with the 3 cm beam: detuning scan for the acceptance mean (report Table 5.5).")


def cu1_trajectories(labels):
    z = np.load(TMP / "cu1_trajectories" / "cu1_trajectories.npz")
    gidx = z["ground_idx"]
    save(J1, "trajectories_cu1_vx0_19p0_20p2MHz_4p8cm", z["pops"],
         axes={"laser_detuning_MHz": (z["detunings"], "MHz"), "start_level_index": (gidx, "index"),
               "time_us": (z["t"] * 1e6, "us"), "level_index": (np.arange(z["pops"].shape[-1]), "index")},
         quantity="population_of_each_level_vs_time",
         fixed_parameters={**CU1_COMMON, "beam_height_4sigma_cm": 4.8, "eom": "none", "vx_m_per_s": 0.0,
                           "y0_mm": 0.0, "level_labels": labels,
                           "start_level_spb_related": [bool(s) for s in z["spb"]],
                           "photons_per_molecule_spb_start": z["photons_spb"],
                           "photons_per_molecule_uniform_start": z["photons_uniform"],
                           "start": "one trajectory per starting J=1 sublevel (start_level_index)"},
         model={"script": "scripts/cu1_trajectories.py", "levels": 26,
                "source": "tmp/cu1_trajectories/cu1_trajectories.npz"},
         solver={**CU1_SOLVER, "output_when": "saveat", "saveat_step_us": 0.1},
         notes="CU1 on-axis trajectories at +19.0 and +20.2 MHz (report Fig. 5.7).")


def cu1_per_level():
    r = json.loads((TMP / "level_diagrams" / "residuals.json").read_text())["cu1"]
    lv = r["levels"]
    save(J1, "per_level_cu1_20p2MHz_4p8cm", np.array([x["residual"] for x in lv]),
         axes={"start_level_index": (np.array([x["idx"] for x in lv]), "index")},
         quantity="final_population_remaining_in_X_J1",
         fixed_parameters={**CU1_COMMON, "beam_height_4sigma_cm": 4.8, "laser_detuning_MHz": r["detuning_MHz"],
                           "vx_m_per_s": 0.0, "y0_mm": 0.0,
                           "labels": [f"F1={x['F1']}, F={x['F']}, mF={x['mF']:+d}" for x in lv],
                           "spb_related": [x["spb"] for x in lv], "p_mJ": [x["p_mj"] for x in lv]},
         model={"script": "scripts/cleanup_level_diagrams.py (stage compute)", "levels": 26,
                "source": "tmp/level_diagrams/residuals.json"},
         solver=dict(CU1_SOLVER), notes="Remaining J=1 population for a molecule starting in each sublevel, on axis.")


def cu1_per_state():
    lv = json.loads((TMP / "spb_states_doppler" / "spb_states.json").read_text())["cu1"]
    cols = ["remaining", "J0", "J2", "J3", "J4"]
    save(J1, "per_state_cu1_spb_20p2MHz_4p8cm", np.array([[x[c] for c in cols] for x in lv]),
         axes={"state_index": (np.arange(len(lv)), "index"), "quantity_index": (np.arange(len(cols)), "index")},
         quantity="final_population (columns: remaining in J=1, then in J=0, 2, 3, 4)",
         fixed_parameters={**CU1_COMMON, "beam_height_4sigma_cm": 4.8, "laser_detuning_MHz": 20.2,
                           "vx_m_per_s": 0.0, "y0_mm": 0.0, "columns": cols,
                           "states": [{k: x[k] for k in ("F1", "F", "mF", "mJ", "m_Tl", "m_F19", "p_dominant",
                                                          "p_absmJ1")} for x in lv]},
         model={"script": "scripts/cleanup_spb_doppler.py (stage states)", "levels": 26,
                "source": "tmp/spb_states_doppler/spb_states.json"},
         solver=dict(CU1_SOLVER), notes="The 8 SPB-related J=1 levels, one trajectory each, on axis.")


def cu2_trajectories(labels):
    z = np.load(TMP / "cu2_laser_aom_scan" / "optimum_trajectories.npz", allow_pickle=True)
    names = [str(n) for n in z["names"]]
    save(J2, "trajectories_cu2_R1_R2_2p4cm", z["pops"],
         axes={"f_AOM_MHz": (z["points"][:, 1], "MHz"), "time_us": (z["t"] * 1e6, "us"),
               "level_index": (np.arange(z["pops"].shape[-1]), "index")},
         quantity="population_of_each_level_vs_time",
         fixed_parameters={**CU2_COMMON, "settings": names, "delta0_MHz": z["points"][:, 0],
                           "delta1_MHz": z["detuning1"], "Delta_main_MHz": float(z["delta_main"]),
                           "photons_per_molecule": z["photons"], "level_labels": labels,
                           "initial_state": "uniform over the 20 J=2 sublevels"},
         model={"script": "scripts/cu2_optimum_trajectories.py", "levels": 56,
                "source": "tmp/cu2_laser_aom_scan/optimum_trajectories.npz"},
         solver={**CONFIRM_SOLVER, "output_when": "saveat", "saveat_step_us": 0.1},
         notes="CU2 trajectories at the two 200 V/cm optima (report Figs. 6.5, 6.6).")


def cu2_per_level():
    r = json.loads((TMP / "level_diagrams" / "residuals.json").read_text())["cu2"]
    lv, regs = r["levels"], ["R1", "R2"]
    g = np.array([[[x["residual"][reg], x["to_J1"][reg]] for x in lv] for reg in regs])
    save(J2, "per_level_cu2_R1_R2_2p4cm", g,
         axes={"f_AOM_MHz": (np.array([r["points_MHz"][reg][1] for reg in regs]), "MHz"),
               "start_level_index": (np.array([x["idx"] for x in lv]), "index"),
               "quantity_index": (np.arange(2), "index")},
         quantity="final_population (columns: remaining in J=2, transferred to J=1)",
         fixed_parameters={**CU2_COMMON, "delta0_MHz": [r["points_MHz"][reg][0] for reg in regs],
                           "Delta_main_MHz": r["delta_main_MHz"], "columns": ["remaining_J2", "to_J1"],
                           "labels": [f"F1={x['F1']}, F={x['F']}, mF={x['mF']:+d}" for x in lv],
                           "spb_source": [x["spb"] for x in lv],
                           "note": "mF < 0 levels copied from their +|mF| partner (symmetry)"},
         model={"script": "scripts/cleanup_level_diagrams.py (stage compute)", "levels": 56,
                "source": "tmp/level_diagrams/residuals.json"},
         solver=dict(CONFIRM_SOLVER), notes="Remaining J=2 population for a molecule starting in each sublevel, on axis.")


def cu2_per_state():
    lv = json.loads((TMP / "spb_states_doppler" / "spb_states.json").read_text())["cu2"]
    regs = ["R1", "R2"]; cols = ["remaining", "J1", "J3", "J4", "J5"]
    states = [x for x in lv if x["region"] == "R1"]
    g = np.array([[[x[c] for c in cols] for x in lv if x["region"] == reg] for reg in regs])
    save(J2, "per_state_cu2_spb_R1_R2_2p4cm", g,
         axes={"f_AOM_MHz": (np.array([208.5, 191.0]), "MHz"), "state_index": (np.arange(len(states)), "index"),
               "quantity_index": (np.arange(len(cols)), "index")},
         quantity="final_population (columns: remaining in J=2, then in J=1, 3, 4, 5)",
         fixed_parameters={**CU2_COMMON, "delta0_MHz": 11.5, "columns": cols,
                           "states": [{k: x[k] for k in ("F1", "F", "mF", "spin")} for x in states]},
         model={"script": "scripts/cleanup_spb_doppler.py (stage states)", "levels": 56,
                "source": "tmp/spb_states_doppler/spb_states.json"},
         solver=dict(CONFIRM_SOLVER), notes="The 4 SPB source levels (J=2, m_J=0), one trajectory each, on axis.")


def cu2_level_scheme():
    z = np.load(TMP / "cu2_laser_aom_scan" / "level_scheme.npz")
    info = json.loads((TMP / "cu2_laser_aom_scan" / "info.json").read_text())
    common = {**CU2_COMMON, "Delta_main_MHz": float(z["delta_main"]), "main_couplings": info["main_couplings"],
              "frame": "positions on each beam's own detuning axis, beam 1 shifted by -Delta_main"}
    mdl = {"script": "scripts/cu2_level_scheme.py (stage compute)", "levels": 56,
           "source": "tmp/cu2_laser_aom_scan/level_scheme.npz"}
    rc = ["beam", "mF_excited", "P_excited", "position_MHz", "relative_strength", "spread_MHz"]
    save(J2, "level_scheme_cu2_200V_sublevels", z["rows"],
         axes={"row_index": (np.arange(len(z["rows"])), "index"), "column_index": (np.arange(len(rc)), "index")},
         quantity="excited-sublevel table", fixed_parameters={**common, "columns": rc}, model=mdl, solver={},
         notes="Coupling-weighted resonance position, summed |M|^2 and spread of each driven excited sublevel.")
    pc = ["row_index", "position_MHz", "relative_strength", "F1_ground", "F_ground", "mF_ground"]
    save(J2, "level_scheme_cu2_200V_pairs", z["pairs"],
         axes={"pair_index": (np.arange(len(z["pairs"])), "index"), "column_index": (np.arange(len(pc)), "index")},
         quantity="coupled-pair table", fixed_parameters={**common, "columns": pc}, model=mdl, solver={},
         notes="Every coupled (excited, ground) pair: resonance position and |M|^2 (X and Z combined).")


def main():
    cu1_eom_grid(); cu1_settings(); cu1_3cm_detuning(); cu1_per_level(); cu1_per_state()
    _, _, obe1 = build_cu1(return_obe=True)
    cu1_trajectories(labels_of(obe1))
    cu2_per_level(); cu2_per_state(); cu2_level_scheme()
    _, _, _, _, obe2 = build(return_obe=True)
    cu2_trajectories(labels_of(obe2))
    print("\n".join(saved))


if __name__ == "__main__":
    main()
