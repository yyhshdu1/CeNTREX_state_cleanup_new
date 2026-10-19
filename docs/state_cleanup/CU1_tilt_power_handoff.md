# CU1 power, beam height and laser tilt; centrex-tlf 0.3.0 — handoff

Status as of 2026-10-19. This covers three pieces of work done after `CU1_test_handoff.md`:

1. the update of centrex-tlf from 0.2.5 to 0.3.0, with local speed and consistency tests;
2. CU1 at 50 mW with a 2 cm beam, compared with the present 80 mW, 3 cm beam;
3. the effect of a tilted CU1 laser (Doppler shift from the forward velocity), for both beams.

Every number is a model prediction (centrex-tlf OBE, single trajectories at v_z = 184 m/s). Efficiencies are
"mean in acceptance": averaged over v_x (Gaussian, sigma = 0.65 m/s, |v_x| <= 1 m/s) and y0 (Gaussian,
sigma = 10.95 mm at CU1, |y0| <= 8.5 mm, the 1.7 cm detection laser), as in report Sec. 5.10. Start: the 8
SPB-related J=1 levels with equal weights. No EOM. 13 passes, 10 %/pass loss, X/Z switching at 0.5 MHz.

## 1. centrex-tlf 0.3.0

### What was done

- Conda env `CeNTREX_state_cleanup_new` upgraded to centrex-tlf **0.3.0** (PyPI wheel; no GPU on this machine, so
  the new CUDA solvers are not used).
- Upstream `ograsdijk/CeNTREX-TLF` main (v0.3.0) merged into this repo (commit `e7cc035`; only conflict
  `.gitignore`, both sides kept).
- The repo's own `centrex_tlf/` is what `scripts/` import (they put the repo root first on `sys.path`). Its compiled
  extension is not tracked by git and there is no Rust toolchain here, so the 0.3.0 extension was copied from the
  installed wheel:

  ```
  copy <env>\Lib\site-packages\centrex_tlf\centrex_tlf_rust.cp311-win_amd64.pyd  centrex_tlf\
  ```

  **After a fresh clone, or after every centrex-tlf upgrade, repeat this copy**, otherwise the scripts fail to
  import or mix versions. Notebooks import the installed wheel.

### Local tests (same machine, 20 logical cores)

| workload | 0.2.5 | 0.2.10 | 0.3.0 |
|---|---|---|---|
| CU1, 3 cm, +19 MHz, 44 trajectories, 8 threads | 63.5 s | 65.5 s | **56.7 s** |
| same, 20 threads | 48.0 s | 47.6 s | **41.7 s** |
| CU1 result vs stored map | identical | max diff 1.2e-6 | max diff 1.2e-6 |
| static R(2) F'=4 spectrum, 170 V/cm, 41 levels, 241 detunings, dopri5 | 50.1 s | 57.1 s | not run |
| same with `solver="dense_eig"` (0.2.9+) | — | 32.8 s (agrees with dopri5 to 5.5e-12) | not run |

- 0.3.0 makes our time-dependent CU1 runs ~11-13 % faster and leaves the CU1 results unchanged.
- `dense_eig` (the 42x speed-up in the upstream changelog) works only for **time-independent** problems.
  - CU1, CU2 and the probe model are time dependent (pass profile, polarization switching, EOM, Gaussian probe),
    so they cannot use it.
  - On a static spectrum with one initial state and final-only output it gave 1.5x here; the gain grows with
    many initial states and many sample times.

### Physics changes that matter

- **J=2 ground energies changed (0.2.8 fix: lower-J Stark mixing restored in the X Hamiltonian).**
  - The R(2) F'=4 spectrum at 170 V/cm moves by +0.5 to +1.0 MHz: peaks -9.0 -> -8.5, 0 -> +1.0, +19 -> +20,
    +28 -> +29 MHz on a 0.5 MHz grid.
  - Peak heights change by up to +-6 % and the main coupling by -2.4 %; the excited-state energies are unchanged.
  - J=1 (CU1) is unaffected.
  - **All CU2 results in `results/SPB_J2_cleanup_scans.h5` were computed with 0.2.5 and should be re-run** before
    being trusted at the < 1 MHz level. Theory curves for the R2F4 field scan should be recomputed as well.
- **Photon integrals with `output_when="saveat"`**: in 0.2.5 they were trapezoids on the save points (two save
  points = one trapezoid, badly wrong). Fixed in 0.2.6/0.2.7. Our stored CU1/CU2 results use final populations and
  the probe yields use final-only integrals, so they are not affected.

## 2. CU1 at 50 mW, 2 cm beam (sigma_y = 5 mm)

Script `scripts/cu1_conditions_scan.py` (`--power`, `--height`), centrex-tlf 0.3.0.

Detuning scan (coarse grid, mean in acceptance):

| detuning | on axis | mean in acceptance |
|---|---|---|
| +18.5 MHz | 96.32 % | 95.27 % |
| **+19.0 MHz** | 98.78 % | **96.52 %** |
| +19.5 MHz | 99.80 % | 96.44 % |
| +20.5 MHz | 99.998 % | 94.62 % |
| +22.0 MHz | 99.66 % | 96.24 % |

The best detuning is again +19.0 MHz. Full (v_x, y0) map at +19.0 MHz compared with the present beam:

| | 80 mW, 3 cm | 50 mW, 2 cm |
|---|---|---|
| on axis | 99.01 % | 98.78 % |
| **mean in acceptance** | **98.70 %** | **96.52 %** |
| v_x = 0, y0 = 6 / 7.5 / 9 mm | 97.47 / 96.11 / 93.94 % | 93.21 / 87.64 / 78.43 % |
| y0 = 0, v_x = 0.5 / 1 m/s | 99.99 / 94.72 % | 99.99 / 93.89 % |

- The J=1 residual is ~2.7x larger with the 50 mW, 2 cm beam.
- The velocity dependence hardly changes. The loss is at the edges of the acceptance:
  - the peak intensity is similar (0.94x: 50/80 x 7.5/5);
  - at y0 = 8.5 mm, however, the 2 cm beam has 24 % of its peak intensity, against 53 % for the 3 cm beam.
- At 50 mW a 3 cm beam reaches only 96.36 % on axis, so the best height for 50 mW lies between 2 and 3 cm
  (not yet scanned).

Figure: `docs/state_cleanup/figures/cleanup_map_vx_y_CU1_2cm_50mW_vs_3cm_80mW.{png,pdf}`.

## 3. Tilted CU1 laser

### Physics

A tilt theta of the laser from perpendicular, in the plane containing the molecular beam, adds

    k v_z sin(theta) / 2 pi  ~  v_z theta / lambda  =  0.68 MHz per mrad  (~12 MHz per degree)

to every pass. Two cases:

- **Antiparallel passes** (retro-reflecting two-prism multipass, as modelled for CU1):
  - The return pass is antiparallel, so the shift alternates in sign with the transverse Doppler shift:
    k.v = +-k (v_x + v_z tan theta).
  - A tilt is therefore exactly a shift of every molecule's effective transverse velocity by
    **Delta v_x = v_z theta = 0.18 m/s per mrad**, while detection still selects on the true v_x.
  - **Re-tuning the laser cannot remove it.**
- **Parallel passes** (the z component of k has the same sign in every pass, e.g. a plane-mirror zig-zag):
  - A common detuning offset; re-tuning the laser removes it completely.
  - If the laser is not re-tuned, the loss is read off the detuning curve.

The relevant angle is between k and the **mean molecular trajectory**, so a molecular-beam pointing error (mean
v_x != 0, e.g. from lens misalignment) has the same effect, at 0.18 m/s of mean v_x per mrad.

### Results, antiparallel passes (our CU1)

| theta (mrad) | shift per pass | 80 mW, 3 cm, +19 MHz | 80 mW, 3 cm, re-tuned | 50 mW, 2 cm, +19 MHz | 50 mW, 2 cm, re-tuned |
|---|---|---|---|---|---|
| -5 | -3.39 MHz | 88.62 % | >= 92.80 % | 83.95 % | >= 88.56 % |
| -3 | -2.03 MHz | 94.85 % | >= 95.97 % | 91.49 % | 92.71 % |
| -2 | -1.35 MHz | 96.89 % | 97.29 % | 94.12 % | 94.64 % |
| -1 | -0.68 MHz | 98.17 % | 98.21 % | 95.85 % | 95.85 % |
| 0 | 0 | **98.66 %** | 98.66 % | **96.52 %** | 96.52 % |
| +1 | +0.68 MHz | 98.38 % | 98.38 % | 96.16 % | 96.16 % |
| +2 | +1.35 MHz | 97.27 % | 97.27 % | 94.71 % | 94.71 % |
| +3 | +2.03 MHz | 95.41 % | 95.90 % | 92.32 % | 92.69 % |
| +5 | +3.39 MHz | 89.46 % | >= 93.75 % | 85.12 % | >= 89.92 % |

">=": the best detuning sits at the edge of the +17 to +23 MHz scan, so the true optimum may be higher.

- **Tolerance: |theta| <~ 1 mrad (0.06 deg, 0.3 mm over 30 cm).**
  - 80 mW, 3 cm: the loss is 0.3-0.5 % at 1 mrad, 1.4-1.8 % at 2 mrad and 3.3-3.8 % at 3 mrad.
  - 50 mW, 2 cm is somewhat more sensitive: 0.4-0.7 % at 1 mrad and 4.2-5.0 % at 3 mrad.
- Re-tuning helps only beyond ~3 mrad, by ~0.5-1 %.
- Why it hurts: +19.0 MHz relies on the two efficient velocity bands near |v_x| = 0.5 m/s. A tilt slides the whole
  velocity distribution across them.

### Parallel passes (80 mW, 3 cm, laser kept at +19.0 MHz; for reference only)

The loss equals a detuning error of 0.68 MHz/mrad:

| theta | -2 mrad | -1 mrad | +1 mrad | +2 mrad | +3 mrad | +5 mrad |
|---|---|---|---|---|---|---|
| efficiency | 94.22 % | 97.34 % | 98.21 % | 97.27 % | 97.05 % | 98.40 % |

- The rise at +5 mrad is the laser reaching the second optimum near +22.5 MHz.
- The values below -2.9 mrad lie outside the scanned detunings.
- With the laser re-tuned there is no loss. This case was not computed for 50 mW.

### How to measure the tilt

- **Low-power CU1 depletion spectrum, one pass at a time.** Block the beam after the first pass, then send it in
  the opposite direction. The two line centres differ by 2 v_z theta / lambda = 1.35 MHz per mrad.
- **TOF-resolved line centre.** Only the tilt term scales with v_z, so fast and slow molecules show different
  shifts.

Figures: `docs/state_cleanup/figures/CU1_laser_tilt.{png,pdf}` (80 mW, 3 cm, with the parallel case) and
`CU1_laser_tilt_2cm_50mW.{png,pdf}` (antiparallel only, with the 80 mW, 3 cm curve as reference).

## 4. Assumptions and limitations

- **Multipass geometry.** Antiparallel return passes are the model's assumption for CU1. Confirm it on the real
  optics: if the passes are not retro-reflected, the parallel or a mixed case applies.
- **Neglected effects.**
  - The forward-velocity spread, i.e. the spread of v_z theta: small for a few mrad compared with sigma_vx.
  - Changes of pass position with tilt.
  - A tilt out of the v_z plane: its Doppler term ~ v_y phi is negligible.
- **Grids.** Coarse: v_x step 0.2 m/s, y0 = 0, 3, 6, 9 mm, detuning step 0.5 MHz, scan +17 to +23 MHz. Effective
  velocities beyond +-1 m/s were computed for the tilt study (to +-2 m/s).
- **Mixed versions.** The 80 mW, 3 cm cube below |v_x| = 1 m/s is the stored 0.2.5 result; everything else is
  0.3.0. CU1 agrees between versions to 1e-6.
- **Distribution and start state.** The transverse distribution is the lensed J=2 beam simulation, and the start
  is the 8 SPB-related levels.
  - The AP1-based CU1 test of `CU1_test_handoff.md` populates the 4 mainly-m_J=0 levels instead.
  - With the lens off or on, J=1 has its own transverse distribution.

## 5. Open items

1. Re-run CU2 (and any J=2 spectra) with 0.3.0; the J=2 ground energies moved by up to ~1 MHz.
2. Beam-height scan at 50 mW (e.g. 2, 2.5, 3, 3.5 cm) to find the best height for the available power.
3. Check the CU1 multipass geometry, then measure the actual tilt (Sec. 3, "How to measure the tilt").
4. Extend the detuning scan beyond +17 to +23 MHz if tilts above ~3 mrad matter (the ">=" entries).
5. Optional: include the forward-velocity distribution, i.e. average over v_z for the tilt.

## 6. Files and commands

| what | where |
|---|---|
| CU1 at a chosen power and height | `scripts/cu1_conditions_scan.py` (`detuning`, `map D`, `plot D`; `--power`, `--height`) |
| tilt study | `scripts/cu1_tilt_scan.py` (`compute`, `archive`, `analyze`, `plot`; `--power`, `--height`, `--parallel`) |
| 50 mW, 2 cm runs | `results/SPB_J1_cleanup_scans.h5`: `detuning_vx_y_cu1_spb_2cm_50mW_coarse`, `doppler_cu1_spb_2cm_50mW_det19MHz`, `detuning_vx_y_cu1_spb_2cm_50mW_coarse_vx2` |
| 80 mW, 3 cm extended cube | same file: `detuning_vx_y_cu1_spb_3cm_coarse_vx2` |
| index of all runs | `results/scan_results.md` (regenerate with `python scripts/write_results_index.py`) |
| figures | `docs/state_cleanup/figures/cleanup_map_vx_y_CU1_2cm_50mW_vs_3cm_80mW`, `CU1_laser_tilt`, `CU1_laser_tilt_2cm_50mW` (`.png`, `.pdf`) |
| upstream | `git remote upstream` = https://github.com/ograsdijk/CeNTREX-TLF ; changelog in `CHANGELOG.md` |

Reproduce, for example:

```
python scripts/cu1_conditions_scan.py detuning --power 50 --height 2
python scripts/cu1_conditions_scan.py map 19.0 --power 50 --height 2
python scripts/cu1_tilt_scan.py compute --power 50 --height 2
python scripts/cu1_tilt_scan.py plot --power 50 --height 2
python scripts/cu1_tilt_scan.py plot --parallel          # 80 mW, 3 cm, with the parallel case
```

Checkpoints of these runs are in `tmp/` (git-ignored); the results needed for the tables are in the HDF5 file.
