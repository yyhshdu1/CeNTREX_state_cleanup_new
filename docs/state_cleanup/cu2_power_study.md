# CU2 at the real laser power — running log

## Summary so far (updated 2026-11-01; final summary in `cu2_power_summary.md`)

Model: 56-level CU2 OBE, 200 V/cm, 13 passes, 5 % loss per pass, vz = 184 m/s, R(2) F'=3 and F'=4 from one laser
and a single-pass AOM, EOM 17.5 MHz β = 1.5, X/Z switching. Figure of merit: cleanup efficiency averaged over the
detection acceptance (|vx| ≤ 1 m/s with σvx = 0.65 m/s, |y0| ≤ 8.5 mm).

| Question | Answer |
|---|---|
| Best split of 60 mW | At the new set point (full mean): SPB 35 + 25 = 40 + 20 = 99.50 %, 30 + 30 99.44 %; uniform 25 + 35 97.82 %, 30 + 30 97.66 %, 35 + 25 97.26 %. **35 + 25 mW if the start is mostly m_J = 0; 30 + 30 mW if >13 % is in other m_J.** |
| Best set point | **δ0/2π = 12.5 MHz, f_AOM = 207.5 MHz** (laser +1 MHz, AOM −1 MHz from R1; F'=4 beam unchanged). SPB 99.50 %, uniform 97.26 % vs R1 99.22 %, 96.66 %. (13, 206) is equivalent (99.47 / 97.41 %). |
| Robust to the vx width? | Yes: better than R1 for σvx = 0.4 … 1.28 m/s. |
| EOM depth | β = 1.5 is optimal (flat to 1.75). |
| X/Z switching | 0.5 MHz is fine: 0.1–2 MHz changes ≤ 0.1 % (SPB), ≤ 0.4 % (uniform); the phase does not matter. |
| Beam height | 2.4 cm is best (SPB 99.58 %, uniform 97.68 %), 3 cm loses 0.08 / 0.42 %. Not worth a big change. |
| Field tolerance | ±2.5 V/cm costs ≤ 0.1 %; ±5 V/cm up to 1 %. Know the field to ~±2 V/cm or find the set point by an in-situ laser scan. |
| Loss per pass | Matters: 10 % → SPB 99.16 %, uniform 95.74 % (proxy, vs 99.64 / 97.87 % at 5 %). The optimal set point does not move at 10 %. |
| More power | 80 mW (47 + 33) → SPB 99.88 %, uniform 99.13 % (proxy). |
| Old model (80 + 24 mW, 10 %, R1) | SPB 99.57 %, uniform 96.66 %. At 35 + 25 mW, 5 % and the new set point: 99.50 / 97.26 %, i.e. about the same. |

"Proxy" values (6 points: vx = −0.6, 0, 0.6 m/s × y0 = 0, 6 mm) read 0.1–0.7 % higher than the full acceptance mean
but rank settings correctly; compare proxy with proxy only.

Started 2026-10-30. Script: `scripts/cu2_power_check.py`, checkpoints in `tmp/cu2_power/` (all of these runs used v_z = 184 m/s; since 2026-11-09 the models default to the measured 156 m/s and the 184 m/s checkpoints are in `tmp/cu2_power/vz184/`).

## Why

All stored CU2 results (report Sec. CU2, R1/R2 set points, Doppler maps) assume **80 mW (F'=3) + 24 mW (F'=4)**
and **10 % loss per multipass pass**. In the lab the beams are **35 + 25 mW**, at most **60 mW in total**, and the
windows lose only about **5 % per pass**.

## Model (unchanged unless stated)

- 56-level CU2 OBE of `scripts/cu2_laser_aom_scan.py`: R(2) F1'=7/2 F'=3 and F'=4, 200 V/cm, B = 1 mG,
  one 17.5 MHz EOM (β = 1.5) on both beams, X/Z polarization switching at 0.5 MHz, vz = 184 m/s.
- 13 passes, **5 % loss per pass**, 3 cm beam (σy = 7.5 mm), σz = 0.75 mm.
- Set points from the 80 + 24 mW, 2.4 cm optimisation: R1 = (δ0/2π = 11.5 MHz, f_AOM = 208.5 MHz),
  R2 = (11.5 MHz, 191.0 MHz). Δ_main = 211.829 MHz with centrex-tlf 0.3.0, same as before.
- Starts: *uniform* = all 20 J=2 sublevels; *SPB* = the four J=2, m_J=0 levels (singlet 2/3, triplets 1/9).
- Acceptance mean as in `scripts/plot_cleanup_maps_current.py`: |vx| <= 1 m/s, |y0| <= 8.5 mm, Gaussian weights
  σvx = 0.65 m/s, σy = 10.35 mm.

## Step 1 — power split for 60 mW, on axis (2026-10-30)

Cleanup efficiency on the beam axis (vx = 0, y0 = 0).

| P(F'=3) + P(F'=4) | R1 uniform | R2 uniform | R1 SPB | R2 SPB |
|---|---|---|---|---|
| 50 + 10 | 97.66 % | 97.74 % | 99.882 % | 99.885 % |
| 45 + 15 | 98.97 % | 99.02 % | 99.936 % | 99.936 % |
| 40 + 20 | 99.49 % | 99.53 % | 99.958 % | 99.955 % |
| 35 + 25 (present) | 99.71 % | 99.74 % | 99.967 % | 99.963 % |
| 30 + 30 | 99.82 % | 99.82 % | **99.972 %** | 99.966 % |
| 25 + 35 | 99.86 % | 99.85 % | **99.972 %** | 99.963 % |
| 20 + 40 | **99.875 %** | 99.84 % | 99.963 % | 99.949 % |
| 80 + 24 (reference) | 99.84 % | 99.89 % | 99.993 % | 99.994 % |

- Power is better spent on the **F'=4** beam. The residual is dominated by the stretched sublevels
  |F1=5/2, F=3, mF=±3>, which are pumped almost only by F'=4.
- Uniform start: still improving at 20 + 40 mW, where it equals 80 + 24 mW.
- SPB start: best at 25–30 mW on F'=3, 2.9e-4 left, about 4x more than at 80 + 24 mW.
- R1 and R2 are within 0.03 %; R1 is slightly better for the SPB start.

### Step 1b — more weight on F'=4 (2026-10-30)

| P(F'=3) + P(F'=4) | R1 uniform | R2 uniform | R1 SPB | R2 SPB |
|---|---|---|---|---|
| 15 + 45 | 99.84 % | 99.76 % | 99.919 % | 99.896 % |
| 10 + 50 | 99.61 % | 99.46 % | 99.743 % | 99.699 % |

On axis the uniform start peaks at 20 + 40 mW and the SPB start at 25–30 mW on F'=3.

## Step 2 — mean over the detection acceptance, R1 (2026-10-30)

Grid vx = −1…1 m/s (0.2 m/s) × y0 = 0, 3, 6, 9 mm; mean as in `plot_cleanup_maps_current.py`.

| Setting | uniform: axis | uniform: acceptance | SPB: axis | SPB: acceptance |
|---|---|---|---|---|
| 80 + 24 mW, 10 % loss, 0.2.5 (stored report map) | 99.52 % | 96.66 % | 99.974 % | 99.57 % |
| 25 + 35 mW, 5 % | 99.86 % | **97.29 %** | 99.972 % | 99.02 % |
| 20 + 40 mW, 5 % | 99.87 % | 97.25 % | 99.96 % | 98.74 % |
| 30 + 30 mW, 5 % | 99.82 % | 97.11 % | 99.972 % | 99.16 % |
| 35 + 25 mW, 5 % (present lab) | 99.71 % | 96.66 % | 99.97 % | 99.22 % |
| 40 + 20 mW, 5 % | – | – | 99.96 % | **99.24 %** |
| 45 + 15 mW, 5 % | – | – | 99.94 % | 99.22 % |

Efficiency versus vx at y0 = 0 (%):

| vx (m/s) | −1 | −0.6 | −0.2 | 0 | 0.2 | 0.6 | 1 |
|---|---|---|---|---|---|---|---|
| uniform, 80 + 24 / 10 % | 95.46 | 97.10 | 99.24 | 99.52 | 99.26 | 97.23 | 95.52 |
| uniform, 25 + 35 | 94.36 | 97.82 | 99.77 | 99.86 | 99.75 | 97.73 | 94.32 |
| uniform, 30 + 30 | 94.42 | 97.72 | 99.70 | 99.82 | 99.68 | 97.66 | 94.39 |
| SPB, 80 + 24 / 10 % | 99.30 | 99.76 | 99.96 | 99.97 | 99.95 | 99.75 | 99.00 |
| SPB, 25 + 35 | 97.61 | 99.30 | 99.95 | 99.97 | 99.91 | 98.79 | 96.96 |
| SPB, 30 + 30 | 97.97 | 99.41 | 99.96 | 99.97 | 99.93 | 98.95 | 97.26 |

- At 60 mW the beam axis is as good as before (better for the uniform start, because 5 % loss keeps the last passes
  brighter), but the tolerance to transverse velocity is lower. With less power the lines are less power broadened,
  and the ±3.7 MHz Doppler shift at 1 m/s moves molecules off resonance.
- The SPB levels (m_J = 0) are pumped mainly by the F'=3 beam, so the SPB start prefers more F'=3 power and the
  uniform start (stretched m_F = ±3 levels, F'=4 only) prefers more F'=4. 25 + 35 and 30 + 30 differ by only
  0.15–0.2 % in either mean.
- These comparisons mix the power and loss change with centrex-tlf 0.2.5 → 0.3.0. On axis at 80 + 24 mW with
  5 % loss and 0.3.0, the SPB start gives 99.993 % (Step 1).

### Proxy for the set-point search

The full grid costs ~1 h per setting and start. For re-optimising (δ0, f_AOM) a 6-point proxy is used:
vx = −0.6, 0, 0.6 m/s × y0 = 0, 6 mm, Gaussian weights in vx and a plain mean in y0. Against the full grids above
it reads 0.1–0.5 % high and ranks them in the same order:

| | full | proxy |
|---|---|---|
| uniform 80+24 / 25+35 / 30+30 | 96.66 / 97.29 / 97.11 | 97.08 / 97.79 / 97.61 |
| SPB 80+24 / 25+35 / 30+30 | 99.57 / 99.02 / 99.16 | 99.68 / 99.26 / 99.34 |

- Over the acceptance the SPB start keeps improving with more F'=3 power (98.74 → 99.02 → 99.16 → 99.22 % from
  20 to 35 mW on F'=3) and is then flat to 45 mW (99.24 % at 40 + 20, 99.22 % at 45 + 15), while the uniform
  start peaks at 25 + 35 mW. Any split from 35 + 25 to 45 + 15 is equally good for the SPB start.
- Which start matters: the lens focuses J=2, m_J=0, so the J=2 population left after SPB is mostly the SPB-type
  start. The uniform start is a worst case.

## Step 3 — (δ0, f_AOM) at 30 + 30 mW, proxy (2026-10-31)

Proxy efficiency (%), rows δ0/2π (MHz), columns f_AOM (MHz). R1 = (11.5, 208.5), R2 = (11.5, 191).

SPB start:

| δ0 \ f_AOM | 186 | 191 | 197 | 203 | 208.5 | 214 |
|---|---|---|---|---|---|---|
| 5.0 | 94.70 | 96.92 | 96.84 | 95.24 | 97.31 | 97.68 |
| 8.0 | 98.46 | 98.07 | 97.94 | 98.55 | 98.80 | 98.42 |
| 11.5 | 99.02 | 99.08 | 98.73 | 99.23 | **99.34** | 98.92 |
| 15.0 | 97.07 | 96.42 | 97.03 | 97.40 | 96.66 | 97.27 |
| 18.5 | 96.46 | 96.60 | 97.45 | 97.11 | 97.31 | 97.36 |

Uniform start:

| δ0 \ f_AOM | 186 | 191 | 197 | 203 | 208.5 | 214 |
|---|---|---|---|---|---|---|
| 5.0 | 83.32 | 90.09 | 91.79 | 86.09 | 91.59 | 94.89 |
| 8.0 | 93.06 | 89.85 | 91.08 | 94.35 | 93.38 | 93.56 |
| 11.5 | 94.10 | 96.27 | 94.66 | 95.01 | **97.61** | 94.25 |
| 15.0 | 95.62 | 92.80 | 95.51 | 96.80 | 92.86 | 96.07 |
| 18.5 | 92.48 | 91.42 | 92.93 | 94.50 | 93.66 | 92.85 |

- R1 is still the best point of the coarse grid for both starts; R2 is 0.3 % (SPB) and 1.3 % (uniform) lower.
- δ0 is sharp: ±3.5 MHz from 11.5 MHz costs 0.5–2.7 % (SPB) and up to 4.8 % (uniform). The landscape in f_AOM is
  spiky (EOM sidebands every 17.5 MHz), so the coarse grid can miss a narrow better point near R1.

## Step 4 — fine (δ0, f_AOM) grid around R1 at 35 + 25 mW, proxy (2026-10-31)

SPB start:

| δ0 \ f_AOM | 206 | 207.5 | 208.5 | 209.5 | 211 |
|---|---|---|---|---|---|
| 10.0 | 99.18 | 99.30 | 99.37 | 99.40 | 99.25 |
| 11.0 | 99.17 | 99.25 | 99.30 | 99.21 | 99.13 |
| 11.5 | 99.31 | 99.38 | 99.38 (R1) | 99.30 | 99.18 |
| 12.0 | 99.50 | 99.56 | 99.49 | 99.43 | 99.23 |
| 13.0 | **99.59** | 99.50 | 99.41 | 99.26 | 99.05 |

Uniform start:

| δ0 \ f_AOM | 206 | 207.5 | 208.5 | 209.5 | 211 |
|---|---|---|---|---|---|
| 10.0 | 93.66 | 96.14 | 96.79 | 97.07 | 96.12 |
| 11.0 | 95.74 | 96.81 | 97.16 | 96.63 | 96.31 |
| 11.5 | 96.62 | 97.18 | 97.18 (R1) | 96.79 | 95.77 |
| 12.0 | 97.25 | 97.70 | 97.18 | 97.06 | 94.88 |
| 13.0 | **98.13** | 97.63 | 97.28 | 95.99 | 93.41 |

- Both starts improve towards larger δ0 and smaller f_AOM: (13, 206) MHz is +0.2 % (SPB) and +1.0 % (uniform)
  in the proxy over R1. That is δ1 = +7.2 MHz instead of +8.2 MHz.
- The best point is at the corner of the grid, and the coarse grid was poor at δ0 = 15 MHz, so the optimum lies
  between.

## Step 5 — beyond the fine-grid edge, and a full check at (13, 206) (2026-10-31)

Proxy, 35 + 25 mW. SPB start:

| δ0 \ f_AOM | 203 | 204.5 | 206 | 207.5 |
|---|---|---|---|---|
| 12.5 | 99.52 | 99.57 | 99.61 | **99.64** |
| 13.0 | 99.43 | 99.55 | 99.59 | 99.50 |
| 13.5 | 99.21 | 99.32 | 99.37 | 99.19 |
| 14.0 | 98.91 | 98.97 | 98.96 | 98.77 |
| 15.0 | 97.88 | 97.93 | 97.68 | 97.36 |

Uniform start:

| δ0 \ f_AOM | 203 | 204.5 | 206 | 207.5 |
|---|---|---|---|---|
| 12.5 | 95.35 | 96.92 | 97.71 | 97.87 |
| 13.0 | 95.71 | 97.48 | **98.13** | 97.63 |
| 13.5 | 96.13 | 97.58 | 98.09 | 97.15 |
| 14.0 | 96.58 | 97.44 | 97.39 | 96.63 |
| 15.0 | 96.62 | 96.73 | 95.84 | 94.42 |

The optimum is now enclosed: δ0 = 12.5–13 MHz, f_AOM = 206–207.5 MHz (δ1 = +7.2…+8.2 MHz).

Full acceptance grid at (13, 206) MHz, 35 + 25 mW, against R1 at the same power:

| | SPB: axis | SPB: acceptance | uniform: axis | uniform: acceptance |
|---|---|---|---|---|
| R1 (11.5, 208.5) | 99.97 % | 99.22 % | 99.71 % | 96.66 % |
| (13, 206) | 99.64 % | **99.47 %** | 98.99 % | **97.41 %** |

Efficiency versus vx at y0 = 0 (%), SPB start:

| vx (m/s) | −1 | −0.6 | −0.4 | −0.2 | 0 | 0.2 | 0.4 | 0.6 | 1 |
|---|---|---|---|---|---|---|---|---|---|
| R1 | 98.11 | 99.49 | 99.84 | 99.95 | 99.97 | 99.93 | 99.72 | 99.12 | 97.35 |
| (13, 206) | 98.99 | 99.93 | 99.95 | 99.82 | 99.64 | 99.78 | 99.92 | 99.86 | 98.25 |

- Detuning δ0 by +1.5 MHz from R1 gives up 0.3–0.7 % on the axis and gains it back at |vx| ≈ 0.4 m/s
  (1.5 MHz Doppler shift) and beyond: the response in vx becomes double-peaked and wider.
- The gain does not depend on the assumed vx width. Acceptance mean (SPB / uniform) for a vx Gaussian of
  σvx = 0.4, 0.65, 1.0, 1.28 m/s (1.28 m/s = measured J=1 width):
  R1 99.52 / 99.22 / 99.06 / 99.01 % and 97.59 / 96.66 / 96.18 / 96.03 %;
  (13, 206) 99.52 / 99.47 / 99.42 / 99.40 % and 97.88 / 97.41 / 97.08 / 96.97 %.
  So (13, 206) is equal at 0.4 m/s and better for any wider distribution.
- The proxy read 0.12 % (SPB) and 0.72 % (uniform) high at (13, 206), as on the earlier grids.

## Step 6 — full grid at (12.5, 207.5) and EOM depth (2026-10-31)

Full acceptance grids, 35 + 25 mW, 5 % loss:

| Set point (δ0, f_AOM) MHz | δ1 (MHz) | SPB: axis | SPB: acceptance | uniform: axis | uniform: acceptance |
|---|---|---|---|---|---|
| R1 (11.5, 208.5) | +8.17 | 99.97 % | 99.22 % | 99.71 % | 96.66 % |
| (12.5, 207.5) | +8.17 | 99.88 % | **99.50 %** | 99.52 % | 97.26 % |
| (13, 206) | +7.17 | 99.64 % | 99.47 % | 98.99 % | **97.41 %** |

Acceptance mean for vx widths σvx = 0.4 / 0.65 / 1.0 / 1.28 m/s:

| Set point | SPB | uniform |
|---|---|---|
| R1 | 99.52 / 99.22 / 99.06 / 99.01 | 97.59 / 96.66 / 96.18 / 96.03 |
| (12.5, 207.5) | **99.65** / **99.50** / 99.40 / 99.37 | **97.94** / 97.26 / 96.88 / 96.75 |
| (13, 206) | 99.52 / 99.47 / **99.42** / **99.40** | 97.88 / **97.41** / **97.08** / **96.97** |

- (12.5, 207.5) moves only the F'=3 beam, by +1 MHz (laser +1 MHz, AOM −1 MHz; F'=4 beam unchanged). It keeps
  most of the on-axis efficiency (99.88 %) and is the best SPB point for σvx ≤ 0.65 m/s. (13, 206) moves further
  and is better for the uniform start and for wider vx distributions. The two differ by ≤ 0.15 %.
- Laser-frequency tolerance at (12.5, 207.5) from the proxy grids: δ0 ± 0.5 MHz costs ≤ 0.14 % (SPB).

EOM modulation depth (one EOM for both beams), proxy, 35 + 25 mW:

| Set point | β = 1.0 | 1.25 | 1.5 | 1.75 | 2.0 | 2.5 |
|---|---|---|---|---|---|---|
| SPB (12.5, 207.5) | 99.26 | 99.55 | **99.64** | 99.61 | 99.47 | 98.83 |
| SPB (13, 206) | 99.12 | 99.47 | **99.59** | 99.57 | 99.46 | 98.93 |
| uniform (12.5, 207.5) | 96.07 | 97.32 | **97.87** | 97.76 | 96.82 | 93.71 |
| uniform (13, 206) | 96.33 | 97.62 | **98.13** | 97.89 | 97.03 | 94.28 |

- The present β = 1.5 is optimal (flat to 1.75); no change needed.

## Step 7 — split at the new point, solver check, beam height (2026-10-31)

Power split at (12.5, 207.5), proxy (%):

| P(F'=3) + P(F'=4) | 25 + 35 | 30 + 30 | 35 + 25 | 40 + 20 | 45 + 15 |
|---|---|---|---|---|---|
| SPB | 99.46 | 99.58 | 99.64 | **99.65** | 99.60 |
| uniform | **98.37** | 98.24 | 97.87 | 97.21 | 96.09 |

- Same picture as at R1: the SPB start is flat from 30 + 30 to 45 + 15 mW, the uniform start wants more F'=4.
  35 + 25 mW stays a good choice; 30 + 30 mW trades 0.06 % (SPB) for +0.4 % (uniform).

Solver check, tight tolerances (reltol 1e-7, abstol 1e-9) against the scan tolerances (1e-6, 1e-8), y0 = 0:

| Point | SPB tight / scan | uniform tight / scan |
|---|---|---|
| (12.5, 207.5), vx = 0 | 99.881 / 99.88 | 99.521 / 99.52 |
| (12.5, 207.5), vx = 0.6 | 99.684 / 99.70 | 98.267 / 98.30 |
| (13, 206), vx = 0 | 99.645 / 99.64 | 99.030 / 98.99 |
| (13, 206), vx = 0.6 | 99.870 / 99.86 | 99.110 / 99.09 |

Converged to ≤ 0.04 %, well below the differences between set points.

Beam height at (12.5, 207.5), 35 + 25 mW, same total power (full grids):

| Beam height (4σ) | SPB: axis | SPB: acceptance | uniform: axis | uniform: acceptance |
|---|---|---|---|---|
| 1.6 cm (σy = 4 mm) | 99.99 % | 99.13 % | 99.96 % | 96.25 % |
| 2 cm (σy = 5 mm) | 99.98 % | 99.53 % | 99.91 % | 97.49 % |
| 2.4 cm (σy = 6 mm) | 99.96 % | **99.58 %** | 99.82 % | **97.68 %** |
| 3 cm (σy = 7.5 mm, present) | 99.88 % | 99.50 % | 99.52 % | 97.26 % |
| 4 cm (σy = 10 mm) | 99.64 % | 99.24 % | 98.59 % | 95.78 % |

- A smaller beam is better: at 60 mW the gain from higher intensity outweighs the loss at the edge of the
  acceptance. Mean over vx per y0 (SPB): 2.4 cm 99.85 / 99.76 / 99.36 / 97.89 % at y0 = 0 / 3 / 6 / 9 mm,
  3 cm 99.67 / 99.59 / 99.30 / 98.48 %.
- Holds for every vx width: 2.4 cm is +0.05 … +0.10 % (SPB) and +0.3 … +0.5 % (uniform) over 3 cm for
  σvx = 0.4 … 1.28 m/s.

- The optimum is near 2.4 cm (σy = 6 mm); 2 cm is close, 1.6 cm loses the edge of the acceptance
  (SPB 88.2 %, uniform 66.7 % at y0 = 9 mm). Acceptance means for σvx = 0.4 / 0.65 / 1.0 / 1.28 m/s:

| Beam height | SPB | uniform |
|---|---|---|
| 1.6 cm | 99.26 / 99.13 / 99.03 / 99.00 | 96.84 / 96.25 / 95.92 / 95.81 |
| 2 cm | 99.64 / 99.53 / 99.46 / 99.43 | 98.03 / 97.49 / 97.18 / 97.08 |
| 2.4 cm | **99.70 / 99.58 / 99.50 / 99.47** | **98.24 / 97.68 / 97.35 / 97.25** |
| 3 cm | 99.65 / 99.50 / 99.40 / 99.37 | 97.94 / 97.26 / 96.88 / 96.75 |
| 4 cm | 99.43 / 99.24 / 99.12 / 99.09 | 96.71 / 95.78 / 95.27 / 95.10 |

- Caveat: the best height is tied to the detection acceptance |y0| ≤ 8.5 mm (detection laser 1.7 cm high).
  A taller detection region would favour a taller beam.
- Edge check: an extra y0 = 7.5 mm row (all vx) for 2, 2.4 and 3 cm changes the acceptance means by ≤ 0.01 %
  (SPB) and ≤ 0.08 % (uniform, 2 cm); the log interpolation between 6 and 9 mm is adequate.
- Set point at 2.4 cm (proxy, δ0 = 11.5–13 MHz × f_AOM = 206–208.5 MHz): best SPB point still (12.5, 207.5)
  (99.73 %), best uniform point still (13, 206) (98.58 %). The set point does not depend on the beam height.

## Step 8 — E-field drift: line shifts (2026-10-31)

Stark shift of R(2) lines near 200 V/cm (reduced X and B Hamiltonians, B = 1 mG), per +10 V/cm:

| Line | shift per +10 V/cm |
|---|---|
| F'=3 main: X F=3, mF=2 → B F'=3, mF'=3 | −2.19 MHz |
| X F=3, mF=3 → B F'=4, mF'=4 (stretched) | −2.11 MHz |
| X F1=5/2, F=2, mF=0 → B F'=3, mF'=1 (SPB level) | −0.78 MHz |
| X F1=3/2, F=1, mF=0 → B F'=3, mF'=1 (SPB level) | −0.75 MHz |
| F'=4 main: X F=3, mF=0 → B F'=4, mF'=1 | −0.57 MHz |

The lines shift by different amounts, so a field change cannot be fully compensated by re-tuning the laser.
With δ0 ± 0.5 MHz costing ≤ 0.14 % (Step 6), field drifts of a few V/cm should be harmless; ±5 and ±10 V/cm are
being simulated directly.

## Step 9 — E-field drift and forward velocity, simulated (2026-11-01)

Model rebuilt at each field; (12.5, 207.5) MHz, 35 + 25 mW, 3 cm beam; proxy (%). At 200 V/cm: SPB 99.64,
uniform 97.87. "Laser fixed" = laser and AOM frequency unchanged while the field is different (δ0 → δ0 − shift of
the F'=3 main line; Δ_main from the rebuilt model agrees with the computed shifts to 1 kHz). "Re-tuned" = same δ0
relative to the F'=3 main line at the new field.

| Field (V/cm) | SPB, laser fixed | SPB, re-tuned | uniform, laser fixed | uniform, re-tuned |
|---|---|---|---|---|
| 190 | 99.40 | 98.76 | 96.56 | 95.07 |
| 195 | 99.44 | 99.03 | 97.29 | 96.49 |
| 197.5 | 99.53 | 99.46 | 97.70 | 97.42 |
| **200** | **99.64** | **99.64** | **97.87** | **97.87** |
| 202.5 | 99.56 | 99.54 | 97.96 | 97.64 |
| 205 | 98.65 | 99.06 | 97.64 | 96.74 |
| 210 | 99.18 | 99.29 | 97.11 | 96.26 |

SPB start, on axis (vx = 0, y0 = 0), laser fixed: 190 V/cm 99.93 %, 195 99.94 %, 200 99.88 %, 205 99.05 %,
210 98.82 %.

- A field error of 5 V/cm costs 0.2–1.0 % (SPB) and up to 0.6 % (uniform) with the laser fixed. The response is
  not monotonic (205 V/cm is worse than 210 V/cm for the SPB start), as expected from a spiky (δ0, f_AOM) landscape.
- Re-tuning the laser to the F'=3 main line does not help: the lines that pump the SPB levels shift by only
  −0.78 MHz per 10 V/cm against −2.19 MHz for the F'=3 main line, so "re-tuned" is mostly mis-tuned for them.
- Field too high is worse than too low on axis, because (12.5, 207.5) already sits 1 MHz above R1; a higher field
  moves the SPB-level lines further from the laser (δ → δ + 0.78 MHz per 10 V/cm).
- **Practical consequence:** the CU2 field should be known to ~±2 V/cm, or the set point should be found by
  scanning the laser in situ at the actual field. ±2.5 V/cm costs ≤ 0.11 % (SPB) and ≤ 0.17 % (uniform).

Forward velocity at the same set point (proxy):

| vz (m/s) | SPB | uniform |
|---|---|---|
| 150 | 99.82 | 98.85 |
| 184 (model) | 99.64 | 97.87 |
| 220 | 99.38 | 96.65 |

Slower molecules spend longer in the beam and are cleaned better; the change is smooth, so averaging over the
forward-velocity distribution changes the mean little.

## Step 10 — loss per pass and more power (2026-11-01)

(12.5, 207.5) MHz, 3 cm beam, proxy (%).

| Loss per pass | 3 % | **5 %** | 8 % | 10 % | 15 % |
|---|---|---|---|---|---|
| SPB | 99.76 | **99.64** | 99.39 | 99.16 | 98.41 |
| uniform | 98.48 | **97.87** | 96.70 | 95.74 | 92.84 |

| Power (F'=3 + F'=4) | 35 + 25 (60 mW) | 47 + 33 (80 mW) | 58 + 42 (100 mW) | 80 + 24 (104 mW) |
|---|---|---|---|---|
| SPB | 99.64 | 99.88 | 99.95 | 99.93 |
| uniform | 97.87 | 99.13 | 99.64 | 99.08 |

- The loss matters as much as the power: each extra % of loss per pass costs ~0.08 % (SPB) and ~0.3 % (uniform).
  At 10 % the SPB start is back to ~99.2 %; at the 15 % measured earlier, 98.4 %.
- 20 mW more power (80 mW total, same ratio) would gain +0.24 % (SPB) and +1.3 % (uniform). At equal total
  power the 7:5 ratio beats 80 + 24 for the uniform start.

## Step 11 — does the set point move at 10 % loss? (2026-11-01)

Proxy (%), 35 + 25 mW, 3 cm beam, 10 % loss per pass; rows δ0/2π (MHz), columns f_AOM (MHz). In brackets the
same point at 5 % (Steps 4–5).

SPB start:

| δ0 \ f_AOM | 206.0 | 207.5 | 208.5 |
|---|---|---|---|
| 11.5 | 98.80 (99.31) | 98.88 (99.38) | 98.87 (99.38) |
| 12.0 | 99.02 (99.50) | 99.11 (99.56) | 98.98 (99.49) |
| 12.5 | 99.13 (99.61) | **99.16** (**99.64**) | 99.00 (—) |
| 13.0 | 99.05 (99.59) | 98.92 (99.50) | 98.77 (99.41) |

Uniform start:

| δ0 \ f_AOM | 206.0 | 207.5 | 208.5 |
|---|---|---|---|
| 11.5 | 94.02 (96.62) | 94.72 (97.18) | 94.78 (97.18) |
| 12.0 | 94.82 (97.25) | 95.53 (97.70) | 94.68 (97.18) |
| 12.5 | 95.44 (97.71) | 95.74 (97.87) | 94.85 (—) |
| 13.0 | **96.08** (**98.13**) | 95.36 (97.63) | 94.83 (97.28) |

- **The optimum does not move.** At 10 % loss the SPB start still peaks at (12.5, 207.5) and the uniform start
  at (13, 206), exactly as at 5 %. The whole map drops by ~0.5 % (SPB) and ~2–2.5 % (uniform).
- The gain of the new point over R1 = (11.5, 208.5) is slightly larger at 10 % loss: +0.29 % (SPB) and +0.96 %
  (uniform), vs +0.26 % and +0.69 % at 5 %. So the recommendation holds for any loss between 5 and 10 %.

## Step 12 — power split at the new set point, full acceptance mean (2026-11-01)

(12.5, 207.5) MHz, 3 cm beam, 5 % loss; full grid (44 points) as in Step 2.

| F'=3 + F'=4 (mW) | 25 + 35 | 30 + 30 | **35 + 25** | 40 + 20 |
|---|---|---|---|---|
| SPB, on axis | 99.75 | 99.85 | 99.88 | 99.88 |
| SPB, acceptance mean | 99.34 | 99.44 | **99.50** | **99.50** |
| uniform, on axis | 99.64 | 99.64 | 99.52 | 99.22 |
| uniform, acceptance mean | **97.82** | 97.66 | 97.26 | 96.53 |

- Same picture as the proxy (Step 7): the SPB start is flat from 35 + 25 to 40 + 20 mW and loses only 0.06 % at
  30 + 30; the uniform start wants more F'=4 power (+0.40 % at 30 + 30, +0.56 % at 25 + 35).
- Which split is best depends on how much of the J=2 population at CU2 is outside m_J = 0. With a fraction u of
  the population distributed like the uniform start, 30 + 30 beats 35 + 25 for u ≳ 13 % and 25 + 35 beats it for
  u ≳ 22 %. **Keep 35 + 25 mW if the lens really delivers mostly m_J = 0; 30 + 30 mW is the safer choice** if the
  m_J ≠ 0 fraction is not known (it costs at most 0.06 %).
- The set point was optimised at 35 + 25; whether 30 + 30 or 25 + 35 prefers a different (δ0, f_AOM) is the next
  run.

## Step 13 — does the set point depend on the split? (2026-11-01)

Proxy (%), 3 cm beam, 5 % loss; rows δ0/2π (MHz), columns f_AOM (MHz). 35 + 25 mW from Steps 4–5 for comparison.

| δ0 \ f_AOM | 206.0 | 207.5 | 208.5 |
|---|---|---|---|
| **SPB, 30 + 30** | | | |
| 11.5 | 99.28 | 99.33 | 99.34 |
| 12.0 | 99.46 | 99.51 | 99.45 |
| 12.5 | 99.57 | **99.58** | 99.48 |
| 13.0 | 99.50 | 99.42 | 99.32 |
| **SPB, 25 + 35** | | | |
| 11.5 | 99.21 | 99.25 | 99.26 |
| 12.0 | 99.38 | 99.42 | 99.37 |
| 12.5 | 99.45 | **99.46** | 99.36 |
| 13.0 | 99.33 | 99.24 | 99.12 |
| **uniform, 30 + 30** | | | |
| 11.5 | 97.14 | 97.61 | 97.61 |
| 12.0 | 97.73 | 98.06 | 97.64 |
| 12.5 | 98.17 | 98.24 | 97.74 |
| 13.0 | **98.45** | 97.97 | 97.63 |
| **uniform, 25 + 35** | | | |
| 11.5 | 97.38 | 97.77 | 97.79 |
| 12.0 | 97.94 | 98.18 | 97.85 |
| 12.5 | 98.35 | 98.37 | 97.90 |
| 13.0 | **98.53** | 98.10 | 97.72 |

- **The set point does not depend on the split.** For every split from 25 + 35 to 35 + 25 mW the SPB start peaks
  at (12.5, 207.5) and the uniform start at (13, 206), the same as in Steps 4–6. So the split can be chosen
  independently of the set point, and the comparison in Step 12 (all at (12.5, 207.5)) is fair.
- At (13, 206) the uniform-start gain of 30 + 30 over 35 + 25 is +0.32 % (98.45 vs 98.13).

## Step 14 — X/Z switching frequency (2026-11-01)

(12.5, 207.5) MHz, 35 + 25 mW, 3 cm beam, 5 % loss, proxy (%). Switching phase 0 at t = 0 as in every run so far.

| Switching frequency (MHz) | 0.1 | 0.25 | **0.5** | 1 | 2 |
|---|---|---|---|---|---|
| SPB | 99.70 | 99.70 | **99.64** | 99.73 | 99.67 |
| uniform | 97.51 | 97.46 | **97.87** | 97.63 | 97.74 |

- The 0.5 MHz values reproduce Steps 7 and 10 exactly, so the new runtime parameter is wired correctly.
- **The switching frequency hardly matters**: ±0.05 % (SPB) and ±0.2 % (uniform) over a factor of 20, with no
  trend. 0.5 MHz is the best value for the uniform start and within 0.1 % of the best for the SPB start. No reason
  to change it.
- The scatter without a trend suggests it comes from the switching phase at the molecule's arrival, which is
  fixed (0) in the model but random in the lab. With σz = 0.75 mm one pass lasts ~10 µs (intensity FWHM), i.e. only ~2.4 periods
  at 0.25 MHz. Next run: average over the phase.

## Step 15 — X/Z switching phase (2026-11-01)

(12.5, 207.5) MHz, 35 + 25 mW, proxy (%). Phase of the X/Z square wave at t = 0 (molecule at z = 0, 3.5 mm
before the first pass), which is random in the lab.

| Frequency | phase 0 | π/2 | π | 3π/2 | mean |
|---|---|---|---|---|---|
| SPB, 0.25 MHz | 99.70 | 99.71 | 99.70 | 99.71 | 99.71 |
| SPB, 0.5 MHz | 99.64 | 99.65 | 99.64 | 99.63 | 99.64 |
| SPB, 1 MHz | 99.73 | 99.73 | 99.73 | 99.73 | 99.73 |
| uniform, 0.25 MHz | 97.46 | 97.48 | 97.46 | 97.48 | 97.47 |
| uniform, 0.5 MHz | 97.87 | 97.90 | 97.87 | 97.87 | 97.88 |
| uniform, 1 MHz | 97.63 | 97.61 | 97.63 | 97.62 | 97.62 |

- **The phase does not matter** (≤ 0.03 %), so the fixed phase 0 used in every run is not a bias.
- The small frequency dependence of Step 14 is therefore real, not a phase artefact, but it is ≤ 0.1 % (SPB) and
  ≤ 0.4 % (uniform) and 0.5 MHz is already the best value for the uniform start. Keep 0.5 MHz.

Study stopped here on 2026-11-01 (user request). Summary: `docs/state_cleanup/cu2_power_summary.md`.

Resumed 2026-11-02 for the question "what changes if the CU2 field is 155 V/cm instead of 200 V/cm?". On
2026-11-02 the user also measured that the total UV power can reach **75 mW** (more green input), so from Step 19
on 75 mW is the operating power.

## Step 16 — 155 V/cm: line shifts and the old set point (2026-11-02)

Model rebuilt at 155 V/cm (`e155` option of `scripts/cu2_power_check.py`; every file name carries `_e155`).
Δ_main = 204.520 MHz (211.829 MHz at 200 V/cm). Shifts of the R(2) lines from 200 to 155 V/cm (`tmp/cu2_stark155.py`):

| Line | Shift |
|---|---|
| F'=3 main (X F=3 mF=2 → F'=3 mF'=3) | +9.80 MHz |
| F'=4 stretched (X F=3 mF=3 → F'=4 mF'=4) | +9.51 MHz |
| SPB level X F1=5/2 F=2 mF=0 → F'=3 | +3.45 MHz |
| SPB level X F1=3/2 F=1 mF=0 → F'=3 | +3.35 MHz |
| F'=4 main (X F=3 mF=0 → F'=4 mF'=1) | +2.49 MHz |

The SPB-level lines move 6.4 MHz less than the F'=3 main line, so the set point cannot simply be translated. At
155 V/cm the four SPB levels are 0.80–0.86 m_J = 0 (0.87–0.92 at 200 V/cm).

(12.5, 207.5) MHz, 35 + 25 mW, proxy (%):

| Case | SPB | uniform |
|---|---|---|
| 200 V/cm (reference) | 99.64 | 97.87 |
| 155 V/cm, laser and AOM unchanged | 99.06 | 89.55 |
| 155 V/cm, laser on the moved F'=3 main line, AOM keeps δ1 (f_AOM 200.2 MHz) | 99.42 | 93.21 |

## Step 17 — 155 V/cm: set point at 60 mW (2026-11-02)

35 + 25 mW, proxy (%). Rows δ0/2π, columns f_AOM (MHz). δ0 is relative to the 155 V/cm F'=3 main line.

Coarse grid `opt155`:

| SPB | 192 | 196 | 200 | 204 | 208 | 212 |
|---|---|---|---|---|---|---|
| 3 | 99.19 | 99.05 | 99.35 | 99.45 | 99.23 | 99.26 |
| 6 | 99.10 | 99.23 | 99.46 | 99.21 | 99.32 | 99.12 |
| 9 | 98.65 | 99.19 | 98.85 | 99.03 | 98.70 | 99.14 |
| 12.5 | 99.49 | 99.31 | 99.40 | 99.12 | 99.51 | 99.49 |
| 16 | 98.91 | 99.01 | 98.64 | 98.95 | 99.16 | 98.58 |

| uniform | 192 | 196 | 200 | 204 | 208 | 212 |
|---|---|---|---|---|---|---|
| 3 | 91.34 | 88.08 | 91.61 | 93.95 | 91.68 | 92.44 |
| 6 | 89.37 | 88.92 | 93.55 | 90.28 | 93.21 | 87.49 |
| 9 | 84.75 | 93.24 | 89.41 | 93.42 | 87.43 | 92.80 |
| 12.5 | 94.24 | 92.87 | 94.73 | 90.98 | 94.24 | 94.63 |
| 16 | 90.88 | 92.19 | 89.10 | 90.08 | 93.06 | 87.02 |

Fine grid `fine155`:

| SPB | 198 | 200 | 202 | 208 | 210 | 212 |
|---|---|---|---|---|---|---|
| 11.5 | 99.24 | 99.40 | 99.37 | 99.28 | 99.49 | 99.46 |
| 12.5 | 99.26 | 99.40 | 99.24 | 99.51 | 99.45 | 99.49 |
| 13.5 | 99.32 | 99.26 | 99.06 | 99.37 | 99.39 | 99.22 |

| uniform | 198 | 200 | 202 | 208 | 210 | 212 |
|---|---|---|---|---|---|---|
| 11.5 | 91.94 | 94.09 | 94.52 | 91.65 | 94.29 | 94.30 |
| 12.5 | 92.27 | 94.73 | 93.02 | 94.24 | 93.51 | 94.63 |
| 13.5 | 93.79 | 94.32 | 92.24 | 94.09 | 94.34 | 91.75 |

- δ0 ≈ 12.5 MHz is best for both starts, as at 200 V/cm. Best uniform start 94.73 % (12.5, 200), ~3 % below
  200 V/cm; the SPB start loses only ~0.1 %.
- At 60 mW the landscape is rough along f_AOM: 2 MHz steps change the uniform start by 1.5–2.5 % without a trend,
  unlike the smooth fine grids at 200 V/cm (Steps 4–5). The repeated point (12.5, 200) is identical in both grids,
  so it is not numerical noise. A guess (not checked): two-photon resonances between the F'=3 and F'=4 beams and
  their EOM sidebands, which are Doppler-free for co-propagating beams and so not averaged by the vx spread.

## Step 18 — reduced acceptance grid (2026-11-02, offline)

Tested on the 28 full acceptance grids (11 vx × 4 y0) already stored: keeping every second vx (−1, −0.6, −0.2,
0.2, 0.6, 1 m/s) and all four y0 (24 points, `VX24`) reproduces the 44-point mean to ≤ 0.08 % (rms 0.03 %).
Dropping any y0 row costs 0.13–0.7 % (the y0 interpolation is logarithmic and the 6 mm row matters). Stages
`acc24` (both starts) and `perstate` use this grid.

## Step 19 — 155 V/cm, 75 mW: power split (2026-11-02)

Proxy (%), 155 V/cm.

| F'=3 + F'=4 (mW) | SPB (12.5, 212) | SPB (12.5, 200) | uniform (12.5, 212) | uniform (12.5, 200) |
|---|---|---|---|---|
| 35 + 25 (60 mW) | 99.49 | 99.40 | 94.63 | 94.73 |
| 44 + 31 | 99.67 | 99.61 | 96.49 | 96.54 |
| 38 + 37 | 99.70 | 99.64 | 97.17 | 97.18 |
| **32 + 43** | **99.70** | **99.64** | 97.61 | 97.59 |
| 26 + 49 | 99.66 | 99.60 | 97.85 | 97.79 |
| 20 + 55 | 99.56 | 99.50 | 97.87 | 97.76 |

- 75 mW gains ~0.2 % (SPB) and ~3 % (uniform) over 60 mW at 155 V/cm.
- **32 + 43 mW** keeps the SPB start at its maximum; 26 + 49 is better only if more than ~14 % of the population
  behaves like the uniform start (the same crossover as at 200 V/cm).

## Step 20 — 155 V/cm, 75 mW: set point cuts (2026-11-02)

32 + 43 mW, proxy (%). Along f_AOM at δ0 = 12.5 MHz:

| f_AOM | 197 | 198 | 199 | **200** | 201 | 202 | 203 |
|---|---|---|---|---|---|---|---|
| SPB | 99.46 | 99.50 | 99.62 | **99.64** | 99.60 | 99.49 | 99.40 |
| uniform | 95.18 | 95.78 | 96.99 | **97.59** | 97.04 | 96.14 | 95.70 |

| f_AOM | 207 | 208 | 209 | 210 | 211 | **212** | 213 | 214 |
|---|---|---|---|---|---|---|---|---|
| SPB | 99.47 | 99.68 | 99.69 | 99.60 | 99.65 | **99.70** | 99.59 | 99.42 |
| uniform | 94.89 | 97.26 | 97.44 | 96.35 | 97.36 | **97.61** | 95.73 | 93.91 |

Along δ0 (laser scan, AOM fixed):

| δ0 | 11.5 | 12 | 12.5 | 13 | 13.5 |
|---|---|---|---|---|---|
| f_AOM 200, SPB | 99.58 | 99.67 | 99.64 | 99.55 | 99.49 |
| f_AOM 200, uniform | 96.88 | 97.47 | 97.59 | 97.27 | 96.88 |
| f_AOM 212, SPB | 99.63 | 99.71 | 99.70 | 99.59 | 99.44 |
| f_AOM 212, uniform | 97.30 | 97.89 | 97.61 | 96.51 | 95.13 |

- At 75 mW the f_AOM structure is much smoother than at 60 mW (more power broadening).
- Region A (f_AOM ≈ 200 MHz): one clean peak, −0.6 % (uniform) for ±1 MHz in f_AOM, −0.3…−0.6 % for ±1 MHz in
  δ0. Its F'=4 detuning δ1 = 7.98 MHz equals that of the 200 V/cm set point (8.17 MHz): it is the same
  configuration relative to the two main lines.
- Region B (f_AOM ≈ 208–212 MHz): slightly higher peak but a dip at 210 MHz and a steep edge at 213 MHz
  (−1.9 %); δ0 +1 MHz costs 1.4 % (uniform). 212 MHz is also near the top of the AOM range (~215 MHz).

## Step 21 — 155 V/cm, 75 mW: candidates on the reduced grid, and 200 V/cm at 75 mW (2026-11-03)

Acceptance mean on the 24-point grid (Step 18), 32 + 43 mW, 155 V/cm:

| Set point (δ0, f_AOM) | SPB | uniform |
|---|---|---|
| B (12, 212) | **99.66** | **97.30** |
| A (12, 200) | 99.58 | 96.91 |
| A (12.5, 200) | 99.57 | 97.01 |

200 V/cm, 32 + 43 mW, proxy (%):

| Set point | SPB | uniform |
|---|---|---|
| (12.5, 207.5) | 99.73 | 99.23 |
| (12, 207.5) | 99.65 | 99.07 |
| (12.5, 206) | 99.72 | 99.25 |
| (13, 206) | 99.66 | 99.31 |

- At 200 V/cm, 75 mW raises the uniform start from 97.87 % to 99.23 % (proxy) at the 60 mW set point, with no
  re-tuning.
- At 155 V/cm the best proxy point (12, 212) gives 99.71 / 97.89 %: the SPB start is the same as at 200 V/cm,
  the uniform start ~1.3 % lower. If the field can be chosen, 200 V/cm is better for any m_J ≠ 0 population.

## Step 22 — 155 V/cm, 75 mW: field drift with laser and AOM fixed (2026-11-03)

Laser and AOM held at the 155 V/cm set point while the field moves by ±3 V/cm, so δ0 → δ0 − (line shift) and
δ0 + f_AOM keeps its relation to the laser (shift −0.651 MHz at 158 V/cm, +0.651 MHz at 152 V/cm). Proxy, 32 + 43 mW
(SPB / uniform, %):

| Set point | 152 V/cm | 155 V/cm | 158 V/cm |
|---|---|---|---|
| A (12.5, 200) | 99.63 / 97.48 | 99.64 / 97.59 | 99.65 / 97.64 |
| B (12, 212) | 99.71 / 97.60 | 99.71 / 97.89 | 99.72 / 97.96 |

- Both set points tolerate ±3 V/cm: ≤ 0.01 % (SPB) and ≤ 0.3 % (uniform). B is better than A at all three fields.

## Step 23 — 155 V/cm: 60 mW reference at the 75 mW set point (2026-11-03)

Reduced-grid acceptance mean at B (12, 212), 155 V/cm (%):

| Power | SPB | uniform |
|---|---|---|
| 35 + 25 mW (60 mW) | 99.45 | 94.22 |
| 32 + 43 mW (75 mW) | 99.66 | 97.30 |

- The extra 15 mW, put into the F'=4 beam, gains +0.21 % (SPB) and +3.1 % (uniform). At 155 V/cm the uniform start
  needs the power much more than at 200 V/cm (the 60 mW optimum there was 97.26 % on the full grid).
- (12, 212) is not the 60 mW optimum at 155 V/cm (Step 17: SPB 99.51 % at (12.5, 208), proxy), so this compares
  powers at fixed tuning, not optimum against optimum.

## Step 24 — 155 V/cm, 75 mW: efficiency of each J=2 level (2026-11-03)

Each of the 20 J=2 levels alone as the initial state, at B (12, 212), 32 + 43 mW, 155 V/cm, on the reduced
acceptance grid (288 trajectories for the 12 levels with m_F ≥ 0). The OBEs are linear in ρ0, so the efficiency of
any incoherent J=2 distribution p_i is Σ p_i ε_i. Levels with m_F < 0 are copied from their m_F > 0 partner; a
direct check of all 8 partners at vx = ±0.6 m/s, y0 = 0 agrees to ≤ 1.2 × 10⁻⁵ of J=2. Table with the full vx × y0
grid: `docs/state_cleanup/data/cu2_perstate_155Vcm_75mW_d12_f212.csv`.

Levels labelled by their zero-field parent; P(m_J=0) at 155 V/cm (%):

| Level (F1, F, \|m_F\|) | P(m_J=0) | on axis (vx 0.2, y0 0) | acceptance mean |
|---|---|---|---|
| 3/2, 1, 0 | 0.20 | 97.28 | 92.02 |
| 3/2, 1, 1 | 0.00 | 99.59 | 93.58 |
| 3/2, 2, 0 | 0.14 | 99.62 | 96.62 |
| 3/2, 2, 1 | 0.19 | 98.13 | 93.79 |
| 3/2, 2, 2 | 0.00 | 97.03 | 92.61 |
| 5/2, 2, 0 (SPB, weight 2/3) | 0.86 | 99.96 | 99.66 |
| 5/2, 2, 1 | 0.00 | 99.97 | 99.75 |
| 5/2, 2, 2 | 0.00 | 99.98 | 99.84 |
| 5/2, 3, 0 (SPB, 1/9) | 0.80 | 99.85 | 99.54 |
| 5/2, 3, 1 (SPB, 1/9 each) | 0.80 | 99.95 | 99.71 |
| 5/2, 3, 2 | 0.00 | 99.96 | 99.81 |
| 5/2, 3, 3 | 0.00 | 99.95 | 99.86 |

- Weighted sums: SPB 99.65 %, uniform 97.29 %, against 99.66 / 97.30 % from the direct runs (Step 21).
- All 12 F1 = 5/2 levels are ≥ 99.5 %. The whole uniform-start residual comes from the 8 F1 = 3/2 levels
  (mean 93.6 %): they are cleaned well on axis but fall to 52–80 % at |vx| = 1 m/s, y0 = 9 mm, i.e. their lines are
  narrow (weakly coupled, little power broadening). A start with fraction q in F1 = 3/2 has efficiency
  ≈ 99.77 − 6.2 q % (q = 0.4 for the uniform start).
- The stretched F1 = 5/2, F = 3, |m_F| = 3 levels, which limited the uniform start at 200 V/cm / 60 mW, are at
  99.86 % here.
- Coherences between levels (e.g. within a ±m_F pair) are not covered by the weighted sum.

## Step 25 — 200 V/cm, 75 mW: candidates on the reduced grid (2026-11-03)

Reduced-grid acceptance mean, 32 + 43 mW, against the 155 V/cm best point (%):

| Field, set point | SPB | uniform |
|---|---|---|
| 200 V/cm, (12.5, 207.5) | 99.61 | 98.81 |
| 200 V/cm, (12.5, 206) | 99.60 | **98.91** |
| 155 V/cm, (12, 212) | **99.66** | 97.30 |

- On the same grid and power, 200 V/cm is 1.5–1.6 % better for the uniform start and 0.05 % worse for SPB. The
  200 V/cm points were optimised at 60 mW, and the 32 + 43 split at 155 V/cm, so the 200 V/cm numbers could still
  improve with a 75 mW re-optimisation.
- **The CU2 field is now fixed at 155 V/cm (user, 2026-11-03).** The 200 V/cm per-state run at (12.5, 206) was
  stopped at 12 of 18 chunks and is not used; Steps 21 and 25 at 200 V/cm stay as history.

## Step 26 — 155 V/cm, 75 mW: laser scan with the AOM fixed (2026-11-03)

The AOM was the reason to leave B (12, 212): its band is 200 ± 15 MHz, so 212 MHz is near the edge. With f_AOM fixed
at the band centre, 200 MHz, a laser scan moves both beams together, δ1 = δ0 + f_AOM − Δ_main = δ0 − 4.52 MHz (a 45°
line in the (δ0, δ1) plane). Stage `corr` of `scripts/cu2_power_check.py`, on axis (vx = 0, y0 = 0), δ0 = −60 … 140
MHz in 1 MHz steps, 32 + 43 mW. Figure: `docs/state_cleanup/figures/CU2_corr_scan_e155_f200.{pdf,png}`
(`scripts/plot_cu2_corr_scan.py 200 32 43 e155 12.5`).

| | SPB | uniform |
|---|---|---|
| highest peaks | 5 MHz 99.994 %, 22 MHz 99.982 %, 40 MHz 99.966 % | 5 MHz 99.925 %, 22 MHz 99.889 %, 40 MHz 99.557 % |
| set point δ0 = 12.5 MHz | 99.84 % | 98.88 % |
| ≥ 99 % | −13 … 57 MHz | only at the peaks (5–6, 13, 22–23, 40 MHz) |
| ≥ 98 % | −14 … 63 MHz | 4 … 40 MHz, not continuous |

- The peaks repeat every 17–18 MHz, the EOM frequency. Each uniform-start peak is only 2–3 MHz wide above 98 %;
  between them the efficiency falls to 80–91 %.
- The SPB start has a broad plateau: a laser offset of a few MHz costs little.
- Above 100 MHz both starts fall to a few %, unlike the old 250 V/cm, 80 + 24 mW scan (δ1 − δ0 = −25 MHz there).

## Step 27 — 155 V/cm, 75 mW: do the on-axis peaks survive the acceptance? (2026-11-03)

Proxy (Step 4) at f_AOM = 200 MHz, 32 + 43 mW (%):

| δ0 (MHz) | on axis SPB | on axis uniform | proxy SPB | proxy uniform |
|---|---|---|---|---|
| 5 | 99.99 | 99.92 | 99.44 | 95.54 |
| **12.5** | 99.84 | 98.88 | **99.64** | **97.59** |
| 13 | 99.87 | 99.07 | 99.55 | 97.27 |
| 22 | 99.98 | 99.89 | 99.48 | 95.84 |
| 23 | 99.96 | 99.81 | 99.55 | 96.01 |

- The narrow on-axis peaks lose over the acceptance: |vx| = 0.6 m/s shifts the lines by 2.2 MHz, into the valleys.
  12–12.5 MHz (Step 20: 12 MHz gives 99.67 / 97.47 %) stays the best. **Set point at 155 V/cm: δ0 = 12.5 MHz,
  f_AOM = 200 MHz, 32 + 43 mW**; reduced-grid acceptance mean 99.57 % (SPB) / 97.01 % (uniform), Step 21.

## Step 28 — final population of every J=2 level for the CU2 test start (2026-11-03)

For the first CU2 test (CU1 and SPB off) the J=2 molecules reaching CU2 are the four lens levels (P(m_J=0) > 0.5:
F1=5/2, F=2, m_F=0 and F=3, m_F=0, ±1). Stage `cu2final` of `scripts/cu2_power_check.py` runs each of them alone
(3 by ±m_F reflection) on the reduced grid at (12.5, 200) MHz and saves the final population of every level
(`tmp/cu2_power/vz184/cu2final_loss5_d12p5_f200_32_43_e155_result.npz`, 72 trajectories, 184 m/s).

| start level | J=2 left (acceptance mean) |
|---|---|
| F=2, m_F=0 (singlet) | 0.49 % |
| F=3, m_F=0 | 0.74 % |
| F=3, m_F=±1 | 0.38 % |

- For either lens weighting, 0.49–0.50 % of J=2 is left, **94 % of it in F1 = 3/2** and only 3 % in F=3: decay
  from B feeds F1 = 3/2, which CU2 pumps weakly (Step 24).

## Step 29 — detection model for the CU2 test (2026-11-03)

`scripts/cu2_detection_test.py`. Probe (user): DET2, locked to the main line, no EOM, flat top 1.5 cm × 1.5 cm,
15 mW, single pass, same direction as CU2; E = 0. The DET2 light is σ+ about its direction k; the DET field plates
(E along k) are off for the CU2 test, so the axis is the earth field, taken as 0.5 G vertical (⊥ k, the CU2 field
axis). Molecules arrive in the zero-field parent (F1, F, m_F) of their CU2 eigenstate (adiabatic ramp, no
coherences); acceptance |vx| ≤ 1 m/s and |y| ≤ 7.5 mm (probe height). Output: photons per molecule per J=2 level,
including leakage of the excited level to other rotational states.

Mixed-state |⟨e|d|g⟩| at zero field from each J=2, m_F=0 level (X light, `tmp/q2_me_table.py`):

| line | F1=3/2 F=1 | F1=3/2 F=2 | F1=5/2 F=2 | F1=5/2 F=3 |
|---|---|---|---|---|
| R(2) F1'=7/2 F'=4 | 0 | 0 | 0 | 0.31 |
| R(2) F1'=5/2 F'=3 | 0 | 0.35 | 0.004 | 0.08 |
| Q(2) F1'=3/2 F'=1 | 0.41 | 0.08 | 0.009 | 0 |
| Q(2) F1'=3/2 F'=2 | 0.14 | 0.43 | 0.02 | 0.005 |
| Q(2) F1'=5/2 F'=2 | 0.009 | 0.03 | 0.47 | 0.06 |
| Q(2) F1'=5/2 F'=3 | 0 | 0.003 | 0.09 | 0.48 |

- The Q(2) lines are nearly diagonal: each sees mainly the ground level with the same (F1, F). Lines are ≥ 230 MHz
  apart (B J'=2: F1'=5/2 F'=2 at 0, F'=3 at 230, F1'=3/2 F'=1 at 667, F'=2 at 1249 MHz), so only the locked line is
  modelled.
- R(2) F1'=5/2 F'=3 does not see F1=5/2 F=2 although ΔF = +1: the R-branch ΔF1 = 0 amplitude (+0.035) is cancelled
  by the 13 % J'=2 admixture of the B level (−0.031, Tl hyperfine h1 = 28.8 GHz) (`tmp/why_f2_dark.py`).

> **Superseded by Step 32.** The photon numbers below and the readings in Step 30 used centrex in the vertical basis,
> where its pairwise decay operators are wrong for σ+ light about k. The matrix elements above are unaffected.

Photons per molecule at vx = 0, 15 mW, σ+ (earth field ⊥ k):

| line | F1=3/2 levels | F1=5/2 F=2 | F1=5/2 F=3 |
|---|---|---|---|
| R(2) F1'=7/2 F'=4 | 0 | 0.01 | 2.11–2.13 (saturated: same at 7.5 and 30 mW) |
| Q(2) F1'=5/2 F'=3 | 0.02 | 1.5–2.0 | 2.9–5.2 |
| Q(2) F1'=3/2 F'=2 | 1.7–2.7 | 0.2–0.4 | 0.04–0.07 |

- R(2) F'=4 (F → F+1) has no dark states for any polarization. The Q lines (ΔF = 0) do: with linear light ∥ B every
  m_F = 0 level is dark, and with pure σ+ (field along k, only with the DET plates on) the stretched m_F = +F level.
  With the plates off the light is half π, half σ± about the axis, so there is no full dark state.
- The 0.01–0.02 photons of F1=5/2 F=2 on F'=4 come from F=2/F=3 mixing by the 0.5 G field.
- At |vx| = 1 m/s (Doppler 3.7 MHz) the yields drop to 1.0–1.3 (R(2) F'=4) and 0.9–1.8 (Q(2) F1'=3/2 F'=2).

## Step 30 — expected CU2 on/off readings, recommended detection lines (2026-11-03)

> **Numbers superseded by Step 32** (decay-operator issue); the choice of lines and their roles stand.

Mean photons per J=2 molecule over the acceptance, CU2 at (12.5, 200) MHz, σ+ probe (earth field ⊥ k), for the lens
weightings equal (1/4 each), singlet 2/3 / triplets 1/9 each, and triplet only (1/3 each):

| line (role) | CU2 off | CU2 on | on/off |
|---|---|---|---|
| R(2) F1'=7/2 F'=4 (check) | 1.26 / 0.57 / 1.67 | 2.3 / 2.1 / 2.5 × 10⁻⁴ | 1.8 / 3.7 / 1.5 × 10⁻⁴ |
| Q(2) F1'=5/2 F'=3 (normalization) | 2.43 / 1.89 / 2.75 | ~5 × 10⁻⁴ | 0.02–0.03 % |
| Q(2) F1'=3/2 F'=2 (residual) | 0.04 / 0.06 / 0.02 | 6.8 / 6.5 / 7.0 × 10⁻³ | large (CU2-off signal tiny) |

True J=2 left 0.49–0.50 %, F1=3/2 left 0.46–0.48 %.

- **R(2) F'=4 sees only F=3**, where 3 % of the residual is: it reads 13–27× below the true residual. It checks that
  CU2 empties F=3 (frequency, power), and it is saturated and free of dark states.
- **Q(2) F1'=5/2 F'=3** sees all four lens levels: the CU2-off signal counts the J=2 molecules. Its CU2-off signal
  over that of R(2) F'=4 measures the lens weighting (0.52 equal, 0.30 for 2/3, 1/9).
- **Q(2) F1'=3/2 F'=2** sees all 8 F1=3/2 levels and hardly the lens levels: its CU2-on signal is the residual,
  1.4–1.45 photons per residual molecule, nearly background-free.
- Lines 2 and 3 together: S_on[Q(2) F1'=3/2 F'=2] / S_off[Q(2) F1'=5/2 F'=3] = 0.28 / 0.34 / 0.25 %, i.e. the
  F1=3/2 residual times 1/1.68, 1/1.34, 1/1.87. With a field component along k (pure σ+, mF-averaged bound) the
  factors are 1.66, 1.16, 1.97. One model factor 1.2–2.0 covers every weighting and field direction.
- For comparison, linear X light ⊥ B gives R(2) F'=4 1.9 / 3.9 × 10⁻⁴ and factors 1.87 / 1.53 (equal / 2/3).
- **RC on P(2) F'=2** (F'=2 decays only to J=0 F=1, so the lens start becomes triplet-rich): R(2) F'=4 CU2-off ×3
  (0.57 → 1.67), CU2-on only +17 %, on/off 3.7 → 1.5 × 10⁻⁴ and nearly independent of the weighting (9 % between 0
  and 10 % singlet). CU2 efficiency unchanged (0.49 → 0.50 % left).
- No bias magnetic field can be added in the DET region (user), so the earth-field direction is taken as given.

Schematic (slides): `docs/state_cleanup/figures/slides/CU2_test_onoff_obe_e155_d12p5_f200.{png,pdf}`
(`python scripts/cu2_detection_test.py plot 12.5,200 32 43 e155`). Test plan: `docs/state_cleanup/CU2_test_plan.md`.

## Step 31 — 155 V/cm, 75 mW: vx × y0 map at the set point (2026-11-03)

Stage `map` of `scripts/cu2_power_check.py`: 0.2 m/s × 3 mm grid, vx −1 … 1 m/s, y0 0 … 15 mm (mirrored), at
(12.5, 200) MHz, 32 + 43 mW, both starts; the 24 reduced-grid points are reused, 42 new per start. Figure
`docs/state_cleanup/figures/cleanup_map_vx_y_CU2_155Vcm.{pdf,png}` (`scripts/plot_cleanup_maps_current.py cu2_155`).

| | on axis | full-grid mean in acceptance | reduced grid (Step 21) |
|---|---|---|---|
| SPB | 99.88 % | 99.56 % | 99.57 % |
| uniform | 99.18 % | 96.92 % | 97.01 % |

- The uniform start has a shallow dip exactly at vx = 0 (99.18 % against 99.48 / 99.51 % at vx = ∓0.2 m/s); the SPB
  start a smaller one (99.88 against 99.93 / 99.94 %). Both are best at |vx| ≈ 0.2 m/s.
- Along vx the uniform start falls to 91 % at |vx| = 1 m/s (y0 = 0), the SPB start to 98.7 %. Along y0 both stay
  near their on-axis values to ~3 mm and fall to 86 % (uniform) and 98 % (SPB) at 9 mm, vx = 0.2 m/s.

## Step 32 — CU2 test detection recomputed with the axis along k (2026-11-04)

**Why.** centrex builds spontaneous decay from one single-jump operator √(Γ·BR) |g⟩⟨e| per (excited, ground) pair
(`centrex_tlf/couplings/collapse.py`), and its prepared solvers accept nothing else
(`lindblad/plan_static.py::_lower_structured_jumps`). Decay therefore never transfers excited coherences into
ground coherences, and the result depends on the quantization axis. Test: an isotropic X J=2 F1=5/2 F=3 ensemble on
Q(2) F1'=5/2 F'=3 (15 mW, vx = 0, B = 1 mG) must give the same photon number for every polarization:

| light | centrex (pairwise decay) | coherent decay C_q (numpy, dense expm) |
|---|---|---|
| linear z | 4.75 | 4.66 |
| linear x | 4.13 | 4.66 |
| σ+ about k = y (vertical axis, used in Steps 29–30) | 4.01 | 2.77 |
| σ+ about z (axis along k) | 2.72 | 2.77 |

The numpy integrator reproduces centrex exactly with the pairwise operators. With the axis along the light's symmetry
axis (z for linear z, k for σ+) pairwise decay is within ~2 % of coherent decay; otherwise per-level errors reach
~75 % (e.g. vertical m_F = ±3 on Q(2) F1'=5/2 F'=3: 5.16 instead of 2.97).

**Method (default `pol k` in `scripts/cu2_detection_test.py`).** Axis along k, where the DET2 light is pure σ+
(centrex run `spk`, B along k). The CU2 populations (m_F about the vertical CU2 field) are rotated onto k:
n(F1, F, m_F) = Σ_m' |d^F_{m'm}(π/2)|² n_k(F1, F, m'). The rotation creates coherences Δm' = ±2, ±4, ±6, which do
not contribute while the probe is symmetric about k; the 0.5 G earth field (⊥ k, plates off) turns a molecule by
≤ 0.4 rad in the 81.5 µs probe and is neglected. Example: the three F=3 triplets (m_F = 0, ±1, 1/3 each) map onto
k as m' = 0, ±1, ±2, ±3 → 0.125, 0.073, 0.104, 0.260 each (⟨J_k²⟩ = 17/3, against ⟨J_z²⟩ = 2/3 vertically).
Checks: `pol spc` = exact (vertical axis, earth field 0.5 G vertical ⊥ k, coherent decay; `yields_coherent`),
`pol kiso` = alignment scrambled between CU2 and the probe (m'-averaged).

Photons per molecule at vx = 0, 15 mW (m_F about the vertical):

| line | F1=3/2 levels | F1=5/2 F=2 | F1=5/2 F=3 |
|---|---|---|---|
| R(2) F1'=7/2 F'=4 | 0 | 0.01–0.02 | 2.12–2.13 |
| Q(2) F1'=5/2 F'=3 | 0.01 | 0.92–1.05 | 2.54–2.86 |
| Q(2) F1'=3/2 F'=2 | 1.45–1.91 | 0.22–0.33 | 0.04–0.08 |

Readings over the acceptance, CU2 at (12.5, 200) MHz, 32 + 43 mW, lens weighting equal / 2/3 (singlet) / triplet
only; axis along k (exact check in brackets where it differs):

| line (role) | S_off, photons/molecule | S_on, photons/molecule | S_on/S_off |
|---|---|---|---|
| R(2) F'=4 (check) | 1.35 / 0.61 / 1.80 | 2.6 / 2.3 / 2.7 × 10⁻⁴ | 1.9 / 3.8 / 1.5 × 10⁻⁴ |
| Q(2) F1'=5/2 F'=3 (normalization) | 2.02 / 1.41 / 2.39 (2.05 / 1.38 / 2.44) | 4.6 × 10⁻⁴ | 0.02–0.03 % |
| Q(2) F1'=3/2 F'=2 (residual) | 0.04 / 0.06 / 0.03 | 5.8 / 5.7 / 5.9 × 10⁻³ | large |

- R(2) F'=4 is saturated and has no dark states, so it hardly changed (on/off was 1.8 / 3.7 / 1.5 × 10⁻⁴): it reads
  13–33× below the true residual.
- Q(2) F1'=5/2 F'=3: the singlet (F1=5/2 F=2) gives 0.92 photons, not 1.82, so S_off is 13–25 % lower than in
  Step 30. S_off[R(2) F'=4] / S_off[Q(2) F1'=5/2 F'=3] = 0.67 (equal), 0.43 (2/3), 0.75 (triplet only).
- Q(2) F1'=3/2 F'=2: 1.24 photons per residual F1=3/2 molecule (was 1.41–1.45).
- Lines 2 and 3 together: S_on[Q(2) F1'=3/2 F'=2] / S_off[Q(2) F1'=5/2 F'=3] = 0.29 / 0.41 / 0.25 % (axis along k),
  0.28 / 0.41 / 0.24 % (exact), 0.28 / 0.40 / 0.24 % (alignment scrambled), against a true F1=3/2 residual of
  0.47 / 0.46 / 0.48 %: model factor 1.6 / 1.1 / 1.9–2.0. One factor 1.1–2.0 covers every weighting and treatment;
  the R(2) F'=4 / Q(2) F1'=5/2 F'=3 CU2-off ratio pins the weighting and with it the factor.
- RC on P(2) F'=2 (triplet-rich start): R(2) F'=4 S_off ×3 (0.61 → 1.80), S_on +18 %, on/off 3.8 → 1.5 × 10⁻⁴.
- Not yet checked for the same issue: the CU2 and CU1 cleanup OBEs (X/Z light with the axis along E, where ±m_F are
  degenerate and decay-transferred coherences do not dephase).

Schematic updated: `docs/state_cleanup/figures/slides/CU2_test_onoff_obe_e155_d12p5_f200.{png,pdf}`.
