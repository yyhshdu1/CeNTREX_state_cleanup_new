# CU2 at the real laser power — summary (2026-11-01)

Simulation study of CU2 with the laser power actually available (at most 60 mW in total, 35 + 25 mW now) and
5 % window loss per pass, 2026-10-30 to 2026-11-01. Full running log with every table:
`docs/state_cleanup/cu2_power_study.md` (Steps 1–15). Script: `scripts/cu2_power_check.py`.

## Recommendation

| Parameter | Recommended | Before (report, 80 + 24 mW / 10 % loss) |
|---|---|---|
| Power split, F'=3 + F'=4 | **35 + 25 mW** (present); 30 + 30 mW if more than ~13 % of the J=2 population at CU2 is outside m_J = 0 | 80 + 24 mW |
| Laser detuning δ0/2π (F'=3 beam) | **12.5 MHz** | 11.5 MHz (R1) |
| AOM frequency f_AOM | **207.5 MHz** (F'=4 beam detuning unchanged, δ1 = +8.17 MHz) | 208.5 MHz (R1) |
| EOM depth β (17.5 MHz) | **1.5** (unchanged) | 1.5 |
| X/Z switching | **0.5 MHz** (unchanged) | 0.5 MHz |
| Beam height | 3 cm is fine; 2.4 cm would gain +0.08 % (SPB) / +0.42 % (uniform) | 3 cm |

The new set point is R1 with the laser +1 MHz and the AOM −1 MHz, i.e. only the F'=3 beam moves.

## Expected cleanup efficiency

Mean over the detection acceptance (|vx| ≤ 1 m/s with σvx = 0.65 m/s, |y0| ≤ 8.5 mm), 3 cm beam, 13 passes,
200 V/cm, vz = 184 m/s. *SPB start* = the four J=2, m_J = 0 levels after SPB (the realistic case, since the lens
focuses m_J = 0); *uniform start* = all 20 J=2 sublevels (worst case).

| Case | SPB start | uniform start |
|---|---|---|
| Old model: 80 + 24 mW, 10 % loss, R1 | 99.57 % | 96.66 % |
| 35 + 25 mW, 5 % loss, R1 | 99.22 % | 96.66 % |
| **35 + 25 mW, 5 % loss, new set point** | **99.50 %** | **97.26 %** |
| 30 + 30 mW, 5 % loss, new set point | 99.44 % | 97.66 % |
| 25 + 35 mW, 5 % loss, new set point | 99.34 % | 97.82 % |
| 40 + 20 mW, 5 % loss, new set point | 99.50 % | 96.53 % |

With the new set point, 60 mW at 5 % loss does about as well as the old 104 mW at 10 % loss. On the beam axis
(vx = 0, y0 = 0) the new set point gives 99.88 % (SPB) and 99.52 % (uniform).

## Main findings

1. **Power split (Steps 1, 2, 7, 12).** The SPB levels are pumped mainly by the F'=3 beam, while the stretched
   m_F = ±3 levels (which dominate the uniform-start residual) are pumped almost only by F'=4. So the SPB start
   is flat from 35 + 25 to 40 + 20 mW, and the uniform start prefers more F'=4 power. The choice between 35 + 25
   and 30 + 30 mW is worth ≤ 0.06 % (SPB) or 0.40 % (uniform). If a fraction u of the population is like the
   uniform start, 30 + 30 is better for u ≳ 13 % and 25 + 35 is better for u ≳ 22 %.
2. **Set point (Steps 3–6, 13).** R1 was optimised on axis at 80 + 24 mW. At 60 mW the lines are less power
   broadened, so the efficiency drops faster with transverse velocity. Moving δ0 up by 1 MHz gives up a little on
   axis and gains it back at |vx| ≈ 0.4 m/s and beyond: +0.28 % (SPB) and +0.60 % (uniform) over the acceptance.
   (13, 206) MHz is equivalent (99.47 / 97.41 %) and slightly better for wide vx distributions. The new point is
   better than R1 for every vx width from σvx = 0.4 to 1.28 m/s. The optimum is the same for every split from
   25 + 35 to 35 + 25 mW, for beam heights of 2.4 and 3 cm, and for 5 and 10 % loss.
3. **EOM depth and polarization switching (Steps 6, 14, 15).** β = 1.5 is optimal (flat to 1.75). The switching
   frequency changes the result by ≤ 0.1 % (SPB) and ≤ 0.4 % (uniform) between 0.1 and 2 MHz, with no trend;
   0.5 MHz is best for the uniform start. The switching phase at the molecule's arrival does not matter
   (≤ 0.03 %).
4. **Beam height (Step 7).** At fixed power a smaller beam is brighter but loses the edge of the acceptance. The
   optimum is 2.4 cm (SPB 99.58 %, uniform 97.68 %): 2 cm gives 99.53 / 97.49 %, 1.6 cm 99.13 / 96.25 % and
   4 cm 99.24 / 95.78 %. This is tied to the detection acceptance |y0| ≤ 8.5 mm.
5. **Solver accuracy (Step 7).** Tighter tolerances change single points by ≤ 0.04 %.

## Sensitivities (proxy, at the new set point, 35 + 25 mW)

Proxy = 6 points (vx = −0.6, 0, 0.6 m/s × y0 = 0, 6 mm). It reads 0.1–0.7 % higher than the full acceptance mean
but ranks settings in the same order; at the new set point it gives 99.64 % (SPB) and 97.87 % (uniform).

| Change | SPB start | uniform start |
|---|---|---|
| Laser δ0 ± 0.5 MHz | ≤ −0.14 % | ≤ −0.24 % |
| E field ± 2.5 V/cm (laser fixed) | ≤ −0.11 % | ≤ −0.17 % |
| E field ± 5 V/cm (laser fixed) | −0.2 … −1.0 % | ≤ −0.6 % |
| Loss per pass 3 / 8 / 10 / 15 % (instead of 5 %) | +0.12 / −0.25 / −0.48 / −1.23 % | +0.61 / −1.17 / −2.13 / −5.03 % |
| vz = 150 / 220 m/s (instead of 184 m/s) | +0.18 / −0.26 % | +0.98 / −1.22 % |
| Total power 80 mW (47 + 33) / 100 mW (58 + 42) | +0.24 / +0.31 % | +1.26 / +1.77 % |

- **E field.** The lines shift by different amounts (F'=3 main line −2.19 MHz, SPB-level lines −0.78 MHz per
  +10 V/cm), so re-tuning the laser to the main line does not compensate a field error. The CU2 field should be
  known to about ±2 V/cm, or the set point should be found by scanning the laser at the actual field.
- **Loss per pass.** Each extra % of loss costs ~0.08 % (SPB) and ~0.3 % (uniform), comparable to the effect of
  the power split. At 10 % loss the optimal set point is unchanged, and its gain over R1 grows to +0.29 % (SPB) and
  +0.96 % (uniform).

Expected shape of an in-situ laser scan (SPB start, proxy, f_AOM = 207.5 MHz):

| δ0/2π (MHz) | 10 | 11 | 11.5 | 12 | **12.5** | 13 | 13.5 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|
| efficiency (%) | 99.30 | 99.25 | 99.38 | 99.56 | **99.64** | 99.50 | 99.19 | 98.77 | 97.36 |

The maximum is shallow (12–13 MHz within 0.15 %) and falls off faster on the blue side.

## What would most improve the prediction

1. **The m_J ≠ 0 fraction of J=2 at CU2** (from the lens simulation or a measurement). It decides between 35 + 25
   and 30 + 30 mW and whether the true efficiency is near 99.5 % or 97 %.
2. **The real loss per pass.** 5 % vs 10 % is worth 0.5 % (SPB) and 2 % (uniform).
3. **The CU2 field to ±2 V/cm**, or an in-situ laser scan at the operating field.

## Model

56-level OBE (`scripts/cu2_laser_aom_scan.py`). R(2) F1'=7/2 F'=3 and F'=4 from one laser and a single-pass AOM.
200 V/cm, B = 1 mG, one 17.5 MHz EOM on both beams, X/Z switching in antiphase. 13 passes spaced 3.5 mm,
σz = 0.75 mm, 3 cm beam (σy = 7.5 mm), Doppler sign flip between passes. Acceptance mean as in
`scripts/plot_cleanup_maps_current.py` (σvx = 0.65 m/s, σy = 10.35 mm). centrex-tlf 0.3.0, Δ_main = 211.829 MHz.

Code changes:

- `scripts/cu2_power_check.py` contains every stage (`split`, `acceptance`, the set-point grids, `beta`,
  `efield`, `vz`, `loss`, `power`, `polfreq`, `polphase`, `summary`).
- `scripts/cu2_laser_aom_scan.py` gained a runtime parameter `pol_phase` (default 0, so earlier results are
  unchanged).
- Checkpoints are in `tmp/cu2_power/`, which is not tracked by git. `python scripts/cu2_power_check.py summary`
  reprints every table from them.
