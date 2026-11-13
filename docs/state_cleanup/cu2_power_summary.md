# CU2 at the real laser power — summary (updated 2026-11-04)

Simulation study of CU2 with the laser power actually available and 5 % window loss per pass,
2026-10-30 to 2026-11-03. **The CU2 field is 155 V/cm** (confirmed 2026-11-03) and the total UV is **75 mW**
(measured). The AOM band is 200 ± 15 MHz. Full running log with every table: `docs/state_cleanup/cu2_power_study.md`
(Steps 1–32). Scripts: `scripts/cu2_power_check.py` (CU2 OBE stages), `scripts/cu2_detection_test.py` (CU2 test
detection).

Figures:

- `docs/state_cleanup/figures/cleanup_map_vx_y_CU2_155Vcm.{pdf,png}`: vx × y0 map at the set point
  (`scripts/plot_cleanup_maps_current.py cu2_155`).
- `docs/state_cleanup/figures/CU2_corr_scan_e155_f200.{pdf,png}`: laser scan with the AOM fixed
  (`scripts/plot_cu2_corr_scan.py`).
- `docs/state_cleanup/figures/slides/CU2_test_onoff_obe_e155_d12p5_f200.{pdf,png}`: CU2 test on/off schematic.
- `docs/state_cleanup/figures/CU2_power_study.{pdf,png}`: the earlier 200 V/cm, 60 mW study (history).

## Recommendation (155 V/cm, 75 mW)

| Parameter | Recommended | Note |
|---|---|---|
| Power split, F'=3 + F'=4 | **32 + 43 mW** | 26 + 49 mW if more than ~14 % of the J=2 population behaves like the uniform start |
| Laser detuning δ0/2π (F'=3 beam) | **12.5 MHz** | 12–12.5 MHz equivalent |
| AOM frequency f_AOM | **200 MHz** (centre of the AOM band) | F'=4 beam δ1 = δ0 + f_AOM − Δ_main = +7.98 MHz, Δ_main = 204.52 MHz |
| EOM (17.5 MHz), X/Z switching | β = 1.5, 0.5 MHz | unchanged |
| Beam height | 3 cm | unchanged |

(12, 212) MHz is 0.09 % (SPB) / 0.3 % (uniform) better on the acceptance grid, but 212 MHz is near the edge of the
AOM band, where the diffraction efficiency may drop.

## Expected cleanup efficiency

Mean over the detection acceptance (|vx| ≤ 1 m/s with σvx = 0.65 m/s, |y0| ≤ 8.5 mm with σy = 10.35 mm), 3 cm
beam, 13 passes, vz = 184 m/s, reduced 24-point grid (≤ 0.08 % from the full grid). *SPB start* = the four J=2,
m_J ≈ 0 lens levels, weighted 2/3 (singlet) and 1/9 each (triplets); *uniform start* = all 20 J=2 sublevels.

| Case (155 V/cm) | SPB start | uniform start |
|---|---|---|
| **32 + 43 mW, (12.5, 200) MHz** | **99.57 %** | **97.01 %** |
| 32 + 43 mW, (12, 212) MHz | 99.66 % | 97.30 % |
| 35 + 25 mW (60 mW), (12, 212) MHz | 99.45 % | 94.22 % |

The full 0.2 m/s × 3 mm map at the set point (Step 31) gives 99.56 % (SPB) and 96.92 % (uniform) over the
acceptance, and 99.88 % / 99.18 % on axis (vx = 0, y0 = 0). For comparison, 200 V/cm at 32 + 43 mW gave
99.60 % / 98.91 % (history, not re-optimised).

## Main findings

1. **Field 155 V/cm (Steps 16, 17).** The F'=4 − F'=3 line separation drops from 211.83 to 204.52 MHz. With the
   200 V/cm set point kept, the uniform start falls to 89.6 % (proxy), so the set point has to be re-optimised.
2. **Power (Steps 19, 23).** The extra 15 mW, put into the F'=4 beam, gains +0.2 % (SPB) and +3.1 % (uniform).
   Split at (12.5, 212), proxy: 32 + 43 mW is best for SPB (99.70 %); 26 + 49 mW is better beyond ~14 %
   uniform-like population.
3. **Where the residual is (Step 24).** All 12 F1 = 5/2 levels are cleaned to ≥ 99.5 %; the uniform-start residual
   comes entirely from the 8 F1 = 3/2 levels (mean 93.6 %), whose lines are narrow. Efficiency ≈ 99.77 − 6.2 q %
   for a fraction q in F1 = 3/2 (uniform start q = 0.4). Per-level table at (12, 212):
   `docs/state_cleanup/data/cu2_perstate_155Vcm_75mW_d12_f212.csv`; any incoherent start is Σ p_i ε_i.
4. **Laser scan with the AOM fixed (Steps 26, 27).** On axis the efficiency has narrow peaks every 17–18 MHz (EOM)
   at δ0 = 5, 22, 40 MHz (uniform 99.9 %, 2–3 MHz wide), but over the acceptance they lose 1.6–2.0 % (uniform) to
   12.5 MHz: Doppler shifts of a few MHz move off-axis molecules into the valleys. The SPB start is ≥ 99 % on axis
   from δ0 = −13 to 57 MHz.
5. **Field drift (Step 22).** ±3 V/cm with laser and AOM fixed costs ≤ 0.01 % (SPB) and ≤ 0.3 % (uniform).

Expected shape of an in-situ laser scan at f_AOM = 200 MHz, proxy (%):

| δ0/2π (MHz) | 5 | 11.5 | 12 | **12.5** | 13 | 13.5 | 22 | 23 |
|---|---|---|---|---|---|---|---|---|
| SPB | 99.44 | 99.58 | 99.67 | **99.64** | 99.55 | 99.49 | 99.48 | 99.55 |
| uniform | 95.54 | 96.88 | 97.47 | **97.59** | 97.27 | 96.88 | 95.84 | 96.01 |

## First CU2 test (CU1 and SPB off): detection

Lens start, CU2 at (12.5, 200) MHz; probe DET2, σ+ about k, 15 mW, flat top 1.5 × 1.5 cm, DET plates off. Computed
with the quantization axis along k (CU2 populations rotated onto k with d^F(π/2)); an exact check with coherent decay
and the earth field ⊥ k agrees to ≤ 3 % (Step 32, which supersedes the numbers of Steps 29–30). True J=2 left
0.49–0.50 %, 94 % of it in F1 = 3/2. Test plan `docs/state_cleanup/CU2_test_plan.md`.

| Line (role) | What it sees | Expected reading |
|---|---|---|
| R(2) F1'=7/2 F'=4 (check) | F=3 only (3 % of the residual) | S_on/S_off = 1.5–3.8 × 10⁻⁴: 13–33× below the true residual |
| Q(2) F1'=5/2 F'=3 (normalization) | all four lens levels | S_off = 1.4–2.4 photons per J=2 molecule |
| Q(2) F1'=3/2 F'=2 (residual) | all 8 F1=3/2 levels | S_on = 5.7–5.9 × 10⁻³ photons per J=2 molecule |

S_on[Q(2) F1'=3/2 F'=2] / S_off[Q(2) F1'=5/2 F'=3] = 0.24–0.41 %, i.e. the F1 = 3/2 residual divided by a model
factor 1.1–2.0 set mostly by the lens weighting, which S_off[R(2) F'=4] / S_off[Q(2) F1'=5/2 F'=3] measures (0.43 for
singlet 2/3, 0.67 equal, 0.75 triplet only). RC on P(2) F'=2 (triplet-rich start) triples the R(2) F'=4 CU2-off
signal, but not the CU2-on signal.

**Model caveat.** centrex's decay operators are pairwise (no coherence transfer in spontaneous emission), so its
results depend on the quantization axis unless the light is pure π or σ about it. The detection numbers above avoid
this; the CU2/CU1 cleanup OBEs (X/Z light, axis along E) have not been checked for it yet.

## What would most improve the prediction

1. **The J=2 population at CU2 outside the four lens levels**, and the singlet/triplet weighting (measurable with
   R(2) F'=4 and Q(2) F1'=5/2 F'=3 at CU2 off).
2. **The real loss per pass.** Each extra % costs ~0.08 % (SPB) and ~0.3 % (uniform) (200 V/cm study).
3. **The AOM diffraction efficiency across 185–215 MHz**, if the (12, 212) point is to be used.

## Model

56-level OBE (`scripts/cu2_laser_aom_scan.py`). R(2) F1'=7/2 F'=3 and F'=4 from one laser and a single-pass AOM,
B = 1 mG, one 17.5 MHz EOM on both beams, X/Z switching in antiphase. 13 passes spaced 3.5 mm, σz = 0.75 mm, 3 cm
beam (σy = 7.5 mm), Doppler sign flip between passes, 5 % loss per pass. centrex-tlf 0.3.0. Checkpoints are in
`tmp/cu2_power/vz184/` (run at v_z = 184 m/s; the models now default to the measured 156 m/s), not tracked by git; `python scripts/cu2_power_check.py summary vz184` reprints the tables.

## Earlier study: 200 V/cm, 60 mW (2026-11-01, history)

Steps 1–15 of the study log, before the field was moved to 155 V/cm and the power to 75 mW.

| Parameter | Then recommended | Report (80 + 24 mW, 10 % loss) |
|---|---|---|
| Power split | 35 + 25 mW (30 + 30 mW beyond ~13 % off m_J = 0) | 80 + 24 mW |
| δ0/2π | 12.5 MHz | 11.5 MHz (R1) |
| f_AOM | 207.5 MHz | 208.5 MHz (R1) |

| Case (200 V/cm, full grid) | SPB start | uniform start |
|---|---|---|
| Old model: 80 + 24 mW, 10 % loss, R1 | 99.57 % | 96.66 % |
| 35 + 25 mW, 5 % loss, R1 | 99.22 % | 96.66 % |
| 35 + 25 mW, 5 % loss, (12.5, 207.5) | 99.50 % | 97.26 % |

Sensitivities at (12.5, 207.5), 35 + 25 mW, proxy (general trends, not re-run at 155 V/cm):

| Change | SPB start | uniform start |
|---|---|---|
| Loss per pass 3 / 8 / 10 / 15 % (instead of 5 %) | +0.12 / −0.25 / −0.48 / −1.23 % | +0.61 / −1.17 / −2.13 / −5.03 % |
| vz = 150 / 220 m/s (instead of 184 m/s) | +0.18 / −0.26 % | +0.98 / −1.22 % |
| EOM β | optimum 1.5, flat to 1.75 | |
| X/Z switching 0.1–2 MHz, switching phase | ≤ 0.1 %, ≤ 0.03 % | ≤ 0.4 % |
| Beam height 2.4 cm instead of 3 cm | +0.08 % | +0.42 % |
