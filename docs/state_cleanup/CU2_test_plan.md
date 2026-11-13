# CU2 test plan

First test of CU2 with **CU1 and SPB off**: compare the J=2 LIF signal with CU2 on and off.

*Simulation of the expected readings:* CU2 final populations from the `cu2final` stage of
`scripts/cu2_power_check.py` (217 V/cm, (14.5, 209) MHz), probe OBE as in `scripts/cu2_detection_test.py`;
spectra `docs/state_cleanup/figures/CU2_test_detplates154_{R2F4,Q2F5o2F3,Q2F3o2F2}_spectrum_cu2_217Vcm.png`.

---

## Common conditions

| | |
|---|---|
| **CU2** | On / Off by shutter; **217 V/cm (measured), δ0 = 14.5 MHz, f_AOM = 209 MHz**, 32 + 43 mW (set point re-optimised for 217 V/cm; was 155 V/cm, (12.5, 200) MHz) |
| **Detection laser** | DET2, σ+, 15 mW, flat top 1.5 × 1.5 cm, no EOM; DET AOM −80 MHz (UV) |
| **Field in the detection region** | 154 V/cm along the laser k (DET plates on: light pure σ+ about the field). The field turns 90° adiabatically from CU2 (vertical) to DET (along k) without passing zero, so every level keeps its label and m |
| **RC** | On P(2) F1'=3/2 F'=2, so the J=2 molecules start in the triplet levels F1=5/2, F=3, mF=0, ±1 (1/3 each) |
| **Measured quantity** | S(CU2 on) / S(CU2 off) |

Simulated true residual: **0.139 % of J=2 left, 84 % of it (0.117 %) in F1=3/2**, which R(2) F'=4 and Q(2) F1'=5/2
F'=3 barely see; the F1=5/2 part is 0.022 %. (At 155 V/cm, (12.5, 200) MHz, it was 0.46 %.)

---

## Laser frequencies

IR frequencies (wavemeter), GHz; UV = 4 × IR. Model values calibrated as in CeNTREX-transition-finder
(R(2) F1'=7/2 F'=3 = 275858.770677 GHz IR at cesium 351730.547 GHz). **No cesium correction applied**: a cesium
reading differing from 351730.547 GHz shifts every line by the difference (rerun
`python scripts/cu2_test_line_frequencies.py --cs <reading in GHz>`). Full tables, including every component:
`docs/state_cleanup/CU2_test_line_frequencies.md`.

**Accuracy.** Against lines measured in the CU1 test, after removing the DET AOM: RC P(2) F'=2 −0.9 MHz IR,
Q(1) F1'=1/2 F'=1 +0.4 MHz IR, but R(1) F'=3 +11 MHz IR (unexplained). The CU1 optimum also sat 9.1 MHz IR above
the transition-finder value. Treat the model numbers as scan centres, good to a few MHz IR.

### RC (no AOM)

| Line | Laser, IR GHz |
|---|---|
| **P(2) F1'=3/2 F'=2** | **275842.0401** (measured 10-19; model 275842.04099) |
| P(2) F1'=3/2 F'=1 | 275841.96179 (model) |

### CU2 (217 V/cm, δ0 = 14.5 MHz, f_AOM = 209 MHz)

The F'=4 beam is the CU2 laser shifted by the single-pass UV AOM (+f_AOM).

| | IR GHz | UV GHz | Note |
|---|---|---|---|
| R(2) F1'=7/2 F'=3, E = 0 | 275858.77068 | 1103435.0827 | calibration line |
| F'=3 main pair at 217 V/cm | 275858.75955 | 1103435.0382 | Stark shift −44.5 MHz UV |
| F'=4 main pair at 217 V/cm | 275858.81319 | 1103435.2528 | Δ_main = 214.55 MHz UV |
| **CU2 laser (F'=3 beam)** | **275858.76318** | 1103435.0527 | main pair + 14.5 MHz UV; −7.50 MHz IR from the E = 0 line |
| F'=4 beam | — | 1103435.2617 | laser + 209 MHz UV; δ1 = +8.95 MHz |

The CU2 optimum is narrow. On axis, with the AOM fixed and the laser scanned, the uniform start stays ≥ 99 % only
for δ0 = 12–14.5 MHz (0.6 MHz IR wide; best at δ0 = 13 MHz = 275858.76281 GHz IR); ±0.5 MHz IR costs 0.5–2.5 %.
Find it by scanning the CU2 laser around 275858.763 GHz in 0.1–0.25 MHz IR steps with f_AOM = 209 MHz. Weaker
optima follow about +2.5 and +4 MHz IR higher (δ0 ≈ 25 and 31 MHz UV), see `figures/CU2_corr_scan_e217_f209.png`.

### DET (AOM −80 MHz UV: laser = line + 80 MHz UV = line + 20 MHz IR)

Plates off (E = 0):

| Line | DET laser, IR GHz |
|---|---|
| **R(2) F1'=7/2 F'=4** | **275858.83564** |
| **Q(2) F1'=5/2 F'=3** | **275848.78954** |
| **Q(2) F1'=3/2 F'=2** | **275849.04442** |
| Q(2) F1'=3/2 F'=1 | 275848.89896 |
| Q(2) F1'=5/2 F'=2 | 275848.73217 |
| R(2) F1'=7/2 F'=3 | 275858.79068 |
| R(2) F1'=5/2 F'=3 | 275859.02111 |
| R(2) F1'=5/2 F'=2 | 275858.93506 |
| R(1) F1'=5/2 F'=3 | 275855.45333 |

Plates on (154 V/cm along k): peaks of the simulated CU2-test spectra (σ+ 15 mW, triplet-only start). Bold:
strongest CU2-off peak.

| Line | CU2 off, IR GHz | CU2 on, IR GHz |
|---|---|---|
| R(2) F1'=7/2 F'=4 | **275858.83482**, 275858.83857, 275858.84019 | 275858.83494, 275858.83332, 275858.83144 |
| Q(2) F1'=5/2 F'=3 | **275848.79007**, 275848.79407, 275848.78082, 275848.78395 | 275848.79082, 275848.78907, 275848.79245 |
| Q(2) F1'=3/2 F'=2 | **275849.04475**, 275849.02350, 275849.06862 | 275849.03412, 275849.03525, 275849.05687, 275849.05825 |

---

## 1. R transition test

**Line:** R(2) F1'=7/2 F'=4

### 1.1 Transition find scan

Scan the detection laser frequency to find the transition. Scan centre (DET laser, IR): 275858.83564 GHz plates off,
275858.83482 GHz plates on (see *Laser frequencies*).

> **Note:** needs an absorption laser; CU2 is used for absorption first.

### 1.2 CU2 on/off ratio measurement

Detection laser locked on the line; record the LIF signal with CU2 on and off.

> **Note:** the CU2 laser is used for CU2 here, not for absorption, since we compare shutter on and off.

**Expected (simulation, 217 V/cm, plates on, triplet start).** With the plates on the line splits into its σ+
m components (spread over about 50 MHz UV), so lock on a peak of the CU2-off spectrum:

| Lock point (offset from the E = 0 line, UV) | DET laser, IR GHz | S(off), photons/molecule | S(on), photons/molecule | S(on)/S(off) |
|---|---|---|---|---|
| −3.3 MHz (strongest; F=3 m=0, −1) | 275858.83482 | 0.67 | 2.2e-4 | 3.3e-4 |
| +11.7 MHz (F=3 m=0) | 275858.83857 | 0.45 | 6.2e-5 | 1.4e-4 |
| +18.2 MHz (F=3 m=+1) | 275858.84019 | 0.44 | 1.3e-4 | 3.0e-4 |

The line sees only F1=5/2, so this checks that CU2 empties the triplet; S(on)/S(off) is 3–10× below the true
residual (0.14 %). Across the whole line the ratio is 1.4–4.8 × 10⁻⁴ wherever S(off) > 0.4.

---

## 2. Q transition test

**Lines:** Q(2) F1'=5/2 F'=3 (normalization), Q(2) F1'=3/2 F'=2 (residual)

### 2.1 Transition find scan

Scan the detection laser frequency to find the transitions. Scan centres (DET laser, IR): Q(2) F1'=5/2 F'=3
275848.78954 GHz (plates off) / 275848.79007 GHz (plates on); Q(2) F1'=3/2 F'=2 275849.04442 GHz / 275849.04475 GHz.

> **Note:** needs an absorption laser; CU2 is used for absorption first.

### 2.2 CU2 on/off ratio measurement

Detection laser locked on each line in turn; record the LIF signal with CU2 on and off.

> **Note:** the CU2 laser is used for CU2 here, not for absorption, since we compare shutter on and off.

**Expected (simulation, 217 V/cm, plates on, triplet start).** Lock points are peaks of the simulated spectra:

| Line, lock point (offset from the E = 0 line, UV) | DET laser, IR GHz | S(off), photons/molecule | S(on), photons/molecule | S(on)/S(off) |
|---|---|---|---|---|
| Q(2) F1'=5/2 F'=3, +2.1 MHz (normalization; strongest) | 275848.79007 | 0.90 | 3.2e-4 | 3.6e-4 |
| Q(2) F1'=3/2 F'=2, −41.2 MHz (F1=3/2 m=0; best residual peak) | 275849.03412 | 0.016 | 5.1e-4 | 3.3e-2 |
| Q(2) F1'=3/2 F'=2, −37.7 MHz (F1=3/2 m=0) | 275849.03500 | 0.036 | 4.8e-4 | 1.3e-2 |
| Q(2) F1'=3/2 F'=2, +49.8 MHz (F1=3/2 m=0) | 275849.05687 | 0.005 | 4.2e-4 | 8.8e-2 |
| Q(2) F1'=3/2 F'=2, +1.3 MHz (F1=3/2 m=−1; strongest CU2-off peak, F=3 m=−1) | 275849.04475 | 0.28 | 3.9e-4 | 1.4e-3 |

S(on)[Q(2) F1'=3/2 F'=2, −41.2 MHz] / S(off)[Q(2) F1'=5/2 F'=3] = 0.057 %, about half the F1=3/2 residual
(0.117 %): with the plates on each Q(2) F1'=3/2 F'=2 peak sees only some m levels of F1=3/2. On the −41 / −38 / +50 /
+55 MHz peaks the CU2-off background (Stark tails of the triplet lines) is small, so S(on) there is mostly F1=3/2
residual. The Q(2) F1'=3/2 F'=2 m=±1 peaks at −84 and +97 MHz (275849.02350 and 275849.06862 GHz) are dominated by
F=3 m=+1 with CU2 off.

*Simulation: centrex probe OBE with the quantization axis along k (light pure σ+ about the DET field), Doppler
average σvx 0.65 m/s, |vx| ≤ 1 m/s; CU2 populations interpolated linearly in vx (0.50 % on the 0.4 m/s grid).
Check of the start: S(off)[R(2) F'=4, −3.3 MHz] / S(off)[Q(2) F1'=5/2 F'=3, +2.1 MHz] = 0.74 for the triplet start.*

---

## 3. Other tests

-
