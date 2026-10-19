# CU1 efficiency test — handoff

Status as of 2026-10-15. Plan and supporting simulations for measuring the cleanup efficiency of CU1
(removal of X, J=1 before SPB). Nothing in this document has been measured yet; every number is a model
prediction unless stated otherwise.

## 1. What CU1 is and what the model predicts

| item | value |
|---|---|
| transition | R(1) F1'=5/2, F'=2, in 150 V/cm |
| set point | +19.0 MHz from the main pair \|J=1, F1=3/2, F=2, mF=+1> -> \|J'=2, F1'=5/2, F'=2, mF'=+2, P'=+> |
| laser | 80 mW, 3 cm beam height (sigma_y = 7.5 mm), 13 passes, 10 %/pass loss, X/Z polarization switching at 0.5 MHz, no EOM |
| molecules | 184 m/s; transverse distribution from the 25 kV lens simulation (sigma_vx = 0.65 m/s, sigma_y = 10.95 mm at CU1) |
| detection acceptance | \|vx\| <= 1 m/s, \|y0\| <= 8.5 mm |

Predicted efficiencies (report `docs/state_cleanup/main.tex`, Secs. 5.9-5.10):

| initial state | beam | on axis | mean in acceptance |
|---|---|---|---|
| 8 SPB-related J=1 levels (mainly \|mJ\|=1) | 3 cm | 99.01 % | 98.70 % |
| same | 4.8 cm | 96.4 % | 97.08 % |
| AP1 triplet levels (see Sec. 2), +19.0 MHz | 4.8 cm | 99.75-99.77 % | not computed |
| AP1 singlet level, +19.0 MHz | 4.8 cm | 99.78 % | not computed |

The J=1 sublevels form two groups 2.4 MHz apart in 150 V/cm: four mainly-mJ=0 levels (upper) and the eight
SPB-related, mainly-|mJ|=1 levels (lower). At +19.0 MHz and vx = 0 the laser is close to the upper group and
1.7-1.9 MHz below the lower one; molecules with |vx| ~ 0.5 m/s see the lower group on resonance in half of the
passes. The hardest levels are |F1=3/2, F=2, mF=+-2> and |F1=1/2, F=0> (weak couplings).

## 2. Test configuration

### Upstream modules

| module | setting | why |
|---|---|---|
| source | on | |
| rotational cooling (RC) | on, optical pumping on **P(2) F1'=3/2 F'=2** | enhances J=0, F=1 (nuclear-spin triplet) as in the AP paper; the normal CeNTREX setting P(2) F'=1 enhances J=0, F=0 (singlet) |
| SPA AP1 (J=0 -> J=1, mJ=0) | **on / off alternated** | produces a large, known J=1 population; on-off difference removes the J=1 left by RC |
| SPA AP2 (J=1 -> J=2) | off | |
| electrostatic lens | on or off, whichever gives more J=1 signal | J=1, mJ=0 is weak-field seeking only below ~15 kV/cm (J=2, mJ=0 below ~44 kV/cm) and at low field has ~4x the Stark shift of J=2, mJ=0: the lens over-focuses or defocuses J=1. Either way the J=1 transverse distribution is not the J=2 beam simulation. |
| CU1 field plates | 150 V/cm | as in the model |
| CU1 laser | **on / off alternated** (shutter) | |
| SPB microwaves | **off** (fields as normal) | CU1 sends ~46 % of the pumped molecules to J=2; SPB would move some of them back to J=1 |
| CU2 laser | off | its J=1 branching is ~1e-5; off for a clean test |

"All modules off" (thermal J=1) was considered first and rejected: too little signal.

### What AP1 puts into J=1

AP1 uses z-polarized microwaves (Delta mJ = 0) and does not act on the nuclear spins, so all four J=0 hyperfine
sublevels are transferred (AP paper, arXiv:2511.16626, Fig. 5):

- triplet S=1 (J=0, F=1, mF=-1,0,+1) -> J=1, mJ=0, S=1 -> at zero field \|J=1, F1=3/2, F=2, mF=-1,0,+1>
- singlet S=0 (J=0, F=0) -> J=1, mJ=0, S=0 -> at zero field \|J=1, F1=3/2, F=1, mF=0>

In CU1's 150 V/cm these are the four mainly-mJ=0 levels of the upper group (|mJ|=1 content <= 0.002), i.e.
**not** the SPB-related levels. The test therefore measures CU1 on the easier group; the residual it leaves is,
however, in the hard levels (molecules pumped and returned there by decay):

| CU1 residual, +19.0 MHz, 4.8 cm, on axis, triplet start | value |
|---|---|
| J=1 left | 0.235 % |
| in \|F=0> | 32 % of it |
| in \|F1=3/2, F=2, mF=+-2> | 64 % of it |
| all other sublevels | ~3-4 % of it |

## 3. Detection

Present probe: 2.5 cm (along the beam) x 1.7 cm (height), single frequency, fixed linear polarization, no phase
EOM, adjustable intensity; zero field. Polarization switching cannot easily be added.

### The two lines

| | R(1) F1'=5/2, F'=3 | Q(1) F1'=1/2, F'=1 |
|---|---|---|
| sees | only F=2 (5 sublevels) | all 12 J=1 sublevels (one frequency, J=1 spans 0.21 MHz at E=0) |
| dark states | none (5 ground vs 7 excited, F'=3 decays inside J=1 only to F=2) | yes (12 ground vs 3 excited, fixed polarization) |
| photons per molecule, saturated | 2.1 for every F=2 sublevel, any polarization (matches the AP paper's 2.10) | 0-4.4 depending on sublevel and on polarization relative to the molecules' axis |
| key property | exact F=2 count | F=0 is a single isotropic sublevel, never dark: ~4.4 photons |
| blind to | F=0 (32 % of the residual), F=1 levels, singlet | pi light: F=2 mF=+-2 and mF=0 of both F=1 levels |

Q(1) photons per molecule (saturated, 50 mW/cm2 peak, probe centre, vx = 0); "par" = polarization along the
molecules' quantization axis, "perp" = perpendicular (perp is approximate, see Sec. 6):

| J=1 sublevel | par | perp |
|---|---|---|
| F1=3/2, F=2, mF=+-2 | 0 | 3.1 |
| F1=3/2, F=2, mF=+-1 | 2.9 | 1.3 |
| F1=3/2, F=2, mF=0 | 4.4 | 0.9 |
| F1=3/2, F=1, mF=+-1 | 2.5 | 1.1 |
| F1=3/2, F=1, mF=0 (singlet) | 0 | 2.9 |
| F1=1/2, F=1, mF=+-1 | 3.1 | 1.4 |
| F1=1/2, F=1, mF=0 | 0 | 3.7 |
| F1=1/2, F=0 | 4.4 | 3.7 |

### Saturation requirement (no EOM)

R(1) F'=3 photons per F=2 molecule (range over mF), probe taken as a Gaussian with 1/e2 diameters 25 x 17 mm;
power = peak intensity x 1.67 cm2:

| peak intensity (power) | centre, vx=0 | vx=1 m/s, y=0 | vx=0.5 m/s, y=8.5 mm | vx=1 m/s, y=8.5 mm |
|---|---|---|---|---|
| 2 mW/cm2 (3 mW) | 2.10 | 0.58-0.88 | 0.29-0.48 | 0.09-0.15 |
| 10 mW/cm2 (17 mW) | 2.10 | 1.77-1.95 | 1.18-1.51 | 0.41-0.65 |
| 50 mW/cm2 (83 mW) | 2.10 | 2.10 | 2.08-2.09 | 1.47-1.75 |

The CU1 residual is concentrated at large |vx| and large |y| (CU1 is 93-95 % at |vx| = 1 m/s), where the probe is
weakest. An unsaturated probe under-counts the residual and over-estimates CU1. Saturation must be shown with a
probe-intensity scan (the on/off ratio must plateau). If the power is not enough, add a slit before the probe
(as in the AP paper) to restrict y and vx.

The probe region must be at low field: Q(1) closure (fraction of decays back to J=1) is 99.7 % at 2 V/cm,
98.3 % at 5 V/cm, 85 % at 20 V/cm and 67 % at 150 V/cm.

## 4. Measurement sequence

0. **Checks.**
   - CU1 at 150 V/cm, 80 mW, 3 cm; record the power at the entrance and after the last pass.
   - Find both probe lines.
   - Measure the probe saturation curve and set the probe intensity on the plateau.
   - Background: probe detuned off the line (or beam blocked), measured with CU1 on and with CU1 off
     (this checks for scattered CU1 light).
   - PMT linearity: S_on is ~0.1-1 % of S_off, so ~100x dynamic range is needed.
1. **Low-power CU1 spectrum.**
   - Reduce CU1 to ~1 mW and scan its frequency, recording the J=1 depletion.
   - The Stark components fix the detuning zero and check the actual field: \|mF'\|=2 near 0 MHz, the
     \|mF'\|=1 groups at +18.5 to +20.9 MHz.
2. **Full-power detuning scan.**
   - 80 mW, +16 to +24 MHz in 0.5 MHz steps, CU1 on/off at each point.
   - Compare with the shape of report Table `tab:j1_3cm`.
3. **Long run at the set point** (+19.0 MHz or the measured optimum): interleave until sigma_eta <~ 0.2 %.
4. **Q(1) at two orthogonal linear polarizations** (half-wave plate), one along the last field direction.
5. **Optional.**
   - CU1 power scan, 10-80 mW.
   - CU1 polarization switching off: the model predicts a clear drop from dark states.
   - TOF slices: the model uses only 184 m/s.

Every quantity uses AP1-on minus AP1-off signals, each with CU1 on and off (4 signal sets per line).

## 5. Analysis

Signals: S = TOF-integrated probe counts minus background, normalized to a source monitor. Do **not** normalize
to a downstream J=2 signal, because CU1 adds J=2. Interleave CU1 on/off shot by shot or in short blocks.

```
Delta S(CU1 x) = S(AP1 on, CU1 x) - S(AP1 off, CU1 x)
R_line         = Delta S(CU1 on) / Delta S(CU1 off)
```

With S = eps * sum_k y_k N_k and the same photon detection efficiency eps for both lines:

- R(1) F'=3 gives the F=2 population exactly: N_F2 = S_R / (2.1 eps). With CU1 off and a triplet start,
  all of J=1 is in F=2, so N_off = S_R^off / (2.1 eps).
- Q(1) gives an upper bound on F=0, because every term of S_Q is >= 0:
  N_F0^on <= S_Q^on / (4.4 eps).
- Combined (eps cancels):

```
(N_F2^on + N_F0^on) / N_off  <=  (S_R^on / 2.1 + S_Q^on / 4.4) / (S_R^off / 2.1)
```

This bounds the F=2 + F=0 part of the residual; the F=1 levels (~3 % of the residual in the model) are not
covered, so for all of J=1 it is "<~". If the singlet is also populated with CU1 off, the R line does not see it,
N_off is under-estimated and the bound stays valid but becomes looser.

Predicted readings (CU1 +19.0 MHz, 4.8 cm, on axis, triplet start; probe saturated, centre, vx = 0):

| quantity | molecules keep the CU1 axis, probe par | keep the axis, probe perp | alignment scrambled |
|---|---|---|---|
| true J=1 residual | 0.235 % | 0.235 % | 0.235 % |
| R(1) F'=3 ratio | 0.15 % | 0.15 % | 0.15 % |
| Q(1) ratio | 0.10 % | 0.64 % | 0.32 % |
| combined bound | 0.232 % | 0.326 % | 0.302 % |

How to read them:
- R(1) alone under-reads the residual by about 1/3 (the F=0 part).
- Q(1) alone is off by a factor 0.4 to 2.7, depending on the molecules' alignment.
- The Q(1) par/perp pair reads about 6x apart if the molecules keep their alignment and the same if it is
  scrambled; this tests both the alignment and the CU1 model.
- Singlet cross-term: if the singlet carries 2/3 of the AP1 population (normal RC), the R(1) ratio rises from
  0.15 % to ~0.44 %. CU1 pumps singlet molecules into F=2, where the R line sees them, while with CU1 off they
  are invisible to it.

Precision: R ~ 0.002, so the error on eta is dominated by S_on and the background. For +-0.2 % the background
must be known to ~0.1 % of S_off.

## 6. Assumptions and known weaknesses

Probe model (`scripts/detection_level_diagrams.py`, functions `build` and `compute`):

- **States.** X, J=1 plus the single probed B hyperfine level (P'=+). Decays out of J=1 count as loss.
  v=0 only: no vibrational leakage, worth a few % at <= 4 photons.
- **Fields and probe.**
  - E = 0, B = 1e-3 G along z. Earth's field (~0.5 G, direction unknown) is not included; it would
    partly remix the Q(1) dark states, by about 0.5-1 rad of Larmor rotation over the probe time.
  - Probe: Gaussian, single pass, single frequency on resonance, zero linewidth, pure linear polarization.
- **Molecules.** Straight line at 184 m/s, ±3 sigma of the probe (204 µs). Transverse velocity enters only as a
  Doppler detuning. Initial state: one zero-field \|F, mF> eigenstate, i.e. the molecules are assumed to keep
  the CU1 field axis up to the probe.
- **Photon counting.** Photons = integral of Gamma * P_excited, with the same collection efficiency for every
  photon. The real PMT solid angle and the pi/sigma emission patterns make the collection efficiency depend on
  sublevel and line, possibly by tens of %. This matters for the relative calibration of 2.1 versus 4.4;
  check it with the measured Q/R ratio with CU1 off (model: ~3.3 for a triplet start).
- **Decay-induced ground-state coherence is missing.** The solver has one collapse operator per (e, g) pair,
  so it leaves out the coherences that decay creates between same-mF levels of different F (14-210 kHz apart).
  - This affects **all Q(1) numbers, both polarizations**, including the F=0 value used in the bound:
    the 4.4 (par) vs 3.7 (perp) difference for an isotropic level is a symptom, so the F=0 yield is uncertain
    at the ~20 % level.
  - R(1) F'=3 is **not** affected: there are no dark states, so its yield is fixed by branching alone.
- **Validation.** The only experimental check of the probe model is R(1) F'=3 = 2.10 photons (AP paper, which
  used a smaller probe with phase modulation).

CU1 model:

- **Coverage of the predictions.**
  - The per-level CU1 results for the AP1 levels exist only at 4.8 cm and on axis (vx = 0, y0 = 0).
  - No 3 cm (vx, y0) maps exist yet for the AP1 levels.
- **Transverse distribution.** The beam-simulation transverse distribution is for lensed J=2. The J=1
  distribution in this test (lens on or off) is not simulated.
- **Field path.** Between SPA and CU1 the molecules pass through zero field. A transfer of F=2,
  mF=-1,0,+1 population into mF=+-2 there would add hard-to-clean molecules (CU1 leaves 9.4 % of mF=+-2 at
  +19.0 MHz). A 1 % admixture would add ~0.09 % residual, comparable with the signal. If the measured residual
  is well above the prediction, suspect this first.

## 7. Open items (in order of importance)

1. CU1 (vx, y0) maps at 3 cm for the four AP1 levels (and a uniform start for the RC remainder), combined with
   the probe yield over (vx, y): the predicted S_on/S_off of each line versus probe intensity. A few minutes
   per map with `scripts/cu1_3cm_scan.py` machinery.
2. Q(1) yields with a full per-q dissipator (decay-induced coherence), Earth's field in an arbitrary direction,
   and a PMT collection-efficiency estimate. The system has < 20 levels, so a dense Lindblad solver is cheap.
   This decides how reliable the 4.4 in the bound is.
3. J=1 transverse distribution at CU1 and at the probe for the chosen lens setting (beam simulation with the
   J=1, mJ=0 Stark curve).
4. Simulated low-power CU1 depletion spectrum for the frequency/field calibration of step 1.
5. Archive `tmp/detection/probe_yields.json` into `results/SPB_J1_cleanup_scans.h5` and list it in
   `results/scan_results.md`; add this test to the LaTeX report.
6. The Q(1) closure-versus-field and Stark-turnover numbers above came from ad-hoc calculations that are not
   in the repo; re-derive them in a script if they are needed in the report.

## 8. Files

| what | where |
|---|---|
| level diagrams of the two probe lines (slides) | `docs/state_cleanup/figures/slides/detection_level_diagram_R1_F3_0Vcm.{png,pdf}`, `..._Q1_F1_0Vcm.{png,pdf}` |
| "what each line tells about CU1" (slides) | `docs/state_cleanup/figures/slides/detection_cu1_information.{png,pdf}` |
| probe model, diagrams, schematic | `scripts/detection_level_diagrams.py` (`compute`, `plot`, `info`) |
| probe photon yields (grid: 3 intensities x 3 vx x 2 y x 2 polarizations x 2 lines) | `tmp/detection/probe_yields.json` (not yet archived) |
| CU1 per-level trajectories used for the residual | `results/SPB_J1_cleanup_scans.h5`, run `trajectories_cu1_vx0_19p0_20p2MHz_4p8cm` |
| CU1 at 3 cm, SPB-related start | same file, runs `doppler_cu1_spb_3cm_det19MHz`, `detuning_vx_y_cu1_spb_3cm_coarse` |
| index of all stored runs | `results/scan_results.md` |
| report | `docs/state_cleanup/main.tex`, Secs. 5.8-5.10 (`sec:j1_traj19`, `sec:j1_present`) |
| SPA and detection reference | Grasdijk et al., "Adiabatic passage of 205TlF with microwaves in a cryogenic beam", arXiv:2511.16626 |
