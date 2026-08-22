"""Create the centrex-tlf 0.2.5/Rust migration of the final J2 cleanup notebook."""
from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "SPB_cleanup_J2_J3F4_opposite_parity.ipynb"

nb = nbf.v4.new_notebook()
nb["metadata"]["kernelspec"] = {
    "display_name": "Python (CeNTREX_state_cleanup_new)",
    "language": "python",
    "name": "centrex_state_cleanup_new",
}
nb["metadata"]["language_info"] = {"name": "python", "version": "3.11"}

cells = []
cells.append(nbf.v4.new_markdown_cell("""# J=2 state cleanup: R(2) F′=3 + F′=4

Migrated from `SPB_cleanup_J2_J3F4 copy.ipynb` to `centrex-tlf 0.2.5`.

- Uses the native Rust Lindblad backend; the Julia extension is no longer required.
- Retains the opposite-parity B-state partners explicitly.
- Preserves the final two-transition scheme: 80 mW on F′=3 and 24 mW on F′=4.
- Both lasers have independent 17.5 MHz EOM phase modulation with β=1.5.
- Frequencies and Rabi rates are angular frequencies in rad/s; E is V/cm and B is Gauss.
"""))
cells.append(nbf.v4.new_code_cell("""from importlib.metadata import version
import site
import sys

# Prefer the verified wheel over a source checkout with the same package name.
sys.path.insert(0, site.getsitepackages()[-1])

import matplotlib.pyplot as plt
import numpy as np

from centrex_tlf import couplings, hamiltonian, lindblad, states, transitions, utils

print("centrex-tlf", version("centrex-tlf"))
assert version("centrex-tlf") == "0.2.5"
"""))
cells.append(nbf.v4.new_code_cell("""plt.rcParams.update({
    "axes.titlesize": 24, "axes.labelsize": 20,
    "xtick.labelsize": 16, "ytick.labelsize": 16,
    "lines.linewidth": 2,
})

def gen_multipass_para(npass, loss, first_center, spacing):
    order = utils.multipass.multipass_prism_order(npass)
    xlocs = [first_center + i * spacing for i in range(npass)]
    ylocs = [0.0] * npass
    relative_power = [(1.0 - loss) ** (pass_number - 1) for pass_number in order]
    return xlocs, ylocs, relative_power
"""))
cells.append(nbf.v4.new_code_cell("""trans = [transitions.R2_F1_7o2_F3, transitions.R2_F1_7o2_F4]
polarizations = [
    [couplings.polarization_X, couplings.polarization_Z],
    [couplings.polarization_X, couplings.polarization_Z],
]

ground_mains = [
    1 * states.CoupledBasisState(
        J=2, F=3, F1=5/2, mF=2, I1=1/2, I2=1/2, Ω=0, P=1,
        electronic_state=states.ElectronicState.X,
    ),
    1 * states.CoupledBasisState(
        J=2, F=3, F1=5/2, mF=0, I1=1/2, I2=1/2, Ω=0, P=1,
        electronic_state=states.ElectronicState.X,
    ),
]
excited_mains = [
    1 * states.CoupledBasisState(
        J=3, F=3, F1=7/2, mF=3, I1=1/2, I2=1/2, Ω=1, P=-1,
        electronic_state=states.ElectronicState.B,
    ),
    1 * states.CoupledBasisState(
        J=3, F=4, F1=7/2, mF=1, I1=1/2, I2=1/2, Ω=1, P=-1,
        electronic_state=states.ElectronicState.B,
    ),
]
transition_selectors = couplings.generate_transition_selectors(
    trans, polarizations, ground_mains=ground_mains, excited_mains=excited_mains
)
"""))
cells.append(nbf.v4.new_code_cell("""# Experimental settings retained from the final notebook.
E_FIELD = np.array([0.0, 0.0, 250.0])       # V/cm
B_FIELD = np.array([0.0, 0.0, 1e-3])        # Gauss; resolves ±mF numerically
vz, vx = 184.0, 0.0                         # m/s
laser_frequency = 1.103e15                  # Hz, for Doppler conversion

power0_default = 80e-3                      # F′=3 laser, W
power1_default = 0.30 * power0_default       # F′=4 laser, W
detuning0_default = 18.5 * 2*np.pi*1e6      # rad/s
detuning1_default = -6.5 * 2*np.pi*1e6      # rad/s
eom0_default = 17.5 * 2*np.pi*1e6           # rad/s
eom1_default = 17.5 * 2*np.pi*1e6           # rad/s
beta0_default = beta1_default = 1.5
polarization_frequency_default = 0.5 * 2*np.pi*1e6

sigma_z_default, sigma_y_default = 0.75e-3, 6e-3
xlocs_default, ylocs_default, pass_power = gen_multipass_para(13, 0.10, 3.5e-3, 3.5e-3)
intensities_per_watt = tuple(
    np.asarray(pass_power) / (2*np.pi*sigma_z_default*sigma_y_default)
)
T_END = 0.05 / vz
"""))
cells.append(nbf.v4.new_code_cell("""%%time
obe_system = lindblad.generate_OBE_system_transitions(
    trans,
    transition_selectors,
    verbose=True,
    qn_compact=True,
    E=E_FIELD,
    B=B_FIELD,
    retain_opposite_parity_levels=True,
    Jmax_X=5,
    Jmax_B=5,
)
print("levels:", len(obe_system.QN))
print("main couplings:", [abs(c.main_coupling) for c in obe_system.couplings])
"""))
cells.append(nbf.v4.new_code_cell("""%%time
params = lindblad.LindbladParameters()
t = params.time()

power0 = params.real("power0", power0_default)
power1 = params.real("power1", power1_default)
detuning0 = params.real("detuning0", detuning0_default)
detuning1 = params.real("detuning1", detuning1_default)
eom0 = params.real("eom0", eom0_default)
eom1 = params.real("eom1", eom1_default)
beta0 = params.real("beta0", beta0_default)
beta1 = params.real("beta1", beta1_default)
pol_frequency = params.real("pol_frequency", polarization_frequency_default)
y0 = params.real("y0", 0.0)
sigma_z = params.real("sigma_z", sigma_z_default)
sigma_y = params.real("sigma_y", sigma_y_default)
intensity_profile = params.real("intensities_per_watt", intensities_per_watt)
xlocs = params.real("xlocs", tuple(xlocs_default))
ylocs = params.real("ylocs", tuple(ylocs_default))

drives = []
for k, (power, beta, eom) in enumerate(((power0, beta0, eom0), (power1, beta1, eom1))):
    main_coupling = params.real(
        f"main_coupling{k}", float(abs(obe_system.couplings[k].main_coupling))
    )
    spatial = lindblad.multipass_2d_rabi(
        vz*t, y0, intensity_profile, xlocs, ylocs,
        sigma_z, sigma_y, main_coupling, 2.6675506e-30,
    )
    drives.append(power**0.5 * lindblad.phase_modulation(t, beta, eom) * spatial)

doppler_per_mps = float(utils.detuning.velocity_to_detuning(1.0, laser_frequency))
for selector, drive, detuning in zip(
    transition_selectors, drives, (detuning0, detuning1)
):
    params.bind(selector.Ω, drive, finalize=False)
    params.bind(selector.δ, detuning + vx*doppler_per_mps, finalize=False)

# Complementary X/Z switching: when laser 0 is X, laser 1 is Z, and vice versa.
z_gate = lindblad.square_wave(t, pol_frequency, 0.0)
px0, pz0 = transition_selectors[0].polarization_symbols
px1, pz1 = transition_selectors[1].polarization_symbols
params.bind(px0, 1.0-z_gate, finalize=False)
params.bind(pz0, z_gate, finalize=False)
params.bind(px1, z_gate, finalize=False)
params.bind(pz1, 1.0-z_gate, finalize=False)
params._finalize()

prepared = lindblad.prepare_lindblad_problem(
    obe_system, params, backend="rust", hamiltonian_representation="decomposed"
)
"""))
cells.append(nbf.v4.new_code_cell("""ground_idx = np.asarray(states.QuantumSelector(
    J=2, electronic=states.ElectronicState.X
).get_indices(obe_system.QN), dtype=int).ravel()
excited_idx = np.asarray(states.QuantumSelector(
    electronic=states.ElectronicState.B
).get_indices(obe_system.QN), dtype=int).ravel()

rho0 = np.zeros(obe_system.H_symbolic.shape, dtype=np.complex128)
rho0[ground_idx, ground_idx] = 1.0 / ground_idx.size
print("initial J=2 states:", ground_idx.size, "trace:", np.trace(rho0).real)
"""))
cells.append(nbf.v4.new_code_cell("""saveat = np.append(np.arange(0.0, T_END, 1e-6), T_END)
solver_options = dict(
    solver="dopri5", execution_mode="expanded_sparse",
    dt=2e-9, reltol=1e-7, abstol=1e-9, maxiters=500_000,
)

results = lindblad.solve_lindblad(
    prepared, rho0, (0.0, T_END), saveat=saveat,
    output="populations", output_when="saveat", **solver_options,
)
populations = results.values
"""))
cells.append(nbf.v4.new_code_cell("""remaining_j2 = populations[-1, ground_idx].sum().real
photons = np.trapezoid(populations[:, excited_idx].sum(axis=1), x=results.t) * hamiltonian.Γ
print("final trace =", populations[-1].sum().real)
print("remaining J=2 population =", remaining_j2)
print("cleanup efficiency =", 1.0-remaining_j2)
print("photons per molecule =", photons.real)

fig, ax = plt.subplots(figsize=(10, 7))
for idx in ground_idx:
    ax.plot(results.t*vz, populations[:, idx], label=str(obe_system.QN[idx].largest))
ax.set(xlabel="Distance / m", ylabel="Population", title="X, J=2 ground-state population")
ax.grid(alpha=0.3)
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
plt.show()
"""))
cells.append(nbf.v4.new_markdown_cell("""## Parameter scans

The prepared Rust problem is reused. Scan values below are angular frequencies. The returned value is the final population remaining in the 20 initial J=2 states; cleanup efficiency is `1 - remaining`.
"""))
cells.append(nbf.v4.new_code_cell("""# Coarse scan of both carrier detunings. Reduce/increase point counts as needed.
detuning0_scan = 2*np.pi*1e6*np.linspace(-40, 40, 41)
detuning1_scan = 2*np.pi*1e6*np.linspace(-40, 40, 41)

scan = lindblad.grid_scan(
    prepared, rho0, (0.0, T_END),
    scan={"detuning0": detuning0_scan, "detuning1": detuning1_scan},
    output="populations", output_indices=ground_idx,
    output_when="final", parallel=True, **solver_options,
)
remaining = scan.values.real.sum(axis=-1).reshape(len(detuning0_scan), len(detuning1_scan))
cleanup = 1.0-remaining
"""))
cells.append(nbf.v4.new_code_cell("""fig, ax = plt.subplots(figsize=(8, 6))
mesh = ax.pcolormesh(
    detuning1_scan/(2*np.pi*1e6), detuning0_scan/(2*np.pi*1e6),
    100*cleanup, shading="auto",
)
fig.colorbar(mesh, ax=ax, label="Cleanup efficiency / %")
ax.set(xlabel="F′=4 detuning / MHz", ylabel="F′=3 detuning / MHz")
plt.show()
"""))

nb["cells"] = cells
nbf.write(nb, OUT)
print(OUT)
