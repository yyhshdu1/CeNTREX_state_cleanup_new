"""Replace the scan sections in the migrated J=1 and J=2 notebooks."""
from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]


def md(text: str):
    return nbf.v4.new_markdown_cell(text)


def code(text: str):
    return nbf.v4.new_code_cell(text)


def update_j1() -> None:
    path = ROOT / "SPB_cleanup_J1-J2F2.ipynb"
    nb = nbf.read(path, as_version=4)
    # Current J=1 apparatus: no phase EOM and 0.5 MHz polarization switching.
    nb.cells[6].source = nb.cells[6].source.replace(
        "phase_frequency_default = 0.0 * 2 * np.pi * 1e6\n"
        "phase_depth_default = 1.5\n"
        "polarization_frequency_default = 1.0 * 2 * np.pi * 1e6",
        "polarization_frequency_default = 0.5 * 2 * np.pi * 1e6",
    )
    nb.cells[8].source = nb.cells[8].source.replace(
        'phase_frequency = params.real("phase_frequency", phase_frequency_default)\n'
        'phase_depth = params.real("phase_depth", phase_depth_default)\n',
        "",
    ).replace(
        "drive = laser_power**0.5 * lindblad.phase_modulation(\n"
        "    t, phase_depth, phase_frequency\n"
        ") * spatial_rabi",
        "# No phase EOM is used in the current J=1 cleanup configuration.\n"
        "drive = laser_power**0.5 * spatial_rabi",
    )
    imports = [
        "from scripts.scan_hdf import load_scan_hdf, save_scan_hdf",
        "from scripts.scan_progress import grid_scan_with_progress",
    ]
    for scan_import in imports:
        if scan_import not in nb.cells[1].source:
            nb.cells[1].source += f"\n{scan_import}\n"
    nb.cells = nb.cells[:21] + [
        md("""## Structured parameter scans

All scan results below are the **final population remaining in the initial J=1 manifold**. Therefore `cleanup_efficiency = 1 - remaining_population`.

The array dimensions always follow the insertion order of `scan_axes`. HDF5 stores that order, every axis and unit, fixed experimental conditions, model construction settings, solver tolerances, and the package version. Frequencies are converted to MHz before storage for readability."""),
        code("""scan_solver_options = dict(
    solver="dopri5",
    execution_mode="expanded_sparse",
    output="populations",
    output_when="final",
    dt=2e-9,
    reltol=1e-7,
    abstol=1e-9,
    maxiters=500_000,
    parallel=True,
)

model_metadata = {
    "notebook": "SPB_cleanup_J1-J2F2.ipynb",
    "centrex_tlf_version": version("centrex-tlf"),
    "transition": str(trans[0]),
    "retain_opposite_parity_levels": True,
    "number_of_levels": len(obe_system.QN),
    "initial_manifold": "X, J=1, uniform population",
    "initial_state_count": int(ground_idx.size),
}

fixed_conditions = {
    "electric_field_V_per_cm": [0.0, 0.0, 150.0],
    "magnetic_field_G": [0.0, 0.0, 1e-3],
    "interaction_length_m": 0.05,
    "molecular_vz_m_per_s": vz,
    "molecular_vx_m_per_s": vx,
    "laser_frequency_Hz": laser_frequency,
    "phase_EOM_used": False,
    "polarization_switch_frequency_MHz": polarization_frequency_default/(2*np.pi*1e6),
    "sigma_z_mm": sigma_z_default*1e3,
    "sigma_y_mm": sigma_y_default*1e3,
    "multipass_count": len(xloctions),
    "multipass_relative_power": list(pow),
    "main_coupling": float(abs(obe_system.couplings[0].main_coupling)),
}
"""),
        md("""### Power × carrier-detuning scan

The current J=1 cleanup configuration does not use a phase EOM. This scan therefore varies only optical power and carrier placement. The stored result has dimensions `[power, detuning]`; no transpose is applied."""),
        code("""# Use coarse axes for testing; increase the point counts for the final production scan.
power_mW = np.linspace(50.0, 110.0, 3)
detuning_MHz = np.linspace(-16.0, 30.0, 23)

scan_axes = {
    "laser_power": power_mW*1e-3,
    "laser_detuning": detuning_MHz*2*np.pi*1e6,
}

frequency_scan = grid_scan_with_progress(
    prepared, rho0, (0.0, t_end), scan=scan_axes, **scan_solver_options
    , chunk_size=32, description="J=1 power/detuning"
)
frequency_populations = frequency_scan.values.reshape(
    len(power_mW), len(detuning_MHz), len(obe_system.QN)
)
frequency_remaining = frequency_populations[..., ground_idx].sum(axis=-1).real
frequency_cleanup = 1.0-frequency_remaining
print("remaining-population shape [power, detuning] =", frequency_remaining.shape)
"""),
        code("""frequency_file = save_scan_hdf(
    "SPB_J1_cleanup_scans.h5",
    "power_detuning_no_phase_eom",
    frequency_remaining,
    axes={
        "laser_power_mW": (power_mW, "mW"),
        "laser_detuning_MHz": (detuning_MHz, "MHz"),
    },
    quantity="final_population_remaining_in_X_J1",
    fixed_parameters=fixed_conditions,
    model=model_metadata,
    solver=scan_solver_options,
    notes="Uniform initial population over all retained X,J=1 states.",
    overwrite=False,  # Protect an existing long-running scan; change run_name for a new run.
)
print("saved:", frequency_file)
"""),
        code("""fig, ax = plt.subplots(figsize=(8, 5))
mesh = ax.pcolormesh(
    detuning_MHz, power_mW, 100*frequency_cleanup, shading="auto"
)
fig.colorbar(mesh, ax=ax, label="Cleanup efficiency / %")
ax.set(xlabel="Carrier detuning / MHz", ylabel="Laser power / mW")
plt.show()
"""),
        md("""### Vertical position × beam diameter scan

The runtime parameter is Gaussian `sigma_y`, while the experimental scan axis is the quoted full diameter `4*sigma_y`. The HDF5 file stores the diameter in mm and documents this conversion."""),
        code("""y_position_mm = np.linspace(-12.0, 12.0, 50)
beam_diameter_mm = np.linspace(10.0, 50.0, 5)

geometry_scan = grid_scan_with_progress(
    prepared, rho0, (0.0, t_end),
    scan={
        "y0": y_position_mm*1e-3,
        "sigma_y": beam_diameter_mm*1e-3/4.0,
    },
    **scan_solver_options,
    chunk_size=32,
    description="J=1 beam geometry",
)
geometry_populations = geometry_scan.values.reshape(
    len(y_position_mm), len(beam_diameter_mm), len(obe_system.QN)
)
# Preserve scan-axis order: [vertical position, beam diameter].
geometry_remaining = geometry_populations[..., ground_idx].sum(axis=-1).real
geometry_cleanup = 1.0-geometry_remaining
print("remaining-population shape [y position, beam diameter] =", geometry_remaining.shape)
"""),
        code("""geometry_file = save_scan_hdf(
    "SPB_J1_cleanup_scans.h5",
    "vertical_position_beam_diameter",
    geometry_remaining,
    axes={
        "vertical_position_mm": (y_position_mm, "mm"),
        "beam_diameter_4sigma_mm": (beam_diameter_mm, "mm"),
    },
    quantity="final_population_remaining_in_X_J1",
    fixed_parameters={
        **fixed_conditions,
        "laser_power_mW": power_default*1e3,
        "carrier_detuning_MHz": laser_detuning_default/(2*np.pi*1e6),
        "phase_EOM_used": False,
    },
    model=model_metadata,
    solver=scan_solver_options,
    notes="beam_diameter_4sigma_mm is converted to runtime sigma_y by sigma_y=diameter/4.",
    overwrite=False,  # Protect an existing long-running scan; change run_name for a new run.
)
print("saved:", geometry_file)
"""),
        code("""fig, ax = plt.subplots(figsize=(8, 5))
mesh = ax.pcolormesh(
    beam_diameter_mm, y_position_mm, 100*geometry_cleanup, shading="auto"
)
fig.colorbar(mesh, ax=ax, label="Cleanup efficiency / %")
ax.set(xlabel="Beam diameter (4σ) / mm", ylabel="Vertical position / mm")
plt.show()
"""),
        md("""### HDF5 readback example

Loading does not require recreating the OBE system. It returns the result, axes, units, fixed conditions, model information and solver settings together."""),
        code("""saved_run = load_scan_hdf("SPB_J1_cleanup_scans.h5", "power_detuning_no_phase_eom")
print(saved_run["results"].shape)
print(saved_run["units"])
print(saved_run["fixed_parameters"])
"""),
    ]
    nbf.write(nb, path)


def update_j2() -> None:
    path = ROOT / "SPB_cleanup_J2_J3F4_opposite_parity.ipynb"
    nb = nbf.read(path, as_version=4)
    imports = [
        "from scripts.scan_hdf import load_scan_hdf, save_scan_hdf",
        "from scripts.scan_progress import grid_scan_with_progress",
    ]
    for scan_import in imports:
        if scan_import not in nb.cells[1].source:
            nb.cells[1].source += f"\n{scan_import}\n"
    nb.cells = nb.cells[:10] + [
        md("""## Structured parameter scans

The final J=2 scheme has two independent transitions, `R(2), F′=3` and `R(2), F′=4`. The most direct scan varies both carrier detunings while keeping powers, both EOM settings, fields, multipass geometry and polarization switching fixed.

The saved quantity is final population remaining in the initial J=2 manifold; cleanup efficiency is `1 - remaining`. Frequencies are stored in MHz, although the solver receives rad/s."""),
        code("""scan_solver_options = {
    **solver_options,
    "output": "populations",
    "output_when": "final",
    "parallel": True,
}

model_metadata = {
    "notebook": "SPB_cleanup_J2_J3F4_opposite_parity.ipynb",
    "centrex_tlf_version": version("centrex-tlf"),
    "transitions": [str(item) for item in trans],
    "retain_opposite_parity_levels": True,
    "number_of_levels": len(obe_system.QN),
    "initial_manifold": "X, J=2, uniform population",
    "initial_state_count": int(ground_idx.size),
}

fixed_conditions = {
    "electric_field_V_per_cm": E_FIELD.tolist(),
    "magnetic_field_G": B_FIELD.tolist(),
    "interaction_length_m": 0.05,
    "molecular_vz_m_per_s": vz,
    "molecular_vx_m_per_s": vx,
    "laser_frequency_Hz": laser_frequency,
    "Fprime3_power_mW": power0_default*1e3,
    "Fprime4_power_mW": power1_default*1e3,
    "Fprime3_eom_frequency_MHz": eom0_default/(2*np.pi*1e6),
    "Fprime4_eom_frequency_MHz": eom1_default/(2*np.pi*1e6),
    "Fprime3_modulation_depth_rad": beta0_default,
    "Fprime4_modulation_depth_rad": beta1_default,
    "polarization_switch_frequency_MHz": polarization_frequency_default/(2*np.pi*1e6),
    "polarization_scheme": "complementary X/Z switching",
    "sigma_z_mm": sigma_z_default*1e3,
    "sigma_y_mm": sigma_y_default*1e3,
    "multipass_count": len(xlocs_default),
    "multipass_relative_power": list(pass_power),
    "main_couplings": [float(abs(c.main_coupling)) for c in obe_system.couplings],
}
"""),
        md("""### F′=3 detuning × F′=4 detuning

Array shape is `[F′=3 detuning, F′=4 detuning]`. This order is retained in memory and in HDF5."""),
        code("""# Coarse axes for validation; increase resolution for the final scan.
detuning0_MHz = np.linspace(-40.0, 40.0, 41)
detuning1_MHz = np.linspace(-40.0, 40.0, 41)

detuning_scan = grid_scan_with_progress(
    prepared, rho0, (0.0, T_END),
    scan={
        "detuning0": detuning0_MHz*2*np.pi*1e6,
        "detuning1": detuning1_MHz*2*np.pi*1e6,
    },
    **scan_solver_options,
    chunk_size=32,
    description="J=2 F′=3/F′=4 detuning",
)
detuning_populations = detuning_scan.values.reshape(
    len(detuning0_MHz), len(detuning1_MHz), len(obe_system.QN)
)
detuning_remaining = detuning_populations[..., ground_idx].sum(axis=-1).real
detuning_cleanup = 1.0-detuning_remaining
print("remaining-population shape [F′=3 detuning, F′=4 detuning] =", detuning_remaining.shape)
"""),
        code("""scan_file = save_scan_hdf(
    "SPB_J2_cleanup_scans.h5",
    "Fprime3_Fprime4_detuning",
    detuning_remaining,
    axes={
        "Fprime3_detuning_MHz": (detuning0_MHz, "MHz"),
        "Fprime4_detuning_MHz": (detuning1_MHz, "MHz"),
    },
    quantity="final_population_remaining_in_X_J2",
    fixed_parameters=fixed_conditions,
    model=model_metadata,
    solver=scan_solver_options,
    notes="Two independent EOM-modulated lasers; opposite-parity B partners retained.",
    overwrite=False,  # Protect an existing long-running scan; change run_name for a new run.
)
print("saved:", scan_file)
"""),
        code("""fig, ax = plt.subplots(figsize=(8, 6))
mesh = ax.pcolormesh(
    detuning1_MHz, detuning0_MHz, 100*detuning_cleanup, shading="auto"
)
fig.colorbar(mesh, ax=ax, label="Cleanup efficiency / %")
ax.set(xlabel="F′=4 carrier detuning / MHz", ylabel="F′=3 carrier detuning / MHz")
plt.show()
"""),
        md("""### HDF5 readback example"""),
        code("""saved_run = load_scan_hdf("SPB_J2_cleanup_scans.h5", "Fprime3_Fprime4_detuning")
print(saved_run["results"].shape)
print(saved_run["units"])
print(saved_run["fixed_parameters"])
"""),
    ]
    nbf.write(nb, path)


if __name__ == "__main__":
    update_j1()
    update_j2()
    print("Updated J=1 and J=2 scan sections")
