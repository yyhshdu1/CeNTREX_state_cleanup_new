"""Structured HDF5 storage for CeNTREX parameter scans."""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Mapping

import h5py
import numpy as np


SCHEMA_VERSION = "centrex-cleanup-scan-v1"


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def save_scan_hdf(
    filename: str | Path,
    run_name: str,
    results: np.ndarray,
    *,
    axes: Mapping[str, tuple[np.ndarray, str]],
    quantity: str,
    fixed_parameters: Mapping[str, Any],
    model: Mapping[str, Any],
    solver: Mapping[str, Any],
    notes: str = "",
    overwrite: bool = False,
) -> Path:
    """Save a scan and all information required to interpret or reproduce it.

    ``axes`` is ordered and defines the dimensions of ``results``. Each entry is
    ``name: (values, unit)``. Axis arrays are saved in the units stated here.
    """
    filename = Path(filename)
    result_array = np.asarray(results)
    expected_shape = tuple(len(np.asarray(values)) for values, _ in axes.values())
    if result_array.shape != expected_shape:
        raise ValueError(
            f"results shape {result_array.shape} does not match axes {expected_shape}"
        )

    with h5py.File(filename, "a") as h5:
        runs = h5.require_group("runs")
        if run_name in runs:
            if not overwrite:
                raise FileExistsError(
                    f"run {run_name!r} already exists in {filename}; choose a new name "
                    "or pass overwrite=True"
                )
            del runs[run_name]

        group = runs.create_group(run_name)
        group.attrs["schema_version"] = SCHEMA_VERSION
        group.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
        group.attrs["quantity"] = quantity
        group.attrs["axis_order_json"] = json.dumps(list(axes))
        group.attrs["notes"] = notes
        group.attrs["fixed_parameters_json"] = json.dumps(_jsonable(fixed_parameters))
        group.attrs["model_json"] = json.dumps(_jsonable(model))
        group.attrs["solver_json"] = json.dumps(_jsonable(solver))

        axes_group = group.create_group("axes")
        for name, (values, unit) in axes.items():
            values = np.asarray(values)
            dataset = axes_group.create_dataset(name, data=values)
            dataset.attrs["unit"] = unit
            if values.size:
                dataset.attrs["minimum"] = float(np.min(values))
                dataset.attrs["maximum"] = float(np.max(values))
                dataset.attrs["points"] = int(values.size)

        result_dataset = group.create_dataset(
            "results", data=result_array, compression="gzip", shuffle=True
        )
        result_dataset.attrs["quantity"] = quantity

    return filename.resolve()


def load_scan_hdf(filename: str | Path, run_name: str) -> dict[str, Any]:
    """Load results, axes, units, and JSON metadata from one saved run."""
    with h5py.File(filename, "r") as h5:
        group = h5[f"runs/{run_name}"]
        axis_order = json.loads(group.attrs["axis_order_json"])
        return {
            "results": group["results"][...],
            "axes": {name: group[f"axes/{name}"][...] for name in axis_order},
            "units": {
                name: group[f"axes/{name}"].attrs["unit"] for name in axis_order
            },
            "quantity": group.attrs["quantity"],
            "created_utc": group.attrs["created_utc"],
            "notes": group.attrs.get("notes", ""),
            "fixed_parameters": json.loads(group.attrs["fixed_parameters_json"]),
            "model": json.loads(group.attrs["model_json"]),
            "solver": json.loads(group.attrs["solver_json"]),
        }
