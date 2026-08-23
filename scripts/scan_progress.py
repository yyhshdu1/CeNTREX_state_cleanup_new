"""Progress-reporting wrappers for centrex-tlf Rust parameter scans."""
from __future__ import annotations

from dataclasses import replace
from typing import Any, Mapping

import numpy as np
from tqdm.auto import tqdm

from centrex_tlf import lindblad


def grid_scan_with_progress(
    prepared,
    rho0,
    t_span,
    *,
    scan: Mapping[str, np.ndarray],
    chunk_size: int = 32,
    description: str = "OBE scan",
    **kwargs: Any,
):
    """Run a Cartesian grid scan in Rust-backed chunks with a real progress bar.

    The returned ``values`` ordering is identical to ``lindblad.grid_scan``:
    dimensions follow the insertion order of ``scan`` and flatten in C order.
    ``chunk_size`` controls the tradeoff between update frequency and call overhead;
    each chunk still uses Rust parallelism when ``parallel=True`` is supplied.
    """
    if not scan:
        raise ValueError("scan must contain at least one parameter")
    if chunk_size < 1:
        raise ValueError("chunk_size must be at least 1")

    slots = list(scan)
    axes = [np.asarray(values, dtype=np.complex128).reshape(-1) for values in scan.values()]
    if any(axis.size == 0 for axis in axes):
        raise ValueError("scan axes must be non-empty")

    # indexing='ij' plus C-order flattening reproduces grid_scan axis ordering.
    parameter_batch = np.stack(
        [grid.reshape(-1) for grid in np.meshgrid(*axes, indexing="ij")], axis=1
    )
    total = parameter_batch.shape[0]
    chunk_results = []

    with tqdm(total=total, unit="trajectory", desc=description) as progress:
        for start in range(0, total, chunk_size):
            stop = min(start + chunk_size, total)
            chunk_results.append(
                lindblad.parameter_scan(
                    prepared,
                    rho0,
                    t_span,
                    parameter_slots=slots,
                    parameter_batch=parameter_batch[start:stop],
                    **kwargs,
                )
            )
            progress.update(stop-start)

    first = chunk_results[0]
    values = np.concatenate([result.values for result in chunk_results], axis=0)
    parameter_values = np.concatenate(
        [result.parameter_values for result in chunk_results], axis=0
    )
    return replace(
        first,
        values=values,
        trajectory_count=total,
        parameter_values=parameter_values,
        solver_stats={"chunks": [result.solver_stats for result in chunk_results]},
        metadata={
            **first.metadata,
            "grid_axis_names": slots,
            "grid_axis_lengths": [int(axis.size) for axis in axes],
            "progress_chunk_size": int(chunk_size),
        },
    )
