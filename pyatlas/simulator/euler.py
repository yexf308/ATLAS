"""Euler-Maruyama simulators used throughout the ATLAS workflow."""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

from ..core.types import DiffusionFunction, DriftFunction, SimulationOutput, SimulationParameters


def _default_rng(rng: np.random.Generator | None) -> np.random.Generator:
    return rng if rng is not None else np.random.default_rng()


def _apply_diffusion(matrix: np.ndarray | float, noise: np.ndarray) -> np.ndarray:
    arr = np.asarray(matrix, dtype=float)
    if arr.ndim == 0:
        return noise * float(arr)
    if arr.ndim == 1:
        return arr * noise
    return arr @ noise


def _compute_statistics(samples: np.ndarray) -> tuple[list[np.ndarray], list[np.ndarray]]:
    covariances: list[np.ndarray] = []
    means: list[np.ndarray] = []
    for slice_ in samples.transpose(1, 0, 2):
        means.append(slice_.mean(axis=0))
        covariances.append(np.cov(slice_, rowvar=False))
    return covariances, means


def simulate_euler(
    drift: DriftFunction,
    diffusion: DiffusionFunction,
    x0: np.ndarray,
    params: SimulationParameters,
) -> SimulationOutput:
    rng = _default_rng(params.rng)
    dt = params.dt
    t_max = params.t_max
    n_paths = params.n_samples

    x0 = np.asarray(x0, dtype=float)
    dim = x0.size
    n_steps = int(round(t_max / dt))
    store = params.store_trajectories
    chunk = max(1, params.chunk_size)

    trajectories = (
        np.zeros((n_paths, n_steps + 1, dim), dtype=float) if store else None
    )
    sum_vectors = np.zeros((n_steps + 1, dim), dtype=float)
    sum_outer = np.zeros((n_steps + 1, dim, dim), dtype=float)

    sqrt_dt = np.sqrt(dt)
    produced = 0
    while produced < n_paths:
        batch = min(chunk, n_paths - produced)
        batch_traj = np.zeros((batch, n_steps + 1, dim), dtype=float)
        batch_traj[:, 0, :] = x0

        for path in range(batch):
            current = x0.copy()
            for step in range(n_steps):
                noise = rng.standard_normal(dim)
                drift_step = drift(current) * dt
                diffusion_matrix = diffusion(current)
                diffusion_step = _apply_diffusion(diffusion_matrix, noise) * sqrt_dt
                current = current + drift_step + diffusion_step
                batch_traj[path, step + 1, :] = current

        if store and trajectories is not None:
            trajectories[produced : produced + batch, :, :] = batch_traj

        sum_vectors += batch_traj.sum(axis=0)
        for step in range(n_steps + 1):
            data = batch_traj[:, step, :]
            sum_outer[step] += data.T @ data

        produced += batch

    means = [sum_vectors[step] / n_paths for step in range(n_steps + 1)]
    covariances: list[np.ndarray] = []
    for step in range(n_steps + 1):
        mean = means[step]
        outer = sum_outer[step]
        centered = outer - n_paths * np.outer(mean, mean)
        denom = max(n_paths - 1, 1)
        covariances.append(centered / denom)

    return SimulationOutput(trajectories=trajectories, covariances=covariances, means=means)


def simulate_euler_parallel(
    drift: DriftFunction,
    diffusion: DiffusionFunction,
    x0: np.ndarray,
    params: SimulationParameters,
) -> SimulationOutput:
    return simulate_euler(drift, diffusion, x0, params)


def simulate_nonlinear_transformed(
    drift: DriftFunction,
    diffusion: DiffusionFunction,
    forward_transform: Callable[[np.ndarray], np.ndarray],
    inverse_transform: Callable[[np.ndarray], np.ndarray],
    x0: np.ndarray,
    params: SimulationParameters,
) -> SimulationOutput:
    rng = _default_rng(params.rng)
    dt = params.dt
    t_max = params.t_max
    n_paths = params.n_samples

    x0 = np.asarray(x0, dtype=float)
    dim = x0.size
    n_steps = int(round(t_max / dt))
    trajectories = np.zeros((n_paths, n_steps + 1, dim), dtype=float)
    transformed = np.zeros_like(trajectories)

    transformed[:, 0, :] = inverse_transform(x0)
    trajectories[:, 0, :] = x0

    sqrt_dt = np.sqrt(dt)
    for path in range(n_paths):
        current = transformed[path, 0, :].copy()
        for step in range(n_steps):
            noise = rng.standard_normal(dim)
            drift_step = drift(current) * dt
            diffusion_step = diffusion(current) @ noise * sqrt_dt
            current = current + drift_step + diffusion_step
            transformed[path, step + 1, :] = current
            trajectories[path, step + 1, :] = forward_transform(current)

    covariances, means = _compute_statistics(trajectories)
    return SimulationOutput(trajectories=trajectories, covariances=covariances, means=means)
