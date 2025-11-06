"""High-level learning workflows that orchestrate simulator and chart updates."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

import numpy as np

from ..core.learning import learn_slow_manifold
from ..core.projection import weighted_projection
from ..core.types import Chart, ModelConfiguration, SimulationParameters
from ..simulator.euler import simulate_euler


@dataclass
class AtlasState:
    charts: List[Chart]
    connectivity: List[List[int]]
    nearest: int


def _sample_long_trajectory(
    config: ModelConfiguration,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    sim_params = SimulationParameters(
        dt=config.simulation["short_simulation"].dt,
        t_max=config.simulation["long_trajectory_steps"] * config.simulation["short_simulation"].dt,
        n_samples=1,
        rng=rng,
    )
    output = simulate_euler(
        config.dynamics.drift,
        config.dynamics.diffusion,
        np.asarray(config.parameters.get("X_int", np.zeros(config.parameters["ambient_dim"]))),
        sim_params,
    )
    if output.trajectories is None:
        raise ValueError("Long trajectory simulation must store trajectories")
    return output.trajectories[0]


def _evenly_sample(points: np.ndarray, count: int, burn_in: int = 0) -> np.ndarray:
    if burn_in > 0:
        points = points[burn_in:]
    if len(points) < count:
        raise ValueError("Not enough points to sample from trajectory")
    indices = np.linspace(0, len(points) - 1, count, dtype=int)
    return points[indices]


def _learn_chart_for_point(
    config: ModelConfiguration,
    point: np.ndarray,
    rng: np.random.Generator | None = None,
) -> Chart:
    sim_params: SimulationParameters = config.simulation["short_simulation"]
    params = SimulationParameters(
        dt=sim_params.dt,
        t_max=sim_params.t_max,
        n_samples=sim_params.n_samples,
        rng=rng,
        chunk_size=sim_params.chunk_size,
        store_trajectories=False,
    )
    output = simulate_euler(
        config.dynamics.drift,
        config.dynamics.diffusion,
        point,
        params,
    )
    data = {
        "dt": sim_params.dt,
        "LowerBound": config.learning["lower_bound"],
        "UpperBound": config.learning["upper_bound"],
    }
    return learn_slow_manifold(
        data,
        output.covariances,
        output.means,
        config.learning["ambient_dim"],
        config.learning["slow_dim"],
        modify=bool(config.parameters.get("modify", 0)),
    )


def learn_initial_charts(
    config: ModelConfiguration,
    *,
    rng: np.random.Generator | None = None,
    burn_in: int = 1000,
) -> List[Chart]:
    trajectory = _sample_long_trajectory(config, rng=rng)
    samples = _evenly_sample(trajectory, config.parameters["num_initial"], burn_in=burn_in)
    charts = []
    for point in samples:
        chart = _learn_chart_for_point(config, point, rng=rng)
        chart.X_int = point
        charts.append(chart)
    return charts


def build_connectivity(
    charts: Sequence[Chart],
    threshold: float,
) -> List[List[int]]:
    centers = np.array([chart.X_int for chart in charts])
    distances = np.linalg.norm(centers[:, None, :] - centers[None, :, :], axis=-1)
    connectivity: List[List[int]] = []
    for i in range(len(charts)):
        neighbors = np.where(distances[i] <= threshold)[0].tolist()
        if i not in neighbors:
            neighbors.append(i)
        connectivity.append(sorted(set(neighbors)))
    return connectivity


def initialise_atlas(
    config: ModelConfiguration,
    *,
    rng: np.random.Generator | None = None,
) -> AtlasState:
    charts = learn_initial_charts(config, rng=rng)
    threshold = config.learning["connectivity_threshold"] * np.sqrt(config.learning["t0"])
    connectivity = build_connectivity(charts, threshold)
    return AtlasState(charts=charts, connectivity=connectivity, nearest=0)


def project_point(
    state: AtlasState,
    point: np.ndarray,
    config: ModelConfiguration,
    *,
    option: int | None = None,
    mode: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    option_value = option if option is not None else config.parameters.get("option", 1)
    chi_p = config.learning["chi_p"]
    t0 = config.learning["t0"]
    threshold = config.parameters["threshold"]
    neighbors = state.connectivity[state.nearest]
    X_proj, b, Lambda_hat, _, _, nearest, neigh = weighted_projection(
        point,
        state.charts,
        neighbors,
        state.nearest,
        state.connectivity,
        t0,
        chi_p,
        config.parameters["ambient_dim"],
        config.parameters["slow_dim"],
        threshold,
        option_value,
        mode,
    )
    state.nearest = nearest
    state.connectivity[state.nearest] = list(neigh)
    return X_proj, b, Lambda_hat
