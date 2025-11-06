"""Peanut model configuration translated from the MATLAB reference implementation."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from ..core.types import ModelConfiguration, ModelDynamics, SimulationParameters


@dataclass
class PeanutParameters:
    ambient_dim: int = 3
    slow_dim: int = 2
    chi_p: float = 5.991
    num_initial: int = 100
    num_short_trajectories: int = 400_000
    epsilon: float = 0.005
    dt: float = 5e-4
    t0: float = 0.1
    T_max: float = 0.1
    T_one: float = 1_000.0

    threshold: tuple[float, float] = (1.0, 0.001)
    connectivity_threshold: float = 3.0
    explore_threshold: float = 1.05
    option: int = 1
    modify: int = 0
    relearn_option: int = 0
    lower_bound: float = 0.05
    upper_bound: float = 0.10

    a1: float = 4.0
    a2: float = 8.0
    c1: float = 2.0
    c2: float = 0.5
    c3: float = 0.05
    c4: float = 0.4
    c5: float = 0.05
    c6: float = 0.4
def _drift(parameters: PeanutParameters) -> Callable[[np.ndarray], np.ndarray]:
    eps = parameters.epsilon
    a1 = parameters.a1
    a2 = parameters.a2
    c1 = parameters.c1
    c3 = parameters.c3
    c5 = parameters.c5
    c4 = parameters.c4
    c6 = parameters.c6

    def drift(state: np.ndarray) -> np.ndarray:
        x, y, z = np.asarray(state, dtype=float)
        r = np.linalg.norm(state)
        rsin = np.linalg.norm(state[:2])
        r_sq = max(r * r, 1e-12)
        rsin_sq = max(rsin * rsin, 1e-12)

        J = np.array(
            [
                [x / r, x * z / rsin if rsin > 0 else 0.0, -y],
                [y / r, y * z / rsin if rsin > 0 else 0.0, x],
                [z / r, -rsin, 0.0],
            ]
        )

        b = np.array(
            [
                -c1 / eps * (1 - np.sqrt(a1 + a2 * z**2 / r_sq) / r),
                c3 * (4 * z**3 / (r**3) - 3 * z / r) / rsin,
                c5 * (y * z / (rsin * r_sq) + x / r_sq),
            ]
        )

        d1 = c4**2 * rsin_sq / (r_sq**2)
        d2 = c6**2 / r_sq
        ito = 0.5 * np.array([-x * (d1 + d2), -y * (d1 + d2), -z * d1])

        return J @ b + ito

    return drift


def _diffusion(parameters: PeanutParameters) -> Callable[[np.ndarray], np.ndarray]:
    eps = parameters.epsilon
    c2 = parameters.c2
    c4 = parameters.c4
    c6 = parameters.c6

    def diffusion(state: np.ndarray) -> np.ndarray:
        x, y, z = np.asarray(state, dtype=float)
        r = np.linalg.norm(state)
        r_sq = max(r * r, 1e-12)
        rsin = np.linalg.norm([x, y])
        rsin_sq = max(rsin * rsin, 1e-12)
        sqrt_eps = np.sqrt(eps)

        d1 = c2 / (sqrt_eps * r_sq)
        d2 = c4 / r_sq
        d3 = c6 / max(r, 1e-12)

        return np.array(
            [
                [d1 * x, d2 * x * z, -y * d3],
                [d1 * y, d2 * y * z, x * d3],
                [d1 * z, -c4 * rsin_sq / r_sq, 0.0],
            ]
        )

    return diffusion


def peanut_configuration(datapath: str | Path | None = None) -> ModelConfiguration:
    params = PeanutParameters()

    drift = _drift(params)
    diffusion = _diffusion(params)

    theta = np.random.rand() * np.pi
    phi = np.random.rand() * 2 * np.pi
    R = np.sqrt(params.a1 + params.a2 * (np.cos(theta) ** 2))
    x_int = np.array([
        R * np.sin(theta) * np.cos(phi),
        R * np.sin(theta) * np.sin(phi),
        R * np.cos(theta),
    ])

    learning = {
        "t0": params.t0,
        "threshold": np.array(params.threshold),
        "connectivity_threshold": params.connectivity_threshold,
        "explore_threshold": params.explore_threshold,
        "option": params.option,
        "modify": params.modify,
        "relearn_option": params.relearn_option,
        "lower_bound": params.lower_bound,
        "upper_bound": params.upper_bound,
        "chi_p": params.chi_p,
        "ambient_dim": params.ambient_dim,
        "slow_dim": params.slow_dim,
    }

    simulation = {
        "long_trajectory_steps": params.T_one,
        "short_simulation": SimulationParameters(
            dt=params.dt,
            t_max=params.T_max,
            n_samples=params.num_short_trajectories,
            chunk_size=5_000,
            store_trajectories=False,
        ),
    }

    file_config = {}
    if datapath is not None:
        base = Path(datapath)
        file_config = {
            "chart": base / "chart.mat",
            "chart_part": base / "chart_part.mat",
            "tran": base / "TranM.mat",
            "fpt": base / "Peanut_FPT.mat",
        }

    return ModelConfiguration(
        parameters={
            "ambient_dim": params.ambient_dim,
            "slow_dim": params.slow_dim,
            "chi_p": params.chi_p,
            "threshold": np.array(params.threshold),
            "connectivity_threshold": params.connectivity_threshold,
            "explore_threshold": params.explore_threshold,
            "option": params.option,
            "modify": params.modify,
            "relearn_option": params.relearn_option,
            "lower_bound": params.lower_bound,
            "upper_bound": params.upper_bound,
            "num_initial": params.num_initial,
            "long_trajectory_steps": params.T_one,
            "num_short_trajectories": params.num_short_trajectories,
            "X_int": x_int,
        }
        | file_config,
        dynamics=ModelDynamics(drift=drift, diffusion=diffusion),
        learning=learning,
        simulation=simulation,
    )
