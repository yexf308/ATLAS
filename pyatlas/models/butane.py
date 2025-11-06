"""Butane model configuration translated from the MATLAB reference implementation."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from ..core.types import ModelConfiguration, ModelDynamics, SimulationParameters


@dataclass
class ButaneParameters:
    ambient_dim: int = 6
    slow_dim: int = 1
    chi_p: float = 1.96
    num_initial: int = 2
    num_short_trajectories: int = 1_600_000
    dt: float = 1e-6
    t0: float = 10e-6
    T_max: float = 50e-6
    T_one: float = 10.0

    threshold: tuple[float, float] = (1.5, 0.001)
    connectivity_threshold: float = 4.0
    explore_threshold: float = 0.95
    option: int = 1
    modify: int = 0
    relearn_option: int = 1
    lower_bound: float = 40e-6
    upper_bound: float = 50e-6

    l: float = 1.53
    k2: float = 319_225.0
    k3: float = 62_500.0
    theta: float = 112 / 180 * np.pi
    c1: float = 2037.82
    c2: float = 158.52
    c3: float = -3227.70
    beta: float = 4e-3

    @property
    def sigma(self) -> float:
        return np.sqrt(2 / self.beta)


def _drift(parameters: ButaneParameters) -> Callable[[np.ndarray], np.ndarray]:
    l = parameters.l
    k2 = parameters.k2
    k3 = parameters.k3
    theta = parameters.theta
    c1 = parameters.c1
    c2 = parameters.c2
    c3 = parameters.c3

    def drift(state: np.ndarray) -> np.ndarray:
        x1, y1, y3, x4, y4, z4 = np.asarray(state, dtype=float)
        drift_vec = np.zeros(6, dtype=float)

        r1 = np.sqrt(x1**2 + y1**2) + 1e-12
        r4 = np.sqrt(x4**2 + z4**2 + (y3 - y4) ** 2) + 1e-12

        b1 = 1 - l / r1
        b3 = 1 - l / r4

        theta1 = theta - np.arccos(np.clip(y1 * np.sign(y3) / r1, -1.0, 1.0))
        theta2 = theta - np.arccos(
            np.clip((y3 - y4) * np.sign(y3) / r4, -1.0, 1.0)
        )

        denom = (x4**2 + z4**2) ** (2.5) + 1e-12
        c_term = z4 * (
            (c1 * (x4**2 + z4**2) + 3 * c3 * x4**2) + 2 * c2 * x4 * np.sqrt(x4**2 + z4**2)
        ) / denom

        drift_vec[0] = k2 * x1 * b1 - k3 * y1 * np.sign(x1 * y3) / (r1**2) * theta1
        drift_vec[1] = k2 * y1 * b1 + k3 * np.abs(x1) * np.sign(y3) / (r1**2) * theta1
        drift_vec[2] = k2 * ((y3 - y4) * b3 + (-l + np.abs(y3)) * np.sign(y3)) + k3 * np.sqrt(x4**2 + z4**2) * np.sign(y3) / (r4**2) * theta2
        drift_vec[3] = (
            k2 * x4 * b3
            - k3 * x4 * (y3 - y4) / np.sqrt(x4**2 + z4**2) * np.sign(y3) / (r4**2) * theta2
            + z4 * c_term
        )
        drift_vec[4] = -k2 * (y3 - y4) * b3 - k3 * np.sqrt(x4**2 + z4**2) * np.sign(y3) / (r4**2) * theta2
        drift_vec[5] = (
            k2 * z4 * b3
            - k3 * z4 * (y3 - y4) / np.sqrt(x4**2 + z4**2) * np.sign(y3) / (r4**2) * theta2
            - x4 * c_term
        )

        return -drift_vec

    return drift


def _diffusion(parameters: ButaneParameters) -> Callable[[np.ndarray], np.ndarray]:
    sigma = parameters.sigma

    def diffusion(_: np.ndarray) -> np.ndarray:
        return float(sigma)

    return diffusion


def butane_configuration(datapath: str | Path | None = None) -> ModelConfiguration:
    params = ButaneParameters()
    drift = _drift(params)
    diffusion = _diffusion(params)

    x_int = np.array([-1.4461, 0.5578, 1.5361, 1.4268, -2.0593, 0.1615])

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
            chunk_size=20_000,
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
            "fpt": base / "Butane_FPT.mat",
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
