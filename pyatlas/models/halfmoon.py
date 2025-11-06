"""Halfmoon model configuration for the Python ATLAS implementation."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import numpy as np

from ..core.types import ModelConfiguration, ModelDynamics, SimulationParameters


@dataclass
class HalfmoonParameters:
    ambient_dim: int = 20
    slow_dim: int = 1
    chi_p: float = 1.96
    num_initial: int = 10
    num_short_trajectories: int = 8_000_000
    dt: float = 0.05
    long_trajectory_time: int = 800_000

    @property
    def t0(self) -> float:
        return 20 * self.dt

    @property
    def short_time(self) -> float:
        return 25 * self.dt

    @property
    def lower_bound(self) -> float:
        return 20 * self.dt

    @property
    def upper_bound(self) -> float:
        return 25 * self.dt

    @property
    def threshold(self) -> np.ndarray:
        return np.array([0.5, 0.01])

    @property
    def connectivity_threshold(self) -> float:
        return 4.0

    @property
    def explore_threshold(self) -> float:
        return 0.95

    @property
    def option(self) -> int:
        return 1

    @property
    def modify(self) -> int:
        return 0

    @property
    def relearn_option(self) -> int:
        return 1


def _build_transforms(dim: int) -> tuple[np.ndarray, np.ndarray]:
    tran = np.diag(np.ones(dim - 2), 1)
    tran[:, 0] = 1.0
    tran = tran[: dim - 2, :]

    tran_inv = np.diag(np.ones(dim - 2), 1)
    tran_inv[:, 0] = -1.0
    tran_inv = tran_inv[: dim - 2, :]
    return tran, tran_inv


def _drift_factory(parameters: Dict[str, float]) -> callable:
    a1 = parameters["a1"]
    a2 = parameters["a2"]
    a3 = parameters["a3"]
    b1 = parameters["b1"]
    b3 = parameters["b3"]

    def drift(state: np.ndarray) -> np.ndarray:
        state = np.asarray(state, dtype=float)
        theta = state[0]
        r = state[1]
        r_fast = state[2:]
        result = np.zeros_like(state)
        result[0] = a1 + a2 * np.sin(2 * theta) + a3 * np.cos(theta)
        result[1] = b1 * (1 - r)
        result[2:] = -b3 * r_fast
        return result

    return drift


def _diffusion_factory(dim: int, parameters: Dict[str, float]) -> callable:
    a4 = parameters["a4"]
    b2 = parameters["b2"]
    b4 = parameters["b4"]

    def diffusion(state: np.ndarray) -> np.ndarray:
        _ = state
        diag = np.full(dim, b4)
        diag[0] = a4
        diag[1] = b2
        return np.diag(diag)

    return diffusion


def _nonlinear_transform_factory(parameters: Dict[str, float], tran: np.ndarray) -> callable:
    t = parameters["t"]

    def transform(state: np.ndarray) -> np.ndarray:
        arr = np.asarray(state, dtype=float)
        if arr.ndim == 1:
            theta = arr[0]
            r = arr[1]
            result = np.zeros_like(arr)
            result[0] = r * np.cos(theta + r + t)
            result[1] = r * np.sin(theta + r + t)
            result[2:] = tran @ arr[1:]
            return result
        if arr.ndim == 2:
            theta = arr[0, :]
            r = arr[1, :]
            result = np.zeros_like(arr)
            result[0, :] = r * np.cos(theta + r + t)
            result[1, :] = r * np.sin(theta + r + t)
            result[2:, :] = tran @ arr[1:, :]
            return result
        raise ValueError("Unsupported shape for transform")

    return transform


def _nonlinear_inverse_factory(parameters: Dict[str, float], tran_inv: np.ndarray) -> callable:
    t = parameters["t"]

    def inverse(state: np.ndarray) -> np.ndarray:
        arr = np.asarray(state, dtype=float)
        if arr.ndim == 1:
            x, y = arr[0], arr[1]
            rest = arr[2:]
            angle = np.mod(np.arctan2(y, x), 2 * np.pi)
            r = np.sqrt(x**2 + y**2)
            theta = angle - r - t
            if theta > np.pi:
                theta -= 2 * np.pi
            result = np.zeros_like(arr)
            result[0] = theta
            result[1] = r
            result[2:] = tran_inv @ np.concatenate(([r], rest))
            return result
        if arr.ndim == 2:
            x = arr[0, :]
            y = arr[1, :]
            rest = arr[2:, :]
            angle = np.mod(np.arctan2(y, x), 2 * np.pi)
            r = np.sqrt(x**2 + y**2)
            theta = angle - r - t
            theta = np.where(theta > np.pi, theta - 2 * np.pi, theta)
            stacked = np.vstack((r, rest))
            fast = tran_inv @ stacked
            result = np.vstack((theta, r, fast))
            return result
        raise ValueError("Unsupported shape for transform inverse")

    return inverse


def halfmoon_configuration(datapath: str | Path | None = None) -> ModelConfiguration:
    params = HalfmoonParameters()
    dim = params.ambient_dim

    tran, tran_inv = _build_transforms(dim)
    parameter_values = {
        "a1": 0.0,
        "a2": 5e-3,
        "a3": 2.5e-3,
        "a4": 0.06,
        "eps": 0.01,
        "b1": 0.04 / 0.01,
        "b2": 0.035 / np.sqrt(0.01),
        "b3": 0.05 / 0.01,
        "b4": 0.02 / np.sqrt(0.01),
        "t": -1.0,
        "Tran": tran,
        "Tran_inv": tran_inv,
    }

    drift = _drift_factory(parameter_values)
    diffusion = _diffusion_factory(dim, parameter_values)
    nonlin = _nonlinear_transform_factory(parameter_values, tran)
    nonlin_inv = _nonlinear_inverse_factory(parameter_values, tran_inv)

    x_int = np.ones(dim)
    x_int[1] = 0.0

    learning_config = {
        "t0": params.t0,
        "threshold": params.threshold,
        "connectivity_threshold": params.connectivity_threshold,
        "explore_threshold": params.explore_threshold,
        "option": params.option,
        "modify": params.modify,
        "relearn_option": params.relearn_option,
        "lower_bound": params.lower_bound,
        "upper_bound": params.upper_bound,
        "chi_p": params.chi_p,
        "slow_dim": params.slow_dim,
        "ambient_dim": params.ambient_dim,
    }

    simulation_config = {
        "long_trajectory_steps": params.long_trajectory_time,
        "short_simulation": SimulationParameters(
            dt=params.dt,
            t_max=params.short_time,
            n_samples=params.num_short_trajectories,
            chunk_size=10_000,
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
            "well": base / "well.mat",
            "fpt": base / "Halfmoon_FPT.mat",
        }

    return ModelConfiguration(
        parameters={
            "ambient_dim": params.ambient_dim,
            "slow_dim": params.slow_dim,
            "chi_p": params.chi_p,
            "threshold": params.threshold,
            "connectivity_threshold": params.connectivity_threshold,
            "explore_threshold": params.explore_threshold,
            "option": params.option,
            "modify": params.modify,
            "relearn_option": params.relearn_option,
            "lower_bound": params.lower_bound,
            "upper_bound": params.upper_bound,
            "num_initial": params.num_initial,
            "long_trajectory_steps": params.long_trajectory_time,
            "num_short_trajectories": params.num_short_trajectories,
            "X_int": x_int,
        }
        | file_config,
        dynamics=ModelDynamics(
            drift=drift,
            diffusion=diffusion,
            nonlinear_forward=nonlin,
            nonlinear_inverse=nonlin_inv,
        ),
        learning=learning_config,
        simulation=simulation_config,
    )
