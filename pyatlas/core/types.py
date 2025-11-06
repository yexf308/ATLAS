"""Core data structures for the Python ATLAS implementation."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional

import numpy as np


@dataclass
class Chart:
    """Representation of a learned local chart on the slow manifold."""

    U: np.ndarray
    sigma: np.ndarray
    Lambda: np.ndarray
    b: np.ndarray
    X_int: np.ndarray
    WU: np.ndarray
    V: np.ndarray | None = None
    sigma_fast: np.ndarray | None = None
    learning_time: float | None = None

    def project_coefficients(self, point: np.ndarray, *, fast_mode: bool = True) -> np.ndarray:
        """Return coordinates of ``point`` in the chart basis."""
        center = self.X_int
        direction = self.WU if fast_mode else self.U
        return (np.asarray(point) - center) @ direction

    def mahalanobis_time(self, point: np.ndarray, chi_p: float) -> float:
        coeff = self.project_coefficients(point)
        sigma_sq = np.square(self.sigma)
        return float(np.sum(np.square(coeff) / sigma_sq) / chi_p)


ChartList = List[Chart]


@dataclass
class SimulationOutput:
    trajectories: np.ndarray | None
    covariances: List[np.ndarray]
    means: List[np.ndarray]


DriftFunction = Callable[[np.ndarray], np.ndarray]
DiffusionFunction = Callable[[np.ndarray], np.ndarray]


@dataclass
class SimulationParameters:
    dt: float
    t_max: float
    n_samples: int = 1
    rng: np.random.Generator | None = None
    chunk_size: int = 1
    store_trajectories: bool = True


@dataclass
class ModelDynamics:
    drift: DriftFunction
    diffusion: DiffusionFunction
    nonlinear_forward: Optional[Callable[[np.ndarray], np.ndarray]] = None
    nonlinear_inverse: Optional[Callable[[np.ndarray], np.ndarray]] = None


@dataclass
class ModelConfiguration:
    parameters: dict
    dynamics: ModelDynamics
    learning: dict = field(default_factory=dict)
    simulation: dict = field(default_factory=dict)
