"""Learning routines for the Python ATLAS implementation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .types import Chart


@dataclass
class LearningWindow:
    dt: float
    lower_bound: float
    upper_bound: float

    @classmethod
    def from_mapping(cls, mapping: dict) -> "LearningWindow":
        return cls(
            dt=float(mapping["dt"]),
            lower_bound=float(mapping["LowerBound"]),
            upper_bound=float(mapping["UpperBound"]),
        )

    def span_indices(self) -> np.ndarray:
        start = int(round(self.lower_bound / self.dt))
        stop = int(round(self.upper_bound / self.dt))
        return np.arange(start, stop + 1)

    def span_times(self) -> np.ndarray:
        indices = self.span_indices()
        times = indices * self.dt
        return times - times.mean()


def _principal_angle(U: np.ndarray, vector: np.ndarray) -> float:
    vector = np.asarray(vector)
    projection = U @ (U.T @ vector)
    if np.allclose(projection, 0):
        return np.pi / 2
    cos_theta = np.clip(
        np.linalg.norm(projection) / np.linalg.norm(vector), -1.0, 1.0
    )
    return float(np.arccos(cos_theta))


def learn_slow_manifold(
    data: dict,
    covariances: Sequence[np.ndarray],
    means: Sequence[np.ndarray],
    ambient_dim: int,
    slow_dim: int,
    modify: bool = False,
) -> Chart:
    window = LearningWindow.from_mapping(data)
    indices = window.span_indices()
    times = window.span_times()

    cov_span = [np.asarray(covariances[i]) for i in indices]
    mean_span = [np.asarray(means[i]) for i in indices]

    sum_cov = np.zeros((ambient_dim, ambient_dim))
    sum_mean = np.zeros(ambient_dim)

    for cov in cov_span:
        sum_cov += cov
    for mean in mean_span:
        sum_mean += mean

    count = len(indices)
    mean_cov = sum_cov / count
    mean_mean = sum_mean / count

    temp_cov = np.zeros_like(mean_cov)
    temp_mean = np.zeros_like(mean_mean)

    for cov, mean, t in zip(cov_span, mean_span, times):
        centered_cov = cov - mean_cov
        centered_mean = mean - mean_mean
        temp_cov += centered_cov * t
        temp_mean += centered_mean * t

    denom = np.linalg.norm(times) ** 2
    lambda_full = temp_cov / denom
    b_full = temp_mean / denom

    midpoint = (window.upper_bound + window.lower_bound) / 2.0
    if modify:
        x_int = mean_mean - b_full * (window.upper_bound - window.lower_bound) / 2.0
    else:
        x_int = mean_mean - b_full * midpoint

    gamma_full = mean_cov - lambda_full * midpoint

    U, S, _ = np.linalg.svd(lambda_full, full_matrices=False)
    U_slow = U[:, :slow_dim]
    S_slow = S[:slow_dim]
    lambda_hat = U_slow @ np.diag(S_slow) @ U_slow.T
    sigma = np.sqrt(S_slow)

    V, S_fast, _ = np.linalg.svd(gamma_full, full_matrices=False)
    diag_vs = S_fast
    candidate_indices = np.where(diag_vs > diag_vs[0] / 3.0)[0]
    fast_vectors = []
    for idx in candidate_indices:
        vector = V[:, idx]
        if _principal_angle(U_slow, vector) / np.pi > 0.2:
            fast_vectors.append(vector)

    if fast_vectors:
        V_index = np.column_stack(fast_vectors)
        E = np.concatenate((U_slow, V_index), axis=1)
    else:
        V_index = np.empty((ambient_dim, 0))
        E = U_slow

    Q, R = np.linalg.qr(E, mode="reduced")
    WU = Q @ np.linalg.pinv(R @ R.T) @ Q.T @ U_slow

    return Chart(
        U=U_slow,
        sigma=sigma,
        Lambda=lambda_hat,
        b=b_full.reshape(-1),
        X_int=x_int.reshape(-1),
        WU=WU,
        V=V_index,
        sigma_fast=S_fast,
        learning_time=None,
    )
