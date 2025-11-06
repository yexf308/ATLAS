"""Projection helpers for the Python ATLAS implementation."""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

from .types import Chart, ChartList


def check_point_in_chart(
    point: np.ndarray,
    chart: Chart,
    t0: float,
    radius_factor: float,
    chi_p: float,
    threshold: float,
    *,
    fast_mode: bool = True,
) -> bool:
    point = np.asarray(point)
    coeff = chart.project_coefficients(point, fast_mode=fast_mode)
    distance = np.sum(np.square(coeff) / np.square(chart.sigma))
    within_radius = distance <= (radius_factor**2) * t0 * chi_p
    close_to_center = np.linalg.norm(point - chart.X_int) < threshold
    return bool(within_radius and close_to_center)


def chart_distance(
    point: np.ndarray,
    chart: Chart,
    t0: float,
    chi_p: float,
    threshold: float,
    mode: int,
    *,
    fast_mode: bool = True,
) -> float:
    point = np.asarray(point)
    if mode == 1 and np.linalg.norm(point - chart.X_int) > threshold:
        return 10.0
    coeff = chart.project_coefficients(point, fast_mode=fast_mode)
    t = np.sum(np.square(coeff) / np.square(chart.sigma)) / chi_p
    return float(min(np.sqrt(t / t0), 10.0))


def weighted_projection(
    point: np.ndarray,
    charts: ChartList,
    neighbors: Sequence[int],
    nearest: int,
    connectivity: Sequence[Sequence[int]],
    t0: float,
    chi_p: float,
    ambient_dim: int,
    slow_dim: int,
    threshold: Sequence[float] | float,
    option: int,
    mode: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, Sequence[int]]:
    point = np.asarray(point)
    neighbors = list(neighbors)
    threshold_value = float(threshold[0] if isinstance(threshold, (list, tuple, np.ndarray)) else threshold)

    while True:
        times = []
        for idx in neighbors:
            chart = charts[idx]
            direction = chart.WU if option == 1 else chart.U
            coeff = (point - chart.X_int) @ direction
            t = np.sum(np.square(coeff) / np.square(chart.sigma)) / chi_p
            times.append(min(np.sqrt(t / t0), 10.0))
        min_idx = int(np.argmin(times))
        candidate = neighbors[min_idx]
        if candidate == nearest:
            break
        nearest = candidate
        neighbors = list(connectivity[nearest])

    Lambda = np.zeros((ambient_dim, ambient_dim))
    b = np.zeros(ambient_dim)
    weight = 0.0
    X_accum = np.zeros((len(neighbors), ambient_dim))
    times = np.zeros(len(neighbors))

    for k, idx in enumerate(neighbors):
        chart = charts[idx]
        center = chart.X_int
        direction = chart.WU if option == 1 else chart.U
        coeff = (point - center) @ direction
        t = np.sum(np.square(coeff) / np.square(chart.sigma)) / chi_p
        t = min(np.sqrt(t / t0), 10.0)
        if mode == 1 and np.linalg.norm(point - center) > threshold_value:
            times[k] = 10.0
            X_accum[k, :] = 0.0
            continue
        times[k] = t
        weight_k = np.exp(-t)
        projected = (coeff @ direction.T) + center
        X_accum[k, :] = projected * weight_k
        Lambda += chart.Lambda * weight_k
        b += chart.b * weight_k
        weight += weight_k

    if weight > 0:
        X_proj = np.sum(X_accum, axis=0) / weight
        Lambda /= weight
        b /= weight
    else:
        X_proj = point.copy()

    U, S, _ = np.linalg.svd(Lambda)
    U = U[:, :slow_dim]
    S = S[:slow_dim]
    Lambda_hat = U @ np.diag(S) @ U.T
    H_hat = U @ np.diag(np.sqrt(S))

    return X_proj, b, Lambda_hat, H_hat, times, nearest, neighbors
