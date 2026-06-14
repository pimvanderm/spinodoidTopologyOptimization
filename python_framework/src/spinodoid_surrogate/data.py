"""Data loading and normalization for the spinodoid surrogate models.

Replicates the normalization scheme from ``trainNeuralNetworks.m``:

* **Inputs** ``[density, theta1, theta2, theta3]`` are min-max scaled to ``[0, 1]``
  using *fixed* physical bounds ``min = [0, 0, 0, 0]``, ``max = [1, 90, 90, 90]``.
* **Outputs** (the nine stiffness components) are min-max scaled to ``[0, 1]``
  using per-column min/max estimated from the *training* set.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np
import torch
from torch.utils.data import Dataset

INPUT_COLUMNS = ("density", "theta1", "theta2", "theta3")
OUTPUT_COLUMNS = ("C11", "C12", "C13", "C22", "C23", "C33", "C44", "C55", "C66")
ALL_COLUMNS = INPUT_COLUMNS + OUTPUT_COLUMNS

N_INPUTS = len(INPUT_COLUMNS)
N_OUTPUTS = len(OUTPUT_COLUMNS)

# Input and output normalization
INPUT_MIN = np.array([0.1, 15, 15, 15], dtype=np.float64)
INPUT_MAX = np.array([0.9, 90.0, 90.0, 90.0], dtype=np.float64)

PathLike = Union[str, Path]


def load_csv(path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Load a spinodoid CSV into ``(theta_raw, stiffness_raw)`` float64 arrays.

    The CSV is expected to have a header row matching :data:`ALL_COLUMNS`; the
    first four columns are design parameters and the remaining nine are the
    stiffness components.
    """
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, dtype=np.float64)
    if data.ndim == 1:  # single-row file
        data = data.reshape(1, -1)
    if data.shape[1] != len(ALL_COLUMNS):
        raise ValueError(
            f"{path} has {data.shape[1]} columns; expected {len(ALL_COLUMNS)} "
            f"({ALL_COLUMNS}). Found header: {header}"
        )
    theta_raw = data[:, :N_INPUTS]
    stiffness_raw = data[:, N_INPUTS:]
    return theta_raw, stiffness_raw


@dataclass
class Normalizer:
    """Min-max normalizer for inputs and outputs, plus optional input z-score.

    Parameters
    ----------
    input_min, input_max:
        Bounds used to map design parameters to ``[0, 1]``. Default to the fixed
        physical bounds from the MATLAB script.
    output_min, output_max:
        Per-column stiffness bounds; populated by :meth:`fit`.
    input_zscore_mean, input_zscore_std:
        Optional standardization of the *already min-max normalized* inputs,
        emulating MATLAB's ``featureInputLayer(..., "Normalization", "zscore")``.
        These are folded into :class:`~spinodoid_surrogate.models.ForwardNet`
        rather than applied here, but are stored for completeness.
    """

    input_min: np.ndarray = None  # type: ignore[assignment]
    input_max: np.ndarray = None  # type: ignore[assignment]
    output_min: Optional[np.ndarray] = None
    output_max: Optional[np.ndarray] = None
    input_zscore_mean: Optional[np.ndarray] = None
    input_zscore_std: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        if self.input_min is None:
            self.input_min = INPUT_MIN.copy()
        if self.input_max is None:
            self.input_max = INPUT_MAX.copy()
        self.input_min = np.asarray(self.input_min, dtype=np.float64)
        self.input_max = np.asarray(self.input_max, dtype=np.float64)

    # -- fitting ----------------------------------------------------------
    @classmethod
    def fit(cls, theta_raw: np.ndarray, stiffness_raw: np.ndarray) -> "Normalizer":
        """Fit output bounds (and z-score stats) from raw *training* arrays."""
        norm = cls()
        norm.output_min = stiffness_raw.min(axis=0)
        norm.output_max = stiffness_raw.max(axis=0)
        theta_mm = norm.normalize_inputs(theta_raw)
        norm.input_zscore_mean = theta_mm.mean(axis=0)
        std = theta_mm.std(axis=0)
        # Guard against zero-variance columns (matches MATLAB's safe divide).
        norm.input_zscore_std = np.where(std == 0.0, 1.0, std)
        return norm

    def _check_fitted(self) -> None:
        if self.output_min is None or self.output_max is None:
            raise RuntimeError("Normalizer.fit(...) must be called before use.")

    # -- inputs -----------------------------------------------------------
    def normalize_inputs(self, theta_raw: np.ndarray) -> np.ndarray:
        return (np.asarray(theta_raw, dtype=np.float64) - self.input_min) / (
            self.input_max - self.input_min
        )

    def denormalize_inputs(self, theta_norm: np.ndarray) -> np.ndarray:
        return np.asarray(theta_norm, dtype=np.float64) * (
            self.input_max - self.input_min
        ) + self.input_min

    # -- outputs ----------------------------------------------------------
    def normalize_outputs(self, stiffness_raw: np.ndarray) -> np.ndarray:
        self._check_fitted()
        return (np.asarray(stiffness_raw, dtype=np.float64) - self.output_min) / (
            self.output_max - self.output_min
        )

    def denormalize_outputs(self, stiffness_norm: np.ndarray) -> np.ndarray:
        self._check_fitted()
        return np.asarray(stiffness_norm, dtype=np.float64) * (
            self.output_max - self.output_min
        ) + self.output_min


class SpinodoidDataset(Dataset):
    """A ``torch`` dataset yielding ``(theta_norm, stiffness_norm)`` float32 pairs."""

    def __init__(
        self,
        theta_raw: np.ndarray,
        stiffness_raw: np.ndarray,
        normalizer: Normalizer,
    ) -> None:
        theta_norm = normalizer.normalize_inputs(theta_raw)
        stiffness_norm = normalizer.normalize_outputs(stiffness_raw)
        self.theta = torch.as_tensor(theta_norm, dtype=torch.float32)
        self.stiffness = torch.as_tensor(stiffness_norm, dtype=torch.float32)

    def __len__(self) -> int:
        return self.theta.shape[0]

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        return self.theta[idx], self.stiffness[idx]
