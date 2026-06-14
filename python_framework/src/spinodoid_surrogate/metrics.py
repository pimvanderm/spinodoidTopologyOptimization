"""Evaluation metrics matching the MATLAB R^2 computation."""

from __future__ import annotations

import numpy as np


def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    """Coefficient of determination, computed per output column.

    Mirrors the MATLAB loop::

        SS_res = sum((y - yhat).^2);
        SS_tot = sum((y - mean(y)).^2);
        R2 = 1 - SS_res / SS_tot;

    Parameters
    ----------
    y_true, y_pred:
        Arrays of shape ``(n_samples, n_outputs)`` (or 1-D for a single output).

    Returns
    -------
    np.ndarray
        One R^2 value per output column.
    """
    y_true = np.asarray(y_true, dtype=np.float64)
    y_pred = np.asarray(y_pred, dtype=np.float64)
    if y_true.ndim == 1:
        y_true = y_true[:, None]
        y_pred = y_pred[:, None]
    ss_res = np.sum((y_true - y_pred) ** 2, axis=0)
    ss_tot = np.sum((y_true - y_true.mean(axis=0)) ** 2, axis=0)
    # Avoid divide-by-zero for constant columns: R^2 is 0 if residual is 0, else -inf.
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = 1.0 - ss_res / ss_tot
    r2 = np.where(ss_tot == 0.0, np.where(ss_res == 0.0, 1.0, -np.inf), r2)
    return r2
