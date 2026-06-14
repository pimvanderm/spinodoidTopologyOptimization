"""PyTorch surrogate models for inverse-designed spinodoid metamaterials.

A faithful port of the MATLAB ``trainNeuralNetworks.m`` reference implementation
based on Kumar et al. (2020), *npj Computational Materials*
(https://doi.org/10.1038/s41524-020-0341-6).

The package exposes:

* :class:`~spinodoid_surrogate.data.Normalizer` -- min-max (and optional z-score)
  normalization matching the MATLAB pipeline.
* :class:`~spinodoid_surrogate.data.SpinodoidDataset` -- a ``torch`` dataset.
* :class:`~spinodoid_surrogate.models.ForwardNet` -- structure -> stiffness (4 -> 9).
* :class:`~spinodoid_surrogate.models.InverseNet` -- stiffness -> structure (9 -> 4).
* :func:`~spinodoid_surrogate.train_surrogate.train_full_pipeline` -- training entry point.
"""

from .config import ForwardConfig, InverseConfig
from .data import (
    INPUT_COLUMNS,
    OUTPUT_COLUMNS,
    Normalizer,
    SpinodoidDataset,
    load_csv,
)
from .metrics import r2_score
from .models import ForwardNet, InverseNet

__all__ = [
    "INPUT_COLUMNS",
    "OUTPUT_COLUMNS",
    "Normalizer",
    "SpinodoidDataset",
    "load_csv",
    "ForwardNet",
    "InverseNet",
    "ForwardConfig",
    "InverseConfig",
    "r2_score",
]

__version__ = "0.1.0"
