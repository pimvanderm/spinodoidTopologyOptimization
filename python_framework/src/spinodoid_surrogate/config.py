"""Hyperparameter containers mirroring the MATLAB ``trainingOptions`` blocks."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence


@dataclass(frozen=True)
class ForwardConfig:
    """Forward model (f-NN) hyperparameters.

    Defaults match ``optionsF`` / ``layersF`` in ``trainNeuralNetworks.m``.
    """

    hidden_layers: Sequence[int] = (128, 128, 64, 64, 32, 32)
    max_epochs: int = 200
    batch_size: int = 64
    learning_rate: float = 1e-3
    lr_drop_period: int = 30      # LearnRateDropPeriod
    lr_drop_factor: float = 0.5   # LearnRateDropFactor
    zscore_inputs: bool = True     # featureInputLayer(..., "Normalization", "zscore")


@dataclass(frozen=True)
class InverseConfig:
    """Inverse model (i-NN) hyperparameters.

    Defaults match the custom training loop in ``trainNeuralNetworks.m``.
    """

    hidden_layers: Sequence[int] = field(default_factory=lambda: tuple([100] * 6))
    max_epochs: int = 200
    batch_size: int = 64
    learning_rate: float = 1e-4
    lambda_init: float = 0.5       # parameter-loss weight
    lambda_drop_epoch: int = 40    # lambda -> 0 after this many epochs
