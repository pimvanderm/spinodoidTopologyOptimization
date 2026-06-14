"""Network architectures: the forward (f-NN) and inverse (i-NN) surrogates.

Both are fully-connected ReLU MLPs. Linear layers use Glorot/Xavier uniform
weight initialization with zero bias, matching MATLAB's default
``fullyConnectedLayer`` initialization.
"""

from __future__ import annotations

from typing import Sequence

import torch
import torch.nn as nn

from .data import N_INPUTS, N_OUTPUTS


def _mlp(in_dim: int, hidden: Sequence[int], out_dim: int) -> nn.Sequential:
    """Build a ReLU MLP with a linear (un-activated) output layer."""
    layers: list[nn.Module] = []
    prev = in_dim
    for width in hidden:
        layers.append(nn.Linear(prev, width))
        layers.append(nn.ReLU())
        prev = width
    layers.append(nn.Linear(prev, out_dim))
    return nn.Sequential(*layers)


def _xavier_init(module: nn.Module) -> None:
    if isinstance(module, nn.Linear):
        nn.init.xavier_uniform_(module.weight)
        if module.bias is not None:
            nn.init.zeros_(module.bias)


class ForwardNet(nn.Module):
    
    def __init__(
        self,
        in_dim: int = N_INPUTS,
        out_dim: int = N_OUTPUTS,
        hidden: Sequence[int] = (64, 64, 64, 64, 64, 64),
    ) -> None:
        super().__init__()
        self.net = _mlp(in_dim, hidden, out_dim)
        self.register_buffer("input_mean", torch.zeros(in_dim))
        self.register_buffer("input_std", torch.ones(in_dim))
        self.apply(_xavier_init)

    def set_input_standardization(
        self, mean: torch.Tensor, std: torch.Tensor, eps: float = 1e-8
    ) -> None:
        """Set the z-score buffers from precomputed training statistics."""
        mean_t = torch.as_tensor(mean, dtype=self.input_mean.dtype)
        std_t = torch.as_tensor(std, dtype=self.input_std.dtype)
        std_t = torch.where(std_t < eps, torch.ones_like(std_t), std_t)
        self.input_mean.copy_(mean_t.to(self.input_mean.device))
        self.input_std.copy_(std_t.to(self.input_std.device))

    def forward(self, theta_norm: torch.Tensor) -> torch.Tensor:
        x = (theta_norm - self.input_mean) / self.input_std
        return self.net(x)


class InverseNet(nn.Module):
    """Stiffness -> structure surrogate (``9 -> 4``)."""

    def __init__(
        self,
        in_dim: int = N_OUTPUTS,
        out_dim: int = N_INPUTS,
        hidden: Sequence[int] = (100, 100, 100, 100, 100, 100),
    ) -> None:
        super().__init__()
        self.net = _mlp(in_dim, hidden, out_dim)
        self.apply(_xavier_init)

    def forward(self, stiffness_norm: torch.Tensor) -> torch.Tensor:
        return self.net(stiffness_norm)
