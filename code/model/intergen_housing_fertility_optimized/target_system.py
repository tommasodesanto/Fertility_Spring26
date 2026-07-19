"""Immutable target-system contract used by calibration and reporting."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class TargetSystem:
    """Ordered SMM target values and weights with a stable fingerprint."""

    name: str
    moment_names: tuple[str, ...]
    target_values: tuple[float, ...]
    weights: tuple[float, ...]

    @classmethod
    def from_mappings(
        cls,
        name: str,
        targets: Mapping[str, float],
        weights: Mapping[str, float],
    ) -> "TargetSystem":
        names = tuple(str(moment) for moment in targets)
        if set(names) != set(weights):
            raise ValueError("target and weight names must match exactly")
        return cls(
            name=str(name),
            moment_names=names,
            target_values=tuple(float(targets[moment]) for moment in names),
            weights=tuple(float(weights[moment]) for moment in names),
        )

    @property
    def count(self) -> int:
        return len(self.moment_names)

    @property
    def fingerprint(self) -> str:
        payload = [
            [name, target, weight]
            for name, target, weight in zip(
                self.moment_names,
                self.target_values,
                self.weights,
            )
        ]
        encoded = json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode()
        return hashlib.sha256(encoded).hexdigest()

    def targets_dict(self) -> dict[str, float]:
        return dict(zip(self.moment_names, self.target_values))

    def weights_dict(self) -> dict[str, float]:
        return dict(zip(self.moment_names, self.weights))

    def loss(self, moments: Mapping[str, float]) -> float:
        """Return the weighted squared-error objective in declared row order."""

        missing = [name for name in self.moment_names if name not in moments]
        if missing:
            raise KeyError(f"missing target moments: {missing}")
        return sum(
            weight * (float(moments[name]) - target) ** 2
            for name, target, weight in zip(
                self.moment_names,
                self.target_values,
                self.weights,
            )
        )

    def require_identified(self, free_parameter_count: int) -> None:
        if self.count < int(free_parameter_count):
            raise ValueError(
                f"target system {self.name!r} is underidentified: "
                f"{self.count} moments for {int(free_parameter_count)} free parameters"
            )
