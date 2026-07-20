"""Named state-axis and family-state layout contracts."""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace


@dataclass(frozen=True)
class StateLayout:
    wealth: int
    tenure: int
    location: int
    age: int
    income: int
    parity: int
    child_state: int

    @classmethod
    def from_parameters(cls, P: SimpleNamespace, wealth_points: int) -> "StateLayout":
        return cls(
            wealth=int(wealth_points),
            tenure=1 + int(P.n_house),
            location=int(P.I),
            age=int(P.J),
            income=int(getattr(P, "Nz", 1)),
            parity=int(P.n_parity),
            child_state=int(P.n_child_states),
        )

    @property
    def markov_shape(self) -> tuple[int, ...]:
        return (
            self.wealth,
            self.tenure,
            self.location,
            self.age,
            self.income,
            self.parity,
            self.child_state,
        )

    @property
    def nonmarkov_shape(self) -> tuple[int, ...]:
        return (
            self.wealth,
            self.tenure,
            self.location,
            self.age,
            self.parity,
            self.child_state,
        )

    @property
    def family_state_count(self) -> int:
        return self.parity * self.child_state

    def encode_family_state(self, parity: int, child_state: int) -> int:
        if not 0 <= int(parity) < self.parity:
            raise IndexError("parity index is outside the state layout")
        if not 0 <= int(child_state) < self.child_state:
            raise IndexError("child-state index is outside the state layout")
        return int(parity) + self.parity * int(child_state)

    def decode_family_state(self, column: int) -> tuple[int, int]:
        if not 0 <= int(column) < self.family_state_count:
            raise IndexError("flattened family-state index is outside the state layout")
        return int(column) % self.parity, int(column) // self.parity

    def require_markov_shape(self, name: str, shape: tuple[int, ...]) -> None:
        if tuple(shape) != self.markov_shape:
            raise ValueError(f"{name} has shape {tuple(shape)}, expected {self.markov_shape}")
