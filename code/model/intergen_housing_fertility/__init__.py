"""No-geography intergenerational housing/fertility model package."""

from .parameters import setup_parameters
from .solver import solve_model

__all__ = ["setup_parameters", "solve_model"]
