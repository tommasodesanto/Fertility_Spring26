"""Intergenerational housing/fertility model package."""

from .parameters import setup_parameters
from .solver import InfeasibleThetaError, run_model_cp_dt

__all__ = ["InfeasibleThetaError", "setup_parameters", "run_model_cp_dt"]
