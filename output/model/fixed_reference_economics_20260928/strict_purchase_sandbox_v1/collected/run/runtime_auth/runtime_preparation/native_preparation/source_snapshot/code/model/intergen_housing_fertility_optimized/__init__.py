"""Intergenerational housing/fertility model package."""

from .parameters import setup_parameters
from .solver import InfeasibleThetaError, run_model_cp_dt
from .target_system import TargetSystem

__all__ = ["InfeasibleThetaError", "TargetSystem", "setup_parameters", "run_model_cp_dt"]
