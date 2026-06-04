"""Intergenerational housing/fertility model package."""

from .parameters import setup_parameters
from .solver import run_model_cp_dt

__all__ = ["setup_parameters", "run_model_cp_dt"]
