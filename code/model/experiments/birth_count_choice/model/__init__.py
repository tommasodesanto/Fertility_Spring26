"""Canonical stationary API. Imports perform no setup, engine load, or solve."""
__all__ = ['DEFAULT_INPUTS','DEFAULT_PARAMETERS','DEFAULT_PRICE','load_inputs','bind_parameters','solve_at_price','solve_stationary_ge']
def __getattr__(name):
    if name in __all__[:5]:
        from . import inputs
        return getattr(inputs,name)
    if name in __all__[5:]:
        from . import equilibrium
        return getattr(equilibrium,name)
    raise AttributeError(name)
