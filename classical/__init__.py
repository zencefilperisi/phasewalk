"""
Phasewalk — classical module.

Classical and stochastic analysis of the Wilson-Cowan model.
"""
from .analysis import (
    FixedPoint,
    classify_fixed_point,
    find_fixed_points,
    find_hopf_threshold,
    jacobian,
    sigmoid_derivative,
)
from .model import WCParams, sigmoid, simulate, wilson_cowan_rhs
from .stochastic import (
    StochasticParams,
    ensemble_simulate,
    euler_maruyama_step,
    lag1_autocorr_in_window,
    simulate_stochastic,
    variance_in_window,
)
from .visualize import phase_portrait, plot_timeseries, plot_trajectory_on_phase

__all__ = [
    # model
    "WCParams",
    "sigmoid",
    "wilson_cowan_rhs",
    "simulate",
    # visualize
    "phase_portrait",
    "plot_timeseries",
    "plot_trajectory_on_phase",
    # analysis
    "FixedPoint",
    "sigmoid_derivative",
    "jacobian",
    "classify_fixed_point",
    "find_fixed_points",
    "find_hopf_threshold",
    # stochastic
    "StochasticParams",
    "euler_maruyama_step",
    "simulate_stochastic",
    "ensemble_simulate",
    "variance_in_window",
    "lag1_autocorr_in_window",
]
