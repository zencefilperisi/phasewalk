"""
Core implementation of the Wilson-Cowan model.

Mean-field equations for two coupled neural populations,
excitatory (E) and inhibitory (I):

    tau_E * dE/dt = -E + S_E(w_EE * E - w_EI * I + P)
    tau_I * dI/dt = -I + S_I(w_IE * E - w_II * I + Q)

Reference
---------
Wilson, H. R., & Cowan, J. D. (1972). Excitatory and inhibitory
interactions in localized populations of model neurons.
Biophysical Journal, 12(1), 1-24.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass
class WCParams:
    """Container for Wilson-Cowan parameters.

    Defaults are chosen near the typical ranges from Wilson & Cowan (1972).
    Sweeping ``P`` is the standard way to observe the transition from a
    healthy fixed-point regime to a limit-cycle (seizure-like) regime.

    Attributes
    ----------
    w_EE, w_EI, w_IE, w_II : float
        Synaptic coupling strengths between the two populations.
    P, Q : float
        External inputs to the E and I populations, respectively.
    tau_E, tau_I : float
        Population time constants (normalized; think milliseconds).
    a_E, theta_E, a_I, theta_I : float
        Slope and threshold of the sigmoid activation for each population.
    """

    # Coupling strengths
    w_EE: float = 16.0   # E -> E (recurrent excitation)
    w_EI: float = 12.0   # I -> E (inhibition onto E)
    w_IE: float = 15.0   # E -> I
    w_II: float = 3.0    # I -> I

    # External drives
    P: float = 1.2       # drive onto excitatory population
    Q: float = 0.0       # drive onto inhibitory population

    # Time constants
    tau_E: float = 1.0
    tau_I: float = 2.0   # inhibitory population is typically slower

    # Sigmoid parameters (separate for E and I)
    a_E: float = 1.3
    theta_E: float = 4.0
    a_I: float = 2.0
    theta_I: float = 3.7


def sigmoid(x, a, theta):
    """Logistic activation function used in the Wilson-Cowan model.

        S(x) = 1 / (1 + exp(-a * (x - theta)))

    Parameters
    ----------
    x : float or array_like
        Input value (total synaptic drive to a population).
    a : float
        Slope. Larger values produce a sharper threshold-like transition.
    theta : float
        Half-activation point: S(theta) = 0.5.

    Returns
    -------
    float or ndarray
        Activation in [0, 1].
    """
    # Clip the exponent to avoid overflow at extreme inputs that the
    # root-finder might probe. The sigmoid saturates anyway beyond ~|z|=500.
    z = np.clip(-a * (x - theta), -500.0, 500.0)
    return 1.0 / (1.0 + np.exp(z))


def wilson_cowan_rhs(t, state, p: WCParams):
    """Right-hand side of the Wilson-Cowan ODE system.

    Compatible with ``scipy.integrate.solve_ivp`` signature.

    Parameters
    ----------
    t : float
        Time (unused for this autonomous system, but required by the
        solver interface).
    state : array_like, shape (2,)
        Current state ``[E, I]`` — instantaneous activities of the
        excitatory and inhibitory populations.
    p : WCParams
        Model parameters.

    Returns
    -------
    list of float
        ``[dE/dt, dI/dt]``.
    """
    E, I = state

    input_E = p.w_EE * E - p.w_EI * I + p.P
    input_I = p.w_IE * E - p.w_II * I + p.Q

    dE_dt = (-E + sigmoid(input_E, p.a_E, p.theta_E)) / p.tau_E
    dI_dt = (-I + sigmoid(input_I, p.a_I, p.theta_I)) / p.tau_I

    return [dE_dt, dI_dt]


def simulate(params: WCParams,
             E0: float = 0.1,
             I0: float = 0.1,
             t_max: float = 100.0,
             dt: float = 0.1,
             method: str = "RK45"):
    """Integrate the Wilson-Cowan system from a given initial condition.

    Parameters
    ----------
    params : WCParams
        Model parameters.
    E0, I0 : float
        Initial conditions, expected in [0, 1].
    t_max : float
        Final integration time.
    dt : float
        Output sampling step (the solver itself uses adaptive steps).
    method : str
        Integration method passed to ``solve_ivp``. ``"RK45"`` is a
        sensible default; switch to ``"LSODA"`` for stiff regimes.

    Returns
    -------
    t : ndarray, shape (n,)
        Time grid.
    E : ndarray, shape (n,)
        Excitatory population activity.
    I : ndarray, shape (n,)
        Inhibitory population activity.

    Raises
    ------
    RuntimeError
        If the solver fails to converge.
    """
    t_eval = np.arange(0, t_max, dt)

    sol = solve_ivp(
        fun=lambda t, y: wilson_cowan_rhs(t, y, params),
        t_span=(0, t_max),
        y0=[E0, I0],
        t_eval=t_eval,
        method=method,
        rtol=1e-8,
        atol=1e-10,
    )

    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")

    return sol.t, sol.y[0], sol.y[1]
