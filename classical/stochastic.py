"""
Stochastic Wilson-Cowan dynamics.

Extends the deterministic Wilson-Cowan model with additive Wiener noise:

    dE = (1/tau_E)[-E + S_E(input_E)] dt + sigma_E dW_E
    dI = (1/tau_I)[-I + S_I(input_I)] dt + sigma_I dW_I

with W_E, W_I independent standard Wiener processes.

Numerical integration uses the Euler-Maruyama scheme, the simplest
consistent SDE integrator:

    E_{n+1} = E_n + f_E(E_n, I_n) dt + sigma_E * sqrt(dt) * xi_n,
    xi_n ~ N(0, 1)

Note the sqrt(dt) factor: it is what distinguishes the SDE update from a
naïve Euler step and what gives the discrete trajectory the correct
mean-square scaling of Brownian motion.
"""
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .model import WCParams, sigmoid


@dataclass
class StochasticParams:
    """Noise amplitudes for the stochastic Wilson-Cowan model.

    Attributes
    ----------
    sigma_E, sigma_I : float
        Standard deviation of the additive Wiener noise on each
        population. Typical exploratory values: 0.01 - 0.1.
        Larger values produce more obvious noise-induced transitions.
    """
    sigma_E: float = 0.05
    sigma_I: float = 0.05


def euler_maruyama_step(state: np.ndarray,
                        p: WCParams,
                        noise: StochasticParams,
                        dt: float,
                        rng: np.random.Generator) -> np.ndarray:
    """Single Euler-Maruyama step for the stochastic Wilson-Cowan system.

    Parameters
    ----------
    state : ndarray, shape (2,)
        Current ``[E, I]``.
    p : WCParams
        Wilson-Cowan parameters.
    noise : StochasticParams
        Noise amplitudes.
    dt : float
        Time step.
    rng : np.random.Generator
        Source of randomness.

    Returns
    -------
    ndarray, shape (2,)
        Updated state, clipped to the physical interval [0, 1].
    """
    E, I = state
    input_E = p.w_EE * E - p.w_EI * I + p.P
    input_I = p.w_IE * E - p.w_II * I + p.Q

    drift_E = (-E + sigmoid(input_E, p.a_E, p.theta_E)) / p.tau_E
    drift_I = (-I + sigmoid(input_I, p.a_I, p.theta_I)) / p.tau_I

    sqrt_dt = np.sqrt(dt)
    xi = rng.standard_normal(2)

    E_new = E + drift_E * dt + noise.sigma_E * sqrt_dt * xi[0]
    I_new = I + drift_I * dt + noise.sigma_I * sqrt_dt * xi[1]

    # Reflect at the physical boundaries: activities cannot leave [0, 1]
    return np.clip([E_new, I_new], 0.0, 1.0)


def simulate_stochastic(params: WCParams,
                        noise: StochasticParams,
                        E0: float = 0.1,
                        I0: float = 0.1,
                        t_max: float = 100.0,
                        dt: float = 0.01,
                        seed: Optional[int] = None
                        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Integrate the stochastic Wilson-Cowan system from a given start.

    Parameters
    ----------
    params : WCParams
        Deterministic part of the model.
    noise : StochasticParams
        Noise amplitudes.
    E0, I0 : float
        Initial conditions.
    t_max : float
        Final integration time.
    dt : float
        Euler-Maruyama step size. For SDEs use a smaller step than for
        ODEs (here dt=0.01 is reasonable; deterministic part used 0.1).
    seed : int or None
        If provided, makes a single trajectory reproducible.

    Returns
    -------
    t : ndarray, shape (n,)
    E : ndarray, shape (n,)
    I : ndarray, shape (n,)
    """
    rng = np.random.default_rng(seed)
    n_steps = int(t_max / dt) + 1

    t = np.linspace(0.0, t_max, n_steps)
    E_arr = np.empty(n_steps)
    I_arr = np.empty(n_steps)
    E_arr[0], I_arr[0] = E0, I0

    state = np.array([E0, I0])
    for k in range(1, n_steps):
        state = euler_maruyama_step(state, params, noise, dt, rng)
        E_arr[k], I_arr[k] = state

    return t, E_arr, I_arr


def ensemble_simulate(params: WCParams,
                      noise: StochasticParams,
                      n_trials: int,
                      E0: float = 0.1,
                      I0: float = 0.1,
                      t_max: float = 100.0,
                      dt: float = 0.01,
                      seed: Optional[int] = None
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run an ensemble of independent stochastic trajectories.

    Parameters
    ----------
    n_trials : int
        Number of independent trajectories.
    seed : int or None
        Seed for the master RNG; each trial uses an independent child
        stream, so the ensemble is reproducible as a whole.

    Returns
    -------
    t : ndarray, shape (n_steps,)
    E_ensemble : ndarray, shape (n_trials, n_steps)
    I_ensemble : ndarray, shape (n_trials, n_steps)
    """
    n_steps = int(t_max / dt) + 1
    t = np.linspace(0.0, t_max, n_steps)
    E_all = np.empty((n_trials, n_steps))
    I_all = np.empty((n_trials, n_steps))

    parent_rng = np.random.default_rng(seed)
    child_seeds = parent_rng.integers(0, 2**31 - 1, size=n_trials)

    for k in range(n_trials):
        _, E_all[k], I_all[k] = simulate_stochastic(
            params, noise, E0=E0, I0=I0,
            t_max=t_max, dt=dt, seed=int(child_seeds[k]),
        )

    return t, E_all, I_all


# --------------------------------------------------------------------------- #
# Critical-slowing-down statistics
# --------------------------------------------------------------------------- #

def variance_in_window(trajectory: np.ndarray,
                        window: int = 500) -> np.ndarray:
    """Rolling variance of a 1-D signal in a sliding window.

    The endpoints (where the window does not fully fit) are filled with
    NaN so that downstream code can detect them.

    Parameters
    ----------
    trajectory : ndarray, shape (n,)
        Time series.
    window : int
        Window size in samples.
    """
    n = trajectory.size
    var = np.full(n, np.nan)
    half = window // 2
    for i in range(half, n - half):
        var[i] = trajectory[i - half:i + half].var()
    return var


def lag1_autocorr_in_window(trajectory: np.ndarray,
                             window: int = 500) -> np.ndarray:
    """Rolling lag-1 autocorrelation, a standard slowing-down indicator.

    A signal that returns to its mean quickly has low lag-1
    autocorrelation; a signal that returns slowly (near a bifurcation)
    has lag-1 autocorrelation approaching 1.
    """
    n = trajectory.size
    ac = np.full(n, np.nan)
    half = window // 2
    for i in range(half, n - half):
        chunk = trajectory[i - half:i + half]
        x = chunk[:-1] - chunk[:-1].mean()
        y = chunk[1:] - chunk[1:].mean()
        denom = np.sqrt((x ** 2).sum() * (y ** 2).sum())
        ac[i] = (x * y).sum() / denom if denom > 0 else np.nan
    return ac
