"""
Fixed-point analysis of the Wilson-Cowan model.

Provides
--------
sigmoid_derivative
    Closed-form derivative S'(x) = a * S(x) * (1 - S(x)).
jacobian
    Analytic 2x2 Jacobian of the Wilson-Cowan system.
classify_fixed_point
    Categorize a fixed point from its Jacobian eigenvalues.
find_fixed_points
    Locate all fixed points by multi-start fsolve.
find_hopf_threshold
    Bisection search for the parameter value where the dominant
    eigenvalue's real part crosses zero.
"""
from dataclasses import dataclass
from typing import Callable, List, Optional, Tuple

import numpy as np
from scipy.optimize import fsolve

from .model import WCParams, sigmoid, wilson_cowan_rhs


# --------------------------------------------------------------------------- #
# Analytic derivatives
# --------------------------------------------------------------------------- #

def sigmoid_derivative(x, a, theta):
    """Closed-form derivative of the logistic sigmoid.

        S'(x) = a * S(x) * (1 - S(x))

    Numerically more stable and faster than finite differences.
    """
    s = sigmoid(x, a, theta)
    return a * s * (1.0 - s)


def jacobian(state, p: WCParams) -> np.ndarray:
    """Analytic Jacobian of the Wilson-Cowan system at the given state.

    Parameters
    ----------
    state : array_like, shape (2,)
        ``[E, I]`` at which to evaluate.
    p : WCParams

    Returns
    -------
    ndarray, shape (2, 2)
        Jacobian matrix. Rows correspond to (dE/dt, dI/dt); columns to
        partial derivatives with respect to (E, I).
    """
    E, I = state

    input_E = p.w_EE * E - p.w_EI * I + p.P
    input_I = p.w_IE * E - p.w_II * I + p.Q

    SpE = sigmoid_derivative(input_E, p.a_E, p.theta_E)
    SpI = sigmoid_derivative(input_I, p.a_I, p.theta_I)

    J = np.array([
        [(-1.0 + p.w_EE * SpE) / p.tau_E,  (-p.w_EI * SpE) / p.tau_E],
        [( p.w_IE * SpI) / p.tau_I,        (-1.0 - p.w_II * SpI) / p.tau_I],
    ])
    return J


# --------------------------------------------------------------------------- #
# Stability classification
# --------------------------------------------------------------------------- #

@dataclass
class FixedPoint:
    """A fixed point together with its stability classification."""
    E: float
    I: float
    eigenvalues: np.ndarray   # complex array, shape (2,)
    classification: str       # e.g. 'stable focus', 'saddle', ...

    @property
    def is_stable(self) -> bool:
        return np.all(self.eigenvalues.real < 0)

    @property
    def max_real_part(self) -> float:
        """Largest real part among the eigenvalues — useful for tracking
        Hopf bifurcations."""
        return float(np.max(self.eigenvalues.real))


def classify_fixed_point(eigenvalues: np.ndarray,
                          tol: float = 1e-8) -> str:
    """Label a fixed point based on the two Jacobian eigenvalues."""
    re = eigenvalues.real
    im = eigenvalues.imag

    is_complex = np.any(np.abs(im) > tol)

    if is_complex:
        if np.all(re < -tol):
            return "stable focus"
        if np.all(re > tol):
            return "unstable focus"
        if np.all(np.abs(re) < tol):
            return "center"
        return "complex (mixed)"
    else:
        if np.all(re < -tol):
            return "stable node"
        if np.all(re > tol):
            return "unstable node"
        if re[0] * re[1] < 0:
            return "saddle"
        return "degenerate"


# --------------------------------------------------------------------------- #
# Fixed-point search
# --------------------------------------------------------------------------- #

def _rhs_for_solver(state, p: WCParams):
    """Wrapper returning the RHS as a numpy array (fsolve expects this)."""
    return np.array(wilson_cowan_rhs(0.0, state, p))


def find_fixed_points(p: WCParams,
                       n_grid: int = 8,
                       dedup_tol: float = 1e-4) -> List[FixedPoint]:
    """Find all fixed points of the Wilson-Cowan system.

    Strategy: run fsolve from many initial guesses on a uniform grid
    over the unit square, then deduplicate solutions that are within
    ``dedup_tol`` of each other.

    Parameters
    ----------
    p : WCParams
        Model parameters.
    n_grid : int
        Per-axis resolution of the multi-start grid (n_grid^2 attempts).
    dedup_tol : float
        Distance below which two solutions are considered identical.

    Returns
    -------
    list of FixedPoint
        All distinct fixed points found, classified by stability.
    """
    candidates: List[Tuple[float, float]] = []

    grid = np.linspace(0.02, 0.98, n_grid)
    for E0 in grid:
        for I0 in grid:
            try:
                sol, info, ier, _ = fsolve(
                    _rhs_for_solver, x0=[E0, I0], args=(p,),
                    full_output=True, xtol=1e-12,
                )
            except Exception:
                continue

            if ier != 1:
                continue

            E_s, I_s = sol
            # Accept only solutions inside the physical [0, 1] box
            if not (0.0 <= E_s <= 1.0 and 0.0 <= I_s <= 1.0):
                continue

            # Verify the residual is small
            if np.linalg.norm(_rhs_for_solver(sol, p)) > 1e-6:
                continue

            candidates.append((E_s, I_s))

    # Deduplicate
    unique: List[Tuple[float, float]] = []
    for cand in candidates:
        if not any(np.hypot(cand[0] - u[0], cand[1] - u[1]) < dedup_tol
                   for u in unique):
            unique.append(cand)

    # Classify each unique fixed point
    fixed_points = []
    for E_s, I_s in unique:
        J = jacobian([E_s, I_s], p)
        evals = np.linalg.eigvals(J)
        fp = FixedPoint(
            E=E_s, I=I_s,
            eigenvalues=evals,
            classification=classify_fixed_point(evals),
        )
        fixed_points.append(fp)

    # Sort by E for reproducible ordering
    fixed_points.sort(key=lambda fp: fp.E)
    return fixed_points


# --------------------------------------------------------------------------- #
# Hopf-bifurcation locator
# --------------------------------------------------------------------------- #

def _max_real_eigenvalue(p: WCParams) -> Optional[float]:
    """Largest real part of any Jacobian eigenvalue, across all fixed
    points. Returns None if no fixed point is found."""
    fps = find_fixed_points(p)
    if not fps:
        return None
    return max(fp.max_real_part for fp in fps)


def find_hopf_threshold(base_params: WCParams,
                         param_name: str = "P",
                         lo: float = 0.0,
                         hi: float = 3.0,
                         tol: float = 1e-5,
                         max_iter: int = 60) -> Optional[float]:
    """Locate a Hopf bifurcation by bisection on a single parameter.

    Searches for the value of ``param_name`` at which the dominant
    eigenvalue's real part crosses zero, while keeping all other
    parameters fixed at ``base_params``.

    Parameters
    ----------
    base_params : WCParams
    param_name : str
        Attribute name on WCParams to vary (e.g., 'P', 'w_EI').
    lo, hi : float
        Bracket. The function assumes the sign of the dominant real part
        changes between ``lo`` and ``hi``; otherwise returns None.
    tol : float
        Stop when ``hi - lo < tol``.
    max_iter : int
        Hard cap on bisection iterations.

    Returns
    -------
    float or None
        Estimated Hopf threshold; None if no sign change is detected.
    """
    def evaluate(value: float) -> Optional[float]:
        params = WCParams(**{**base_params.__dict__, param_name: value})
        return _max_real_eigenvalue(params)

    f_lo = evaluate(lo)
    f_hi = evaluate(hi)
    if f_lo is None or f_hi is None:
        return None
    if f_lo * f_hi > 0:
        # No sign change in the bracket
        return None

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = evaluate(mid)
        if f_mid is None:
            return None
        if f_mid * f_lo < 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
        if hi - lo < tol:
            break

    return 0.5 * (lo + hi)
