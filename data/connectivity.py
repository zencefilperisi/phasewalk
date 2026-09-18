"""
Functional connectivity computation and FCM-to-graph conversion.

Two entry points:
  - compute_fcm(signal, measure="pearson")  : signal -> FCM
  - fcm_to_graph(fcm, layer="positive")     : FCM -> graph-usable matrix

The graph conversion implements ADR 0008:
  - diagonal is always zeroed (definitional, not tunable)
  - layer="positive"  clips negatives to zero  (Layer P, primary)
  - layer="signed"    keeps negatives as-is    (Layer S, robustness)

Currently only "pearson" is implemented for compute_fcm. PLV and other
measures (ADR 0007 sub-decision 2) can be added behind the same
interface later without changing call sites.
"""
from __future__ import annotations

import numpy as np


def compute_fcm(signal: np.ndarray, measure: str = "pearson") -> np.ndarray:
    """Compute a functional connectivity matrix from a multichannel signal.

    Parameters
    ----------
    signal : ndarray, shape (n_channels, n_samples)
        Channel-time array, one row per electrode.
    measure : {"pearson"}
        Connectivity measure to use. Only "pearson" is implemented.
        Additional measures (PLV, coherence) will be added behind this
        same argument as they are needed.

    Returns
    -------
    fcm : ndarray, shape (n_channels, n_channels)
        Symmetric functional connectivity matrix with 1s on the
        diagonal (untouched — pass through fcm_to_graph for graph use).

    Raises
    ------
    ValueError
        If `signal` is not 2-D, or if `measure` is unknown.
    """
    if signal.ndim != 2:
        raise ValueError(
            f"signal must be 2-D (channels, samples); got ndim={signal.ndim}"
        )

    if measure == "pearson":
        # np.corrcoef treats each ROW as a variable, which matches our
        # (n_channels, n_samples) convention.
        return np.corrcoef(signal)

    raise ValueError(
        f"unknown measure {measure!r}; supported: 'pearson'"
    )


def fcm_to_graph(fcm: np.ndarray, layer: str = "positive") -> np.ndarray:
    """Convert an FCM into a matrix usable as a graph adjacency / Hamiltonian.

    Implements ADR 0008:
      1. Always zero the diagonal (definitional).
      2. For layer="positive": clip negative entries to zero.
         For layer="signed":  keep negatives as-is.

    The returned matrix is safe to pass to any function in
    quantum_walks.walks (as the adjacency argument) and to classical
    diffusion routines — provided the caller has selected the layer
    appropriate to the walk model (see ADR 0008: classical diffusion
    is only well-defined on the positive layer).

    Parameters
    ----------
    fcm : ndarray, shape (n, n)
        Symmetric FCM (e.g. from `compute_fcm`).
    layer : {"positive", "signed"}
        Which ADR-0008 layer to produce.

    Returns
    -------
    graph : ndarray, shape (n, n)
        Symmetric matrix with zero diagonal; entries in [0, 1] for
        layer="positive", in [-1, 1] for layer="signed".

    Raises
    ------
    ValueError
        If fcm is not square/2-D, or if `layer` is unknown.
    """
    if fcm.ndim != 2 or fcm.shape[0] != fcm.shape[1]:
        raise ValueError(
            f"fcm must be a square 2-D matrix; got shape={fcm.shape}"
        )

    # Make a copy so we never mutate the caller's array.
    graph = fcm.copy()
    np.fill_diagonal(graph, 0.0)

    if layer == "positive":
        # Clip negatives to zero. np.clip is more explicit than np.maximum(g, 0).
        graph = np.clip(graph, a_min=0.0, a_max=None)
    elif layer == "signed":
        # Keep signs; nothing further to do.
        pass
    else:
        raise ValueError(
            f"unknown layer {layer!r}; supported: 'positive', 'signed'"
        )

    return graph
