"""
Functional connectivity computation and FCM-to-graph conversion.

Two entry points:
  - compute_fcm(signal, measure="pearson"|"plv")  : signal -> FCM
  - fcm_to_graph(fcm, layer="positive")           : FCM -> graph-usable matrix

compute_fcm measures:
  - "pearson" : Pearson correlation (ADR 0007 primary). Values in [-1, 1].
  - "plv"     : phase-locking value (ADR 0007 robustness). Broadband
                (0.5-150 Hz, as filtered in the dataset), a coarse
                phase-synchronisation measure. Values in [0, 1]. Used
                only as a robustness check against the Pearson findings,
                never as the primary measure. Band-resolved PLV is a
                Phase 3 refinement (ADR 0007 note 2026-10-07).

The graph conversion implements ADR 0008:
  - diagonal is always zeroed (definitional, not tunable)
  - layer="positive"  clips negatives to zero  (Layer P, primary)
  - layer="signed"    keeps negatives as-is    (Layer S, robustness)

Note: PLV is already non-negative (in [0, 1]), so for a PLV FCM the
"positive" and "signed" layers differ only trivially (there are no
negatives to clip); the diagonal is still zeroed.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import hilbert


def _plv_matrix(signal: np.ndarray) -> np.ndarray:
    """Phase-locking value matrix from a multichannel signal.

    For each channel the instantaneous phase is extracted via the
    analytic (Hilbert) signal. For channels i, j:

        PLV_ij = | (1/T) sum_t exp( i (phi_i(t) - phi_j(t)) ) |

    which is 1 for a constant phase difference (perfect locking) and
    tends to 0 for a drifting/random phase difference.
    """
    analytic = hilbert(signal, axis=1)          # (n_channels, n_samples), complex
    phase = np.angle(analytic)                  # instantaneous phase
    comp = np.exp(1j * phase)                    # unit phasors
    n_samples = signal.shape[1]
    # M_ij = sum_t comp_i(t) * conj(comp_j(t)) ; |M|/T is the PLV.
    m = comp @ comp.conj().T
    plv = np.abs(m) / n_samples
    return plv


def compute_fcm(signal: np.ndarray, measure: str = "pearson") -> np.ndarray:
    """Compute a functional connectivity matrix from a multichannel signal.

    Parameters
    ----------
    signal : ndarray, shape (n_channels, n_samples)
        Channel-time array, one row per electrode.
    measure : {"pearson", "plv"}
        Connectivity measure. "pearson" (primary, ADR 0007) returns
        values in [-1, 1]; "plv" (robustness) returns values in [0, 1].

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
        # np.corrcoef treats each ROW as a variable, matching our
        # (n_channels, n_samples) convention.
        return np.corrcoef(signal)

    if measure == "plv":
        return _plv_matrix(signal)

    raise ValueError(
        f"unknown measure {measure!r}; supported: 'pearson', 'plv'"
    )


def fcm_to_graph(fcm: np.ndarray, layer: str = "positive") -> np.ndarray:
    """Convert an FCM into a matrix usable as a graph adjacency / Hamiltonian.

    Implements ADR 0008:
      1. Always zero the diagonal (definitional).
      2. For layer="positive": clip negative entries to zero.
         For layer="signed":  keep negatives as-is.

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

    graph = fcm.copy()
    np.fill_diagonal(graph, 0.0)

    if layer == "positive":
        graph = np.clip(graph, a_min=0.0, a_max=None)
    elif layer == "signed":
        pass
    else:
        raise ValueError(
            f"unknown layer {layer!r}; supported: 'positive', 'signed'"
        )

    return graph
