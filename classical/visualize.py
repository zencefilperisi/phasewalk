"""
Phase-space visualization tools for the Wilson-Cowan model.

Functions
---------
phase_portrait
    Draws the vector field together with the two nullclines.
plot_timeseries
    Plots the E and I activities over time.
plot_trajectory_on_phase
    Overlays one or more integrated trajectories on a phase portrait.
"""
import matplotlib.pyplot as plt
import numpy as np

from .model import WCParams, wilson_cowan_rhs


def _compute_field(params: WCParams, n: int):
    """Evaluate (dE/dt, dI/dt) on an n x n grid over the unit square."""
    E_axis = np.linspace(0, 1, n)
    I_axis = np.linspace(0, 1, n)
    EE, II = np.meshgrid(E_axis, I_axis)

    dE = np.zeros_like(EE)
    dI = np.zeros_like(II)

    for i in range(n):
        for j in range(n):
            d = wilson_cowan_rhs(0.0, [EE[i, j], II[i, j]], params)
            dE[i, j] = d[0]
            dI[i, j] = d[1]

    return EE, II, dE, dI


def phase_portrait(params: WCParams,
                   ax=None,
                   show_nullclines: bool = True,
                   show_vectorfield: bool = True,
                   n_dense: int = 300):
    """Draw the full Wilson-Cowan phase portrait.

    The figure layers:
      - vector field as streamlines, coloured by flow magnitude;
      - nullclines as the zero-level contours of dE/dt (red)
        and dI/dt (blue). Their intersections are the fixed points.

    Parameters
    ----------
    params : WCParams
        Model parameters.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. A new figure is created if ``None``.
    show_nullclines, show_vectorfield : bool
        Toggle each layer.
    n_dense : int
        Grid resolution for both the streamplot and contour computation.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    EE, II, dE, dI = _compute_field(params, n_dense)

    if show_vectorfield:
        magnitude = np.sqrt(dE ** 2 + dI ** 2) + 1e-12
        ax.streamplot(
            EE, II, dE, dI,
            color=magnitude, cmap="Greys",
            density=1.4, linewidth=0.8, arrowsize=1.0,
        )

    if show_nullclines:
        ax.contour(EE, II, dE, levels=[0],
                   colors="crimson", linewidths=2.0)
        ax.contour(EE, II, dI, levels=[0],
                   colors="royalblue", linewidths=2.0)

    ax.set_xlabel("E (excitatory activity)")
    ax.set_ylabel("I (inhibitory activity)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")

    # Proxy artists for the legend
    if show_nullclines:
        ax.plot([], [], color="crimson", lw=2,
                label="dE/dt = 0  (E-nullcline)")
        ax.plot([], [], color="royalblue", lw=2,
                label="dI/dt = 0  (I-nullcline)")
        ax.legend(loc="upper right", framealpha=0.9)

    return ax


def plot_timeseries(t, E, I, ax=None, title: str = ""):
    """Plot E and I population activities versus time.

    Parameters
    ----------
    t, E, I : array_like
        Time grid and the two activity series.
    ax : matplotlib.axes.Axes, optional
    title : str
        Optional axes title.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 4))

    ax.plot(t, E, color="crimson", lw=2, label="E (excitatory)")
    ax.plot(t, I, color="royalblue", lw=2, label="I (inhibitory)")
    ax.set_xlabel("time")
    ax.set_ylabel("activity")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")
    if title:
        ax.set_title(title)
    return ax


def plot_trajectory_on_phase(params: WCParams, trajectories, ax=None):
    """Overlay one or more integrated trajectories on a phase portrait.

    Parameters
    ----------
    params : WCParams
    trajectories : iterable of (E_array, I_array) tuples
        Each tuple is one trajectory's (E, I) time series.
    ax : matplotlib.axes.Axes, optional
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    phase_portrait(params, ax=ax)

    for E_arr, I_arr in trajectories:
        ax.plot(E_arr, I_arr, color="black", lw=1.3, alpha=0.75)
        ax.scatter(E_arr[0], I_arr[0], color="green",
                   s=40, zorder=5, edgecolor="black")
        ax.scatter(E_arr[-1], I_arr[-1], color="red",
                   s=40, zorder=5, edgecolor="black")

    return ax
