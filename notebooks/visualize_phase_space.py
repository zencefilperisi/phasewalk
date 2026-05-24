"""
Visualize the Wilson-Cowan phase plane.

This script draws the full phase-space picture for several values of the
external drive P:

  - Streamlines show the local direction of (dE/dt, dI/dt). Their
    intensity (cividis colormap) encodes the magnitude of the flow.
  - The two nullclines (dE/dt = 0 and dI/dt = 0) intersect at the
    fixed points of the system.
  - Sample trajectories from random initial conditions illustrate
    whether the system relaxes to a fixed point or to a closed orbit.

All figures are produced first and shown together at the end.

Run as:
    python notebooks/visualize_phase_space.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from classical.model import WCParams, simulate, wilson_cowan_rhs
from _style import COLORS, SEQ_CMAP, apply_style

apply_style()
np.random.seed(42)


def _evaluate_field(params, n=240):
    """Evaluate (dE/dt, dI/dt) on an n × n grid over the unit square."""
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


def draw_phase_plane(params, ax, n=240):
    EE, II, dE, dI = _evaluate_field(params, n)
    speed = np.sqrt(dE ** 2 + dI ** 2) + 1e-12

    # Streamlines coloured by flow magnitude (cividis)
    ax.streamplot(EE, II, dE, dI,
                  color=speed, cmap=SEQ_CMAP,
                  density=1.25, linewidth=0.7, arrowsize=0.85)

    # Nullclines: thicker, semantic colours
    ax.contour(EE, II, dE, levels=[0],
               colors=COLORS["excitatory"], linewidths=2.0)
    ax.contour(EE, II, dI, levels=[0],
               colors=COLORS["inhibitory"], linewidths=2.0)

    ax.set_xlabel("E (excitatory activity)")
    ax.set_ylabel("I (inhibitory activity)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.grid(False)


def overlay_trajectories(params, ax, n_traj=8, t_max=30):
    for _ in range(n_traj):
        E0, I0 = np.random.uniform(0.05, 0.95, size=2)
        _, E, I = simulate(params, E0=E0, I0=I0, t_max=t_max)
        ax.plot(E, I, color=COLORS["muted"], lw=1.05, alpha=0.85)
        ax.scatter([E[0]], [I[0]], s=24, color=COLORS["accent"],
                   edgecolor="white", linewidth=0.7, zorder=6)
        ax.scatter([E[-1]], [I[-1]], s=24, color=COLORS["unstable"],
                   edgecolor="white", linewidth=0.7, zorder=6)


# Shared legend handles for both figures
LEGEND_HANDLES = [
    Line2D([0], [0], color=COLORS["excitatory"], lw=2.0,
           label=r"$dE/dt = 0$  (E-nullcline)"),
    Line2D([0], [0], color=COLORS["inhibitory"], lw=2.0,
           label=r"$dI/dt = 0$  (I-nullcline)"),
    Line2D([0], [0], marker="o", color="w",
           markerfacecolor=COLORS["accent"], markeredgecolor="white",
           markersize=7, label="trajectory start"),
    Line2D([0], [0], marker="o", color="w",
           markerfacecolor=COLORS["unstable"], markeredgecolor="white",
           markersize=7, label="trajectory end"),
]


# --------------------------------------------------------------------------- #
# Figure 1 — single phase plane at P = 1.2
# --------------------------------------------------------------------------- #
params = WCParams(P=1.2)
fig1, ax = plt.subplots(figsize=(7.6, 7.6))
draw_phase_plane(params, ax)
overlay_trajectories(params, ax, n_traj=8)
ax.legend(handles=LEGEND_HANDLES, loc="upper right")

fig1.suptitle(f"Phase plane of the Wilson-Cowan model at P = {params.P}",
              fontsize=13.5, y=0.995)
fig1.tight_layout()
fig1.savefig(ROOT / "notebooks" / "fig_phase_plane_single.png")


# --------------------------------------------------------------------------- #
# Figure 2 — comparison across three values of P
# --------------------------------------------------------------------------- #
P_values = [1.0, 1.3, 1.6]
fig2, axes = plt.subplots(1, 3, figsize=(17, 6.4))

for ax, P in zip(axes, P_values):
    params = WCParams(P=P)
    draw_phase_plane(params, ax)
    overlay_trajectories(params, ax, n_traj=5, t_max=40)
    ax.set_title(f"$P = {P}$", loc="left")

# Shared legend across the top
fig2.legend(handles=LEGEND_HANDLES, loc="upper center", ncol=4,
            bbox_to_anchor=(0.5, 1.04), frameon=False)
fig2.suptitle("Deformation of the phase plane with the external drive",
              fontsize=13.5, y=1.10)
fig2.tight_layout()
fig2.savefig(ROOT / "notebooks" / "fig_phase_plane_sweep.png")


print("Done. Two figures will open together. Close them to end.")
plt.show()
