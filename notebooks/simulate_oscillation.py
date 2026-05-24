"""
Simulate the Wilson-Cowan model and visualize oscillatory dynamics.

This script demonstrates two foundational behaviours of the Wilson-Cowan
neural mass model:

  1. A single trajectory in the oscillatory regime (P = 1.2), shown both
     as a time series and as a closed orbit in the (E, I) phase plane.

  2. The transition from quasi-stationary behaviour to sustained
     oscillation as the external drive P increases. Each panel uses
     identical initial conditions; only P differs.

All figures are produced first and shown together at the end.

Run as:
    python notebooks/simulate_oscillation.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np

from classical.model import WCParams, simulate
from _style import COLORS, apply_style

apply_style()
np.random.seed(42)


# --------------------------------------------------------------------------- #
# Figure 1 — single trajectory in the oscillatory regime
# --------------------------------------------------------------------------- #
params = WCParams(P=1.2)
t, E, I = simulate(params, E0=0.1, I0=0.05, t_max=50)

fig1, axes = plt.subplots(1, 2, figsize=(13, 4.8),
                          gridspec_kw={"width_ratios": [1.4, 1.0]})

# Left: time series
axes[0].plot(t, E, color=COLORS["excitatory"], lw=1.8, label="E (excitatory)")
axes[0].plot(t, I, color=COLORS["inhibitory"], lw=1.8, label="I (inhibitory)")
axes[0].set_xlabel("time")
axes[0].set_ylabel("population activity")
axes[0].set_ylim(-0.05, 1.05)
axes[0].set_title("(a) Time series", loc="left")
axes[0].legend(loc="upper right", ncol=2)

# Right: phase trajectory
axes[1].plot(E, I, color=COLORS["muted"], lw=1.4, alpha=0.95)
axes[1].scatter([E[0]], [I[0]], s=70, color=COLORS["accent"],
                edgecolor="white", linewidth=1.4, zorder=5, label="t = 0")
axes[1].scatter([E[-1]], [I[-1]], s=70, color=COLORS["unstable"],
                edgecolor="white", linewidth=1.4, zorder=5,
                label=f"t = {t[-1]:.0f}")
axes[1].set_xlabel("E")
axes[1].set_ylabel("I")
axes[1].set_xlim(0, 1)
axes[1].set_ylim(0, 1)
axes[1].set_aspect("equal")
axes[1].set_title("(b) Phase-space trajectory", loc="left")
axes[1].legend(loc="lower right")

fig1.suptitle(f"Sustained oscillation in the Wilson-Cowan model "
              f"at P = {params.P}",
              fontsize=13.5, y=1.02)
fig1.tight_layout()
fig1.savefig(ROOT / "notebooks" / "fig_oscillation_single.png")


# --------------------------------------------------------------------------- #
# Figure 2 — sweep across the bifurcation
# --------------------------------------------------------------------------- #
P_values = [1.0, 1.2, 1.4, 1.6]
fig2, axes = plt.subplots(len(P_values), 1, figsize=(11, 8.5), sharex=True)

for ax, P in zip(axes, P_values):
    params = WCParams(P=P)
    t, E, I = simulate(params, E0=0.1, I0=0.05, t_max=80)
    ax.plot(t, E, color=COLORS["excitatory"], lw=1.5, label="E")
    ax.plot(t, I, color=COLORS["inhibitory"], lw=1.5, label="I")
    ax.set_ylabel("activity")
    ax.set_ylim(-0.05, 1.05)
    ax.text(0.012, 0.86, f"$P = {P}$",
            transform=ax.transAxes, fontsize=11,
            color=COLORS["muted"],
            bbox=dict(facecolor=COLORS["background"],
                      edgecolor=COLORS["grid"],
                      pad=4, boxstyle="round,pad=0.35"))

axes[0].legend(loc="upper right", ncol=2)
axes[-1].set_xlabel("time")

fig2.suptitle("Onset of sustained oscillation as the external drive grows",
              fontsize=13.5, y=0.995)
fig2.tight_layout()
fig2.savefig(ROOT / "notebooks" / "fig_oscillation_sweep.png")


print("Done. Two figures will open together. Close them to end.")
plt.show()
