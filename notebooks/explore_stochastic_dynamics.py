"""
Explore stochastic Wilson-Cowan dynamics.

Three demonstrations:

  (1) A single deterministic trajectory and a single stochastic trajectory
      at the same operating point, side by side.

  (2) An ensemble of stochastic trajectories from identical initial
      conditions — the spread illustrates the loss of determinism.

  (3) Phase-plane scatter of the ensemble end-states, showing where the
      stochastic dynamics concentrates probability mass.

Run as:
    python notebooks/explore_stochastic_dynamics.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np

from classical.model import WCParams, simulate
from classical.stochastic import (
    StochasticParams, ensemble_simulate, simulate_stochastic,
)
from _style import COLORS, apply_style

apply_style()


# --------------------------------------------------------------------------- #
# Common setup
# --------------------------------------------------------------------------- #
params = WCParams(P=0.6)              # below Hopf threshold — stable regime
noise = StochasticParams(sigma_E=0.04, sigma_I=0.04)
E0, I0 = 0.1, 0.05
t_max = 60.0


# --------------------------------------------------------------------------- #
# Figure 1 — deterministic vs single stochastic run
# --------------------------------------------------------------------------- #
t_det, E_det, I_det = simulate(params, E0=E0, I0=I0, t_max=t_max)
t_sto, E_sto, I_sto = simulate_stochastic(
    params, noise, E0=E0, I0=I0, t_max=t_max, dt=0.01, seed=0)

fig1, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

axes[0].plot(t_det, E_det, color=COLORS["excitatory"], lw=1.8, label="E")
axes[0].plot(t_det, I_det, color=COLORS["inhibitory"], lw=1.8, label="I")
axes[0].set_ylabel("activity")
axes[0].set_ylim(-0.05, 1.05)
axes[0].set_title("(a) Deterministic dynamics", loc="left")
axes[0].legend(loc="upper right", ncol=2)

axes[1].plot(t_sto, E_sto, color=COLORS["excitatory"], lw=1.0, label="E")
axes[1].plot(t_sto, I_sto, color=COLORS["inhibitory"], lw=1.0, label="I")
axes[1].set_xlabel("time")
axes[1].set_ylabel("activity")
axes[1].set_ylim(-0.05, 1.05)
axes[1].set_title(
    f"(b) Stochastic dynamics  "
    f"($\\sigma_E = \\sigma_I = {noise.sigma_E}$)", loc="left")
axes[1].legend(loc="upper right", ncol=2)

fig1.suptitle(
    "Adding Wiener noise turns a single deterministic orbit "
    "into a sample from a distribution",
    fontsize=13.5, y=1.00,
)
fig1.tight_layout()
fig1.savefig(ROOT / "notebooks" / "fig_stochastic_vs_deterministic.png")


# --------------------------------------------------------------------------- #
# Figure 2 — ensemble of stochastic trajectories
# --------------------------------------------------------------------------- #
n_trials = 50
t_e, E_ens, I_ens = ensemble_simulate(
    params, noise, n_trials=n_trials,
    E0=E0, I0=I0, t_max=t_max, dt=0.01, seed=42)

fig2, axes = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)

for k in range(n_trials):
    axes[0].plot(t_e, E_ens[k], color=COLORS["excitatory"],
                 lw=0.5, alpha=0.18)
    axes[1].plot(t_e, I_ens[k], color=COLORS["inhibitory"],
                 lw=0.5, alpha=0.18)

# Ensemble mean ± std for both populations
mean_E, std_E = E_ens.mean(axis=0), E_ens.std(axis=0)
mean_I, std_I = I_ens.mean(axis=0), I_ens.std(axis=0)

axes[0].plot(t_e, mean_E, color=COLORS["excitatory"], lw=2.0,
             label="ensemble mean")
axes[0].fill_between(t_e, mean_E - std_E, mean_E + std_E,
                      color=COLORS["excitatory"], alpha=0.18,
                      label="± 1 std")
axes[0].set_ylabel("E (excitatory)")
axes[0].set_ylim(-0.05, 1.05)
axes[0].set_title(f"(a) {n_trials} stochastic runs from identical "
                  f"initial condition", loc="left")
axes[0].legend(loc="upper right")

axes[1].plot(t_e, mean_I, color=COLORS["inhibitory"], lw=2.0,
             label="ensemble mean")
axes[1].fill_between(t_e, mean_I - std_I, mean_I + std_I,
                      color=COLORS["inhibitory"], alpha=0.18,
                      label="± 1 std")
axes[1].set_xlabel("time")
axes[1].set_ylabel("I (inhibitory)")
axes[1].set_ylim(-0.05, 1.05)
axes[1].set_title("(b) Same ensemble, inhibitory population", loc="left")
axes[1].legend(loc="upper right")

fig2.suptitle("Identical inputs, different outcomes: the stochastic ensemble",
              fontsize=13.5, y=1.00)
fig2.tight_layout()
fig2.savefig(ROOT / "notebooks" / "fig_stochastic_ensemble.png")


# --------------------------------------------------------------------------- #
# Figure 3 — phase-plane scatter of ensemble end states
# --------------------------------------------------------------------------- #
fig3, ax = plt.subplots(figsize=(7.5, 7.5))

# Show the deterministic trajectory as context
ax.plot(E_det, I_det, color=COLORS["muted"], lw=1.2, alpha=0.9,
        label="deterministic trajectory")

# End points of all stochastic trials
end_E = E_ens[:, -1]
end_I = I_ens[:, -1]
ax.scatter(end_E, end_I, s=40, color=COLORS["unstable"],
           edgecolor="white", linewidth=0.8, alpha=0.85,
           label=f"end states ({n_trials} runs)")

# Common initial condition
ax.scatter([E0], [I0], s=120, color=COLORS["accent"],
           edgecolor="white", linewidth=1.5, zorder=5,
           label="shared initial condition")

ax.set_xlabel("E (excitatory activity)")
ax.set_ylabel("I (inhibitory activity)")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect("equal")
ax.legend(loc="upper right")
ax.set_title("Ensemble end-states cluster around the deterministic "
             "fixed point", loc="left")

fig3.suptitle("Phase-plane view of the stochastic ensemble at t = "
              f"{t_max:.0f}", fontsize=13.5, y=1.00)
fig3.tight_layout()
fig3.savefig(ROOT / "notebooks" / "fig_stochastic_endstate_scatter.png")


print(f"Done. Three figures will open together. Close them to end.")
print(f"Ensemble end-state mean: E = {end_E.mean():.4f}, "
      f"I = {end_I.mean():.4f}")
print(f"Ensemble end-state std:  E = {end_E.std():.4f}, "
      f"I = {end_I.std():.4f}")
plt.show()
