"""
Noise-induced transitions in the bistability window.

In the bistability window (0.75 < P < 0.95) the deterministic Wilson-Cowan
model has a stable focus near low activity coexisting with a large-amplitude
limit cycle. Without noise, the system stays in whichever basin it started.

With noise, small fluctuations can push the trajectory across the saddle
separating the basins, triggering a sudden jump into the limit cycle —
the mathematical analogue of a seizure that arrives without warning in a
clinically "normal" patient.

This script:
  (1) Runs many stochastic trajectories from the low-activity basin and
      records how many escape into the large-amplitude cycle within a
      fixed time window.
  (2) Plots example traces showing the abrupt transition.
  (3) Quantifies escape probability as a function of noise amplitude.

Run as:
    python notebooks/noise_induced_transitions.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np

from classical.model import WCParams
from classical.stochastic import StochasticParams, ensemble_simulate
from _style import COLORS, SEQ_CMAP, apply_style

apply_style()


# --------------------------------------------------------------------------- #
# Setup — operate inside the bistability window
# --------------------------------------------------------------------------- #
# At P = 0.85 the deterministic system has three fixed points: a
# stable focus near (E*, I*) ≈ (0.04, 0.003), a saddle, and an
# unstable focus surrounded by a large-amplitude limit cycle.
P = 0.85
params = WCParams(P=P)

# Start near the low-activity stable focus
E0, I0 = 0.04, 0.003
T_MAX = 250.0
DT = 0.01
ESCAPE_THRESHOLD = 0.5   # E exceeding this counts as "in the seizure cycle"

# Noise sweep — same physics, just louder
NOISE_LEVELS = [0.005, 0.015, 0.030, 0.050]
N_TRIALS = 40


def fraction_escaped(E_ensemble, threshold=ESCAPE_THRESHOLD):
    """Fraction of trajectories that exceed `threshold` at some point."""
    return np.mean(E_ensemble.max(axis=1) > threshold)


# --------------------------------------------------------------------------- #
# Run ensembles at each noise level
# --------------------------------------------------------------------------- #
print(f"Running noise-induced transition experiments at P = {P}...")
print(f"({len(NOISE_LEVELS)} noise levels × {N_TRIALS} trials each.)\n")

results = []
example_traces = {}
for sigma in NOISE_LEVELS:
    noise = StochasticParams(sigma_E=sigma, sigma_I=sigma)
    t, E_ens, _ = ensemble_simulate(
        params, noise, n_trials=N_TRIALS,
        E0=E0, I0=I0, t_max=T_MAX, dt=DT,
        seed=int(10_000 * sigma))
    frac = fraction_escaped(E_ens)
    results.append((sigma, frac))
    # Keep one escaped + one staying trace if available, otherwise any
    escaped = np.where(E_ens.max(axis=1) > ESCAPE_THRESHOLD)[0]
    stayed = np.where(E_ens.max(axis=1) <= ESCAPE_THRESHOLD)[0]
    ex_escape = E_ens[escaped[0]] if escaped.size else None
    ex_stay = E_ens[stayed[0]] if stayed.size else None
    example_traces[sigma] = (t, ex_escape, ex_stay)

    print(f"  sigma = {sigma:.3f}:  escape probability = {frac:.2%}")


# --------------------------------------------------------------------------- #
# Figure 1 — example traces at each noise level
# --------------------------------------------------------------------------- #
fig1, axes = plt.subplots(len(NOISE_LEVELS), 1, figsize=(12, 9), sharex=True)

for ax, sigma in zip(axes, NOISE_LEVELS):
    t, esc, sta = example_traces[sigma]
    if sta is not None:
        ax.plot(t, sta, color=COLORS["stable"], lw=0.8, alpha=0.9,
                label="trial that stayed")
    if esc is not None:
        ax.plot(t, esc, color=COLORS["unstable"], lw=0.8, alpha=0.95,
                label="trial that escaped")
    ax.axhline(ESCAPE_THRESHOLD, color=COLORS["muted"], lw=0.8, ls=":",
               alpha=0.7)
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel("E")
    ax.text(0.012, 0.82, f"$\\sigma = {sigma}$",
            transform=ax.transAxes, fontsize=11,
            color=COLORS["muted"],
            bbox=dict(facecolor=COLORS["background"],
                      edgecolor=COLORS["grid"],
                      pad=4, boxstyle="round,pad=0.35"))
    ax.legend(loc="upper right", fontsize=9)

axes[-1].set_xlabel("time")
fig1.suptitle(
    f"Noise-induced escape from the low-activity basin "
    f"(P = {P}, inside bistability window)",
    fontsize=13.5, y=1.00,
)
fig1.tight_layout()
fig1.savefig(ROOT / "notebooks" / "fig_noise_traces.png")


# --------------------------------------------------------------------------- #
# Figure 2 — escape probability vs noise amplitude
# --------------------------------------------------------------------------- #
sigmas = np.array([r[0] for r in results])
fracs = np.array([r[1] for r in results])

fig2, ax = plt.subplots(figsize=(9, 5.2))
ax.plot(sigmas, fracs, "o-", color=COLORS["unstable"],
        markersize=9, lw=1.8,
        markerfacecolor=COLORS["unstable"],
        markeredgecolor="white", markeredgewidth=1.2)
ax.set_xlabel("noise amplitude $\\sigma$")
ax.set_ylabel(f"escape probability in $t \\leq {T_MAX:.0f}$")
ax.set_ylim(-0.03, 1.03)
ax.set_title(
    "Larger fluctuations make spontaneous transitions more likely",
    loc="left")

fig2.suptitle(
    f"Escape probability vs noise amplitude  "
    f"($P = {P}$, {N_TRIALS} trials per point)",
    fontsize=13.5, y=1.00,
)
fig2.tight_layout()
fig2.savefig(ROOT / "notebooks" / "fig_noise_escape_probability.png")


print("\nDone. Two figures will open. Close them to end.")
plt.show()
