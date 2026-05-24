"""
Detect critical slowing down as the system approaches the Hopf bifurcation.

When a stochastic dynamical system is driven slowly toward a
bifurcation, its return to equilibrium after perturbation becomes
sluggish. Two signatures of this slowing-down are:

  - increasing variance of the noise-induced fluctuations,
  - increasing lag-1 autocorrelation of the time series.

These are the standard early-warning indicators of an impending
critical transition, and they are exactly what clinical neurology
looks for in the pre-ictal EEG.

This script computes both indicators over an ensemble of stochastic
trajectories at several P values approaching the Hopf threshold
(P_Hopf ≈ 0.753).

Run as:
    python notebooks/detect_critical_slowing.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np

from classical.model import WCParams
from classical.stochastic import (
    StochasticParams, ensemble_simulate,
)
from _style import COLORS, SEQ_CMAP, apply_style

apply_style()


# --------------------------------------------------------------------------- #
# Setup — sweep P values approaching the Hopf threshold from below
# --------------------------------------------------------------------------- #
P_VALUES = [0.30, 0.50, 0.65, 0.72]   # all below P_Hopf ≈ 0.753
N_TRIALS = 30
T_MAX = 200.0
DT = 0.01
NOISE = StochasticParams(sigma_E=0.025, sigma_I=0.025)

# Discard a transient at the start so we measure equilibrium fluctuations
TRANSIENT_FRAC = 0.4

print("Running stochastic ensembles at each P value...")
print("(This takes about 30 seconds total.)")


def equilibrium_statistics(ens):
    """Per-trial variance and lag-1 autocorrelation, averaged over trials."""
    n_steps = ens.shape[1]
    cut = int(n_steps * TRANSIENT_FRAC)
    tails = ens[:, cut:]                          # shape (n_trials, n_tail)

    variances = tails.var(axis=1)                 # per trial
    # Lag-1 autocorrelation per trial
    x = tails[:, :-1] - tails[:, :-1].mean(axis=1, keepdims=True)
    y = tails[:, 1:]  - tails[:, 1:].mean(axis=1, keepdims=True)
    num = (x * y).sum(axis=1)
    denom = np.sqrt((x ** 2).sum(axis=1) * (y ** 2).sum(axis=1))
    ac1 = np.where(denom > 0, num / denom, np.nan)

    return variances, ac1


results = []
example_traces = []   # keep one trace per P for the panel plot
for P in P_VALUES:
    params = WCParams(P=P)
    t, E_ens, _ = ensemble_simulate(
        params, NOISE, n_trials=N_TRIALS,
        E0=0.05, I0=0.02, t_max=T_MAX, dt=DT, seed=int(1000 * P))
    var, ac1 = equilibrium_statistics(E_ens)
    results.append({
        "P": P,
        "variance_mean": var.mean(), "variance_sem": var.std() / np.sqrt(N_TRIALS),
        "ac1_mean": ac1.mean(),       "ac1_sem":      ac1.std() / np.sqrt(N_TRIALS),
    })
    example_traces.append((P, t, E_ens[0]))


# --------------------------------------------------------------------------- #
# Figure 1 — example traces at increasing P
# --------------------------------------------------------------------------- #
# Use a discrete sample of the cividis palette so each P gets a unique
# colour on a perceptually uniform scale.
cmap = plt.get_cmap(SEQ_CMAP)
trace_colors = [cmap(0.15 + 0.7 * k / (len(P_VALUES) - 1))
                for k in range(len(P_VALUES))]

fig1, axes = plt.subplots(len(P_VALUES), 1, figsize=(12, 8.5), sharex=True)
for ax, (P, t, trace), c in zip(axes, example_traces, trace_colors):
    ax.plot(t, trace, color=c, lw=0.8)
    ax.set_ylim(-0.05, 0.35)
    ax.set_ylabel("E")
    ax.text(0.012, 0.82, f"$P = {P}$",
            transform=ax.transAxes, fontsize=11,
            color=COLORS["muted"],
            bbox=dict(facecolor=COLORS["background"],
                      edgecolor=COLORS["grid"],
                      pad=4, boxstyle="round,pad=0.35"))

axes[-1].set_xlabel("time")
fig1.suptitle(
    "Fluctuation amplitude grows as the system approaches "
    f"the Hopf threshold ($P_{{\\mathrm{{Hopf}}}} \\approx 0.753$)",
    fontsize=13.5, y=1.00,
)
fig1.tight_layout()
fig1.savefig(ROOT / "notebooks" / "fig_csd_traces.png")


# --------------------------------------------------------------------------- #
# Figure 2 — variance and lag-1 autocorrelation vs P
# --------------------------------------------------------------------------- #
Ps = np.array([r["P"] for r in results])
var_mean = np.array([r["variance_mean"] for r in results])
var_sem  = np.array([r["variance_sem"]  for r in results])
ac_mean  = np.array([r["ac1_mean"]      for r in results])
ac_sem   = np.array([r["ac1_sem"]       for r in results])

fig2, axes = plt.subplots(1, 2, figsize=(13, 5))

axes[0].errorbar(Ps, var_mean, yerr=var_sem,
                 fmt="o-", color=COLORS["stable"],
                 markersize=8, linewidth=1.6, capsize=4,
                 markerfacecolor=COLORS["stable"],
                 markeredgecolor="white", markeredgewidth=1.2)
axes[0].set_xlabel("$P$  (external drive)")
axes[0].set_ylabel("variance of E (equilibrium)")
axes[0].set_title("(a) Variance grows toward the threshold", loc="left")
axes[0].axvline(0.753, color=COLORS["accent"], ls="--", lw=1.3,
                alpha=0.8, label="$P_{\\mathrm{Hopf}} \\approx 0.753$")
axes[0].legend(loc="upper left")

axes[1].errorbar(Ps, ac_mean, yerr=ac_sem,
                 fmt="s-", color=COLORS["unstable"],
                 markersize=8, linewidth=1.6, capsize=4,
                 markerfacecolor=COLORS["unstable"],
                 markeredgecolor="white", markeredgewidth=1.2)
axes[1].set_xlabel("$P$  (external drive)")
axes[1].set_ylabel("lag-1 autocorrelation of E")
axes[1].set_title("(b) Autocorrelation rises toward the threshold", loc="left")
axes[1].axvline(0.753, color=COLORS["accent"], ls="--", lw=1.3,
                alpha=0.8, label="$P_{\\mathrm{Hopf}} \\approx 0.753$")
axes[1].legend(loc="upper left")

fig2.suptitle(
    "Critical slowing down: two early-warning signals "
    "of an impending Hopf bifurcation",
    fontsize=13.5, y=1.00,
)
fig2.tight_layout()
fig2.savefig(ROOT / "notebooks" / "fig_csd_indicators.png")


# --------------------------------------------------------------------------- #
# Print a tidy summary
# --------------------------------------------------------------------------- #
print("\nEarly-warning indicators across P values:")
print(f"{'P':>6}  {'variance':>14}  {'lag-1 AC':>14}")
for r in results:
    print(f"{r['P']:>6.2f}  "
          f"{r['variance_mean']:>9.5f} ± {r['variance_sem']:.5f}  "
          f"{r['ac1_mean']:>7.3f} ± {r['ac1_sem']:.3f}")

print("\nDone. Two figures will open. Close them to end.")
plt.show()
