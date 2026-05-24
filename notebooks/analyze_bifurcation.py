"""
Locate and characterize the Hopf bifurcation of the Wilson-Cowan model.

This script:

  1. Lists fixed points and their stability classes at four representative
     values of the external drive P.
  2. Locates the Hopf bifurcation threshold by bisection on the dominant
     eigenvalue.
  3. Sweeps P over [0, 2] and produces:
       (a) the bifurcation diagram with branches separated by stability
           class and the limit-cycle envelope from direct simulation;
       (b) the dominant eigenvalue's real part as a function of P,
           also separated by branch.

All figures are produced first and shown together at the end.

Run as:
    python notebooks/analyze_bifurcation.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np

from classical.analysis import find_fixed_points, find_hopf_threshold
from classical.model import WCParams, simulate
from _style import COLORS, apply_style

apply_style()
np.random.seed(42)


# --------------------------------------------------------------------------- #
# Step 1 — fixed-point inventory at selected P values
# --------------------------------------------------------------------------- #
print("=" * 64)
print("Fixed-point inventory at selected P values")
print("=" * 64)

for P in [0.5, 0.9, 1.2, 1.6]:
    params = WCParams(P=P)
    fps = find_fixed_points(params)
    print(f"\nP = {P}:  found {len(fps)} fixed point(s)")
    for k, fp in enumerate(fps):
        print(f"  FP{k}: (E*, I*) = ({fp.E:.4f}, {fp.I:.4f})")
        print(f"        eigenvalues = {fp.eigenvalues}")
        print(f"        type        = {fp.classification}")
        print(f"        max Re(λ)   = {fp.max_real_part:+.4f}")


# --------------------------------------------------------------------------- #
# Step 2 — locate the Hopf threshold by bisection
# --------------------------------------------------------------------------- #
print("\n" + "=" * 64)
print("Hopf bifurcation threshold (bisection on P)")
print("=" * 64)

base = WCParams()
P_hopf = find_hopf_threshold(base, param_name="P", lo=0.0, hi=2.0)
if P_hopf is not None:
    print(f"P_Hopf ≈ {P_hopf:.5f}")
else:
    print("No sign change detected — adjust the bracket.")


# --------------------------------------------------------------------------- #
# Step 3 — sweep P, collect branches and the limit-cycle envelope
# --------------------------------------------------------------------------- #
print("\nSweeping P for the bifurcation diagram...")
print("(Roughly 30 seconds — fixed-point search is the bottleneck.)")

P_values = np.linspace(0.0, 2.0, 81)

# Group fixed points by classification, both for E* (bifurcation diagram)
# and for max Re(λ) (eigenvalue panel). Plotting each class as a
# separate scatter avoids the zig-zag that arises when bistability gives
# multiple fixed points at one P.
branches = {cls: [] for cls in
            ("stable focus", "stable node", "saddle",
             "unstable focus", "unstable node")}
eig_traces = {cls: [] for cls in branches}

cycle_P, cycle_min, cycle_max = [], [], []

for P in P_values:
    params = WCParams(P=P)
    fps = find_fixed_points(params)
    for fp in fps:
        if fp.classification in branches:
            branches[fp.classification].append((P, fp.E))
            eig_traces[fp.classification].append((P, fp.max_real_part))

    # Limit cycle envelope: simulate beyond the Hopf threshold and read
    # the asymptotic min/max of E from the tail of the trajectory.
    if P_hopf is not None and P > P_hopf:
        _, E_traj, _ = simulate(params, E0=0.1, I0=0.05,
                                t_max=200, dt=0.05)
        tail = E_traj[-2000:]
        cycle_P.append(P)
        cycle_min.append(tail.min())
        cycle_max.append(tail.max())

# Detect the bistability window: a P-range where a stable focus or node
# coexists with an unstable focus. Used to add a shaded region to plots.
stable_Ps = set()
for cls in ("stable focus", "stable node"):
    stable_Ps.update(P for P, _ in branches[cls])
unstable_Ps = set()
for cls in ("unstable focus", "unstable node"):
    unstable_Ps.update(P for P, _ in branches[cls])
bistable_Ps = sorted(stable_Ps & unstable_Ps)
bistable_lo = min(bistable_Ps) if bistable_Ps else None
bistable_hi = max(bistable_Ps) if bistable_Ps else None


# --------------------------------------------------------------------------- #
# Marker style table — consistent across all panels
# --------------------------------------------------------------------------- #
class_styles = {
    "stable focus":   dict(color=COLORS["stable"],   marker="o",
                           filled=True,  s=22, label="Stable focus"),
    "stable node":    dict(color=COLORS["stable"],   marker="s",
                           filled=True,  s=22, label="Stable node"),
    "saddle":         dict(color=COLORS["saddle"],   marker="x",
                           filled=True,  s=28, label="Saddle"),
    "unstable focus": dict(color=COLORS["unstable"], marker="o",
                           filled=False, s=22, label="Unstable focus"),
    "unstable node":  dict(color=COLORS["unstable"], marker="s",
                           filled=False, s=22, label="Unstable node"),
}


def _scatter_class(ax, points, style):
    """Plot one stability class on `ax` using its style spec."""
    if not points:
        return
    Ps, ys = zip(*points)
    if style["filled"]:
        ax.scatter(Ps, ys, color=style["color"], marker=style["marker"],
                   s=style["s"], label=style["label"], zorder=4)
    else:
        ax.scatter(Ps, ys, edgecolors=style["color"], facecolors="none",
                   marker=style["marker"], s=style["s"], linewidths=1.2,
                   label=style["label"], zorder=4)


# --------------------------------------------------------------------------- #
# Figure 1 — bifurcation diagram
# --------------------------------------------------------------------------- #
fig1, ax = plt.subplots(figsize=(11, 6.5))

# Bistability window — shaded background
if bistable_lo is not None:
    ax.axvspan(bistable_lo, bistable_hi,
               color=COLORS["saddle"], alpha=0.08, zorder=0)
    ax.text((bistable_lo + bistable_hi) / 2, 0.78,
            "Bistability\nwindow",
            ha="center", va="center", fontsize=10,
            color=COLORS["saddle"],
            bbox=dict(facecolor=COLORS["background"],
                      edgecolor=COLORS["grid"],
                      pad=4, boxstyle="round,pad=0.4"))

# Limit-cycle envelope first, so dots sit on top of the band
if cycle_P:
    ax.fill_between(cycle_P, cycle_min, cycle_max,
                    color=COLORS["limit_cycle"], alpha=0.10,
                    label="Limit-cycle envelope", zorder=1)
    ax.plot(cycle_P, cycle_max, color=COLORS["limit_cycle"], lw=1.4)
    ax.plot(cycle_P, cycle_min, color=COLORS["limit_cycle"], lw=1.4)

# Fixed-point branches as scatter
for cls, points in branches.items():
    _scatter_class(ax, points, class_styles[cls])

# Hopf marker
if P_hopf is not None:
    ax.axvline(P_hopf, color=COLORS["accent"], ls="--", lw=1.3,
               alpha=0.85,
               label=f"Hopf threshold  ($P \\approx {P_hopf:.3f}$)")

ax.set_xlabel("$P$  (external drive to E)")
ax.set_ylabel("$E^*$  (excitatory activity at fixed point or cycle extrema)")
ax.set_xlim(P_values.min(), P_values.max())
ax.set_ylim(-0.05, 1.05)

ax.set_title("Subcritical Hopf bifurcation produces a bistability window "
             "before sustained oscillation",
             loc="left")
ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
          frameon=False, fontsize=9.5)

fig1.suptitle("Bifurcation diagram of the Wilson-Cowan model",
              fontsize=13.5, y=1.00)
fig1.tight_layout()
fig1.savefig(ROOT / "notebooks" / "fig_bifurcation_diagram.png")


# --------------------------------------------------------------------------- #
# Figure 2 — eigenvalue trace, separated by branch
# --------------------------------------------------------------------------- #
print("Building eigenvalue tracking plot...")

fig2, ax = plt.subplots(figsize=(11, 5.2))

# Bistability window again
if bistable_lo is not None:
    ax.axvspan(bistable_lo, bistable_hi,
               color=COLORS["saddle"], alpha=0.08, zorder=0)

ax.axhline(0.0, color=COLORS["muted"], lw=0.9, ls="--", alpha=0.7)

for cls, points in eig_traces.items():
    _scatter_class(ax, points, class_styles[cls])

if P_hopf is not None:
    ax.axvline(P_hopf, color=COLORS["accent"], ls="--", lw=1.3,
               alpha=0.85,
               label=f"Hopf threshold  ($P \\approx {P_hopf:.3f}$)")

ax.set_xlabel("$P$")
ax.set_ylabel(r"max  $\mathrm{Re}(\lambda)$  per fixed point")
ax.set_title(
    "Stable branch loses stability at the Hopf point; "
    "saddle and unstable branches appear together",
    loc="left",
)
ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
          frameon=False, fontsize=9.5)

fig2.suptitle("Eigenvalue trace as the external drive varies",
              fontsize=13.5, y=1.00)
fig2.tight_layout()
fig2.savefig(ROOT / "notebooks" / "fig_eigenvalue_trace.png")


print("\nDone. Two figures will open together. Close them to end.")
plt.show()
