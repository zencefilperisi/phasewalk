"""
Micro-step 13: First real functional connectivity matrix.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib.pyplot as plt

from data.loader import load_seizure


MAT_PATH = r"C:\Users\User\Desktop\SWEC\ID1\Sz1.mat"

sz = load_seizure(MAT_PATH)
print(f"Loaded {sz.patient_id}/{sz.seizure_id}")
print(f"  shape: {sz.eeg.shape}     (electrodes, time samples)")
print(f"  ictal duration: {sz.ictal_duration_s:.1f} s")

ictal = sz.segment("ictal")
n_electrodes, n_samples = ictal.shape
print(f"\nIctal segment: {ictal.shape}  "
      f"({n_samples / sz.fs:.1f} s at {sz.fs} Hz)")

fcm = np.corrcoef(ictal)

print(f"\nFCM shape:           {fcm.shape}")
print(f"FCM symmetric:       {np.allclose(fcm, fcm.T)}")
print(f"FCM diagonal ~= 1:   {np.allclose(np.diag(fcm), 1.0)}")
print(f"FCM values in [-1,1]:{np.all((fcm >= -1.001) & (fcm <= 1.001))}")

mask = ~np.eye(n_electrodes, dtype=bool)
off_diag = fcm[mask]
print(f"\nOff-diagonal correlations:")
print(f"  min:    {off_diag.min():+.3f}")
print(f"  max:    {off_diag.max():+.3f}")
print(f"  mean:   {off_diag.mean():+.3f}")
print(f"  median: {np.median(off_diag):+.3f}")

fig, axes = plt.subplots(1, 2, figsize=(15, 6),
                         gridspec_kw={"width_ratios": [1, 1.05]})

im = axes[0].imshow(fcm, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
axes[0].set_xlabel("channel")
axes[0].set_ylabel("channel")
axes[0].set_title(f"Pearson FCM — {sz.patient_id}/{sz.seizure_id}  "
                  f"(ictal, {sz.ictal_duration_s:.0f} s)")
plt.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04,
             label="correlation")

axes[1].hist(off_diag, bins=40, color="#4C72B0",
             edgecolor="white", alpha=0.9)
axes[1].axvline(0, color="black", lw=1, ls="--")
axes[1].axvline(off_diag.mean(), color="crimson", lw=1.5,
                label=f"mean = {off_diag.mean():+.3f}")
axes[1].set_xlabel("off-diagonal correlation")
axes[1].set_ylabel("count")
axes[1].set_title("Distribution of off-diagonal FCM entries")
axes[1].legend()

plt.tight_layout()
plt.savefig("figures/micro_step_13_graph.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_13_first_fcm.png")
