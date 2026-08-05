"""
Micro-step 10: First look at real iEEG data.
"""
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


MAT_PATH = r"C:\Users\User\Desktop\SWEC\ID1\Sz1.mat"

mat = sio.loadmat(MAT_PATH)

print("Keys inside the .mat file:")
for k in mat.keys():
    if k.startswith("__"):
        continue
    v = mat[k]
    print(f"  {k!r}: type={type(v).__name__}, shape={getattr(v, 'shape', 'N/A')}, "
          f"dtype={getattr(v, 'dtype', 'N/A')}")

EEG = mat['EEG']

print(f"\nEEG array: shape={EEG.shape}, dtype={EEG.dtype}")

n_rows, n_cols = EEG.shape
if n_rows > n_cols:
    n_time, n_channels = n_rows, n_cols
    print(f"\nOrientation looks correct (T x M):")
else:
    print(f"\nOrientation looks TRANSPOSED — swapping axes so time is rows.")
    EEG = EEG.T
    n_time, n_channels = EEG.shape

print(f"  time samples: {n_time}")
print(f"  electrodes:   {n_channels}")

fs = 512
duration_s = n_time / fs
print(f"  duration: {duration_s:.1f} s  (={duration_s/60:.2f} min)")

seizure_onset_sample = fs * 3 * 60
seizure_end_sample = n_time - fs * 3 * 60
seizure_duration_s = (seizure_end_sample - seizure_onset_sample) / fs
print(f"  seizure onset at sample {seizure_onset_sample} (t = 3.00 min)")
print(f"  seizure ends at sample {seizure_end_sample} (t = {seizure_end_sample/fs/60:.2f} min)")
print(f"  seizure duration: {seizure_duration_s:.1f} s")

channels_to_plot = [0, 1, 2, 3]

t = np.arange(n_time) / fs

fig, axes = plt.subplots(len(channels_to_plot), 1,
                          figsize=(13, 8), sharex=True)

for ax, ch in zip(axes, channels_to_plot):
    ax.plot(t, EEG[:, ch], lw=0.5, color="#4C72B0")
    ax.axvline(seizure_onset_sample / fs, color="crimson", lw=1.2,
               ls="--", label="seizure onset" if ch == channels_to_plot[0] else None)
    ax.axvline(seizure_end_sample / fs, color="crimson", lw=1.2,
               ls="--", label="seizure end" if ch == channels_to_plot[0] else None)
    ax.set_ylabel(f"ch {ch}")
    ax.grid(alpha=0.3)

axes[0].legend(loc="upper right")
axes[-1].set_xlabel("time (s)")
fig.suptitle(f"iEEG recording — {MAT_PATH.split(chr(92))[-1]}", fontsize=13)
plt.tight_layout()
plt.savefig("figures/micro_step_10_graph.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_10_first_look.png")
