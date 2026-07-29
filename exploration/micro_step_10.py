"""
Micro-step 10: First look at real iEEG data.

We are about to open, for the first time, a real intracranial EEG
recording from a patient with drug-resistant epilepsy. The point of
this micro-step is NOT to analyse anything yet. It is to:

  (1) load the .mat file successfully,
  (2) inspect what is inside,
  (3) verify the shape and orientation of the EEG array,
  (4) plot a small window of a few channels so we SEE the signal.

Nothing is more important, when starting on real data, than looking
at it before touching it with any analysis.

Fill in the TODOs. Run from anywhere; edit MAT_PATH to point at a
real .mat file on your disk.
"""
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------- #
# TODO 1: set MAT_PATH to the absolute path of ONE seizure file on your disk.
#
# Example (Windows):
#   MAT_PATH = r"C:\Users\User\Desktop\SWEC\ID1\Sz1.mat"
#
# The 'r' before the string is important on Windows so that backslashes
# are not interpreted as escape codes. Alternative: use forward slashes.
#
MAT_PATH = r"C:\Users\User\Desktop\ID1\Sz1.mat"   # <-- TODO 1
# --------------------------------------------------------------------------- #

if MAT_PATH is None:
    raise SystemExit("TODO 1 not done: set MAT_PATH first.")


# --------------------------------------------------------------------------- #
# Load the .mat file. loadmat returns a dict-like object with the
# variables that were saved from MATLAB, plus some private '__..__' keys
# we can ignore.
# --------------------------------------------------------------------------- #
mat = sio.loadmat(MAT_PATH)

print("Keys inside the .mat file:")
for k in mat.keys():
    # skip MATLAB's internal bookkeeping keys
    if k.startswith("__"):
        continue
    v = mat[k]
    print(f"  {k!r}: type={type(v).__name__}, shape={getattr(v, 'shape', 'N/A')}, "
          f"dtype={getattr(v, 'dtype', 'N/A')}")


# --------------------------------------------------------------------------- #
# TODO 2: grab the EEG array.
#
# According to the SWEC-ETHZ documentation, the short-term dataset stores
# the signal in a variable called 'EEG'. Retrieve it from the loaded dict.
#
#   EEG = mat['EEG']
#
EEG = mat['EEG']   # <-- TODO 2
# --------------------------------------------------------------------------- #

if EEG is None:
    raise SystemExit("TODO 2 not done: fetch mat['EEG'].")

print(f"\nEEG array: shape={EEG.shape}, dtype={EEG.dtype}")


# --------------------------------------------------------------------------- #
# Orientation check. Documentation says:
#   short-term:  T x M    with T = time samples, M = electrodes
# Number of electrodes for SWEC-ETHZ short-term is 36-100.
# Number of time samples for a ~6-16 min recording at 512 Hz is
# hundreds of thousands. So the LARGER axis must be time, the SMALLER
# must be electrodes.
# --------------------------------------------------------------------------- #
n_rows, n_cols = EEG.shape
if n_rows > n_cols:
    n_time, n_channels = n_rows, n_cols
    axis_time = 0
    print(f"\nOrientation looks correct (T x M):")
else:
    # Unexpected, but handle it gracefully by transposing.
    print(f"\nOrientation looks TRANSPOSED — swapping axes so time is rows.")
    EEG = EEG.T
    n_time, n_channels = EEG.shape
    axis_time = 0

print(f"  time samples: {n_time}")
print(f"  electrodes:   {n_channels}")

# Sanity: at 512 Hz, how many seconds is this?
fs = 512
duration_s = n_time / fs
print(f"  duration: {duration_s:.1f} s  (={duration_s/60:.2f} min)")

# For SWEC-ETHZ short-term, seizure onset is at sample 512*3*60 = 92160.
seizure_onset_sample = fs * 3 * 60
seizure_end_sample = n_time - fs * 3 * 60
seizure_duration_s = (seizure_end_sample - seizure_onset_sample) / fs
print(f"  seizure onset at sample {seizure_onset_sample} (t = 3.00 min)")
print(f"  seizure ends at sample {seizure_end_sample} (t = {seizure_end_sample/fs/60:.2f} min)")
print(f"  seizure duration: {seizure_duration_s:.1f} s")


# --------------------------------------------------------------------------- #
# TODO 3: pick FOUR channels to plot. Any four indices will do; a
# reasonable pick is the first four: channels 0, 1, 2, 3.
#
# Fill in a list of four integer channel indices:
#
channels_to_plot = [0, 1, 2, 3]   # <-- TODO 3   (e.g. [0, 1, 2, 3])
# --------------------------------------------------------------------------- #

if channels_to_plot is None:
    raise SystemExit("TODO 3 not done: pick four channel indices.")


# --------------------------------------------------------------------------- #
# Plot the four channels across the whole recording, with vertical lines
# marking the seizure onset and end.
# --------------------------------------------------------------------------- #
t = np.arange(n_time) / fs   # time axis in seconds

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
plt.savefig("micro_step_10_first_look.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_10_first_look.png")
print("\nLook at the plot. You should see:")
print("  - relatively calm signal in the first 3 min (pre-ictal)")
print("  - a change in character between the two dashed red lines (ictal)")
print("  - post-ictal signal after the second dashed line")
print("  - some channels may show the seizure more clearly than others")