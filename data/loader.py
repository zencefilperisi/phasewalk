"""
Loader for SWEC-ETHZ short-term seizure files.

One entry point: `load_seizure(path)` reads a .mat file, corrects the
array orientation, verifies that the file is structurally sane, and
returns a `Seizure` dataclass with the signal and its pre-computed
segment boundaries.

Every fixed number in this module has a documented source:
  - fs = 512 Hz            : SWEC-ETHZ short-term docs (fixed sampling rate)
  - onset at 3 min          : SWEC-ETHZ short-term docs (fixed preictal length)
  - post-ictal length 3 min : SWEC-ETHZ short-term docs
  - 36 <= n_electrodes <=100: SWEC-ETHZ short-term docs (implantation range)

None of these values were tuned to ID1 or any other patient.
See ADR 0006 for the parameter-independence rule this follows.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union

import numpy as np
import scipy.io as sio


# --- Constants from the SWEC-ETHZ short-term documentation ---------------- #
# Source: http://ieeg-swez.ethz.ch/ , "Short-term Dataset" section.
# See also ADR 0003 for the choice of the short-term subset.
FS_HZ = 512
PRE_ICTAL_LENGTH_S = 3 * 60         # 3 minutes
POST_ICTAL_LENGTH_S = 3 * 60        # 3 minutes
MIN_ELECTRODES = 36
MAX_ELECTRODES = 100


@dataclass(frozen=True)
class Seizure:
    """A single seizure recording with pre-computed segment boundaries.

    Attributes
    ----------
    eeg : ndarray, shape (n_electrodes, n_samples)
        The iEEG signal, oriented so rows are electrodes and columns
        are time samples. This is the matrix-algebra natural
        orientation (each electrode is a "channel" indexed by row).
    fs : int
        Sampling frequency in Hz.
    n_electrodes, n_samples : int
        Convenient shape aliases.
    duration_s : float
        Total recording length in seconds.
    preictal_range, ictal_range, postictal_range : tuple[int, int]
        (start, end) sample indices for each segment, half-open
        intervals in the Python sense: eeg[:, start:end] gives the
        segment. Using these ranges is the ONLY correct way to slice
        the signal into peri-ictal epochs.
    ictal_duration_s : float
        Length of the ictal segment in seconds.
    patient_id, seizure_id : str
        Identifiers derived from the file path (e.g. "ID1", "Sz3").
    source_path : str
        Absolute path the file was loaded from (for traceability).
    """
    eeg: np.ndarray
    fs: int
    n_electrodes: int
    n_samples: int
    duration_s: float
    preictal_range: tuple
    ictal_range: tuple
    postictal_range: tuple
    ictal_duration_s: float
    patient_id: str
    seizure_id: str
    source_path: str

    def segment(self, name: str) -> np.ndarray:
        """Return the requested segment as (n_electrodes, n_samples_seg).

        `name` must be 'preictal', 'ictal', or 'postictal'.
        """
        ranges = {
            "preictal": self.preictal_range,
            "ictal": self.ictal_range,
            "postictal": self.postictal_range,
        }
        if name not in ranges:
            raise ValueError(
                f"segment name must be one of {list(ranges)}, got {name!r}"
            )
        start, end = ranges[name]
        return self.eeg[:, start:end]


class LoaderError(RuntimeError):
    """Raised when a .mat file cannot be interpreted as a Seizure."""


def _extract_ids(path: Path) -> tuple:
    """Infer (patient_id, seizure_id) from a path like .../ID1/Sz3.mat."""
    seizure_id = path.stem                    # "Sz3"
    patient_id = path.parent.name             # "ID1"
    return patient_id, seizure_id


def _orient_electrodes_first(eeg: np.ndarray) -> np.ndarray:
    """Return the array oriented so rows are electrodes, columns are time.

    The SWEC-ETHZ short-term docs say the on-disk layout is (T, M) with
    T = time samples and M = electrodes. The two axes are easily
    distinguished by size: T is hundreds of thousands, M is 36-100.
    We transpose to (electrodes, time) which is the matrix-algebra
    natural orientation.
    """
    if eeg.ndim != 2:
        raise LoaderError(f"EEG array must be 2-D, got ndim={eeg.ndim}")

    n_rows, n_cols = eeg.shape
    if MIN_ELECTRODES <= n_rows <= MAX_ELECTRODES and n_cols > n_rows:
        # already (electrodes, time)
        return eeg
    if MIN_ELECTRODES <= n_cols <= MAX_ELECTRODES and n_rows > n_cols:
        # (time, electrodes) on disk — transpose
        return eeg.T
    raise LoaderError(
        f"Cannot identify electrode axis: shape={eeg.shape}. "
        f"Expected one axis in [{MIN_ELECTRODES}, {MAX_ELECTRODES}]."
    )


def _health_checks(eeg: np.ndarray, n_samples: int) -> None:
    """Raise LoaderError if the signal is obviously corrupt."""
    if not np.isfinite(eeg).all():
        n_bad = int((~np.isfinite(eeg)).sum())
        raise LoaderError(f"EEG contains {n_bad} non-finite (NaN/Inf) values")

    # A channel that is exactly zero throughout the recording is
    # broken / unplugged, not a valid recording.
    zero_channels = np.where(np.all(eeg == 0.0, axis=1))[0]
    if len(zero_channels) > 0:
        raise LoaderError(
            f"channels are entirely zero: {zero_channels.tolist()}"
        )

    # Minimum sensible recording length: preictal (3 min) + at least
    # a few seconds of ictal + postictal (3 min). Anything shorter
    # cannot even contain the expected structure.
    min_samples = FS_HZ * (PRE_ICTAL_LENGTH_S + 1 + POST_ICTAL_LENGTH_S)
    if n_samples < min_samples:
        raise LoaderError(
            f"recording is too short: {n_samples} samples "
            f"(< minimum expected {min_samples})"
        )


def load_seizure(path: Union[str, Path]) -> Seizure:
    """Load one SWEC-ETHZ short-term seizure recording.

    Parameters
    ----------
    path : str or Path
        Path to a `.mat` file (e.g., ".../ID1/Sz3.mat").

    Returns
    -------
    Seizure

    Raises
    ------
    LoaderError
        If the file is missing, malformed, transposed in an
        unrecognized way, or contains obviously corrupt data.
    """
    path = Path(path).resolve()
    if not path.is_file():
        raise LoaderError(f"file does not exist: {path}")

    try:
        mat = sio.loadmat(str(path))
    except Exception as e:
        raise LoaderError(f"scipy.io.loadmat failed: {e}") from e

    if "EEG" not in mat:
        raise LoaderError(f"no 'EEG' variable in {path.name}")

    eeg = _orient_electrodes_first(mat["EEG"])
    n_electrodes, n_samples = eeg.shape
    _health_checks(eeg, n_samples)

    duration_s = n_samples / FS_HZ
    onset_sample = FS_HZ * PRE_ICTAL_LENGTH_S
    end_sample = n_samples - FS_HZ * POST_ICTAL_LENGTH_S
    ictal_duration_s = (end_sample - onset_sample) / FS_HZ

    if ictal_duration_s <= 0:
        raise LoaderError(
            f"implied ictal length is non-positive ({ictal_duration_s} s) — "
            f"file may not follow the SWEC-ETHZ short-term structure"
        )

    patient_id, seizure_id = _extract_ids(path)

    return Seizure(
        eeg=eeg,
        fs=FS_HZ,
        n_electrodes=n_electrodes,
        n_samples=n_samples,
        duration_s=duration_s,
        preictal_range=(0, onset_sample),
        ictal_range=(onset_sample, end_sample),
        postictal_range=(end_sample, n_samples),
        ictal_duration_s=ictal_duration_s,
        patient_id=patient_id,
        seizure_id=seizure_id,
        source_path=str(path),
    )


def load_patient(folder: Union[str, Path]) -> list:
    """Load all seizures for one patient folder.

    Returns a list of Seizure objects, sorted by seizure_id. Files
    that fail to load are skipped with a warning printed to stderr;
    they do not stop the batch.
    """
    import sys

    folder = Path(folder).resolve()
    if not folder.is_dir():
        raise LoaderError(f"patient folder does not exist: {folder}")

    seizures = []
    for f in sorted(folder.glob("Sz*.mat")):
        try:
            seizures.append(load_seizure(f))
        except LoaderError as e:
            print(f"WARNING: skipping {f.name}: {e}", file=sys.stderr)
    return seizures
