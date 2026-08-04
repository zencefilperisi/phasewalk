"""
Phasewalk — data module.

Loading and structuring intracranial EEG recordings from the
SWEC-ETHZ short-term dataset.
"""
from .loader import (
    FS_HZ,
    LoaderError,
    MAX_ELECTRODES,
    MIN_ELECTRODES,
    POST_ICTAL_LENGTH_S,
    PRE_ICTAL_LENGTH_S,
    Seizure,
    load_patient,
    load_seizure,
)

__all__ = [
    "Seizure",
    "LoaderError",
    "load_seizure",
    "load_patient",
    "FS_HZ",
    "PRE_ICTAL_LENGTH_S",
    "POST_ICTAL_LENGTH_S",
    "MIN_ELECTRODES",
    "MAX_ELECTRODES",
]
