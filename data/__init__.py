"""
Phasewalk — data module.

Loading and structuring intracranial EEG recordings, plus computing
functional connectivity and converting FCMs into graph-usable form.
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
from .connectivity import compute_fcm, fcm_to_graph

__all__ = [
    # loader
    "Seizure",
    "LoaderError",
    "load_seizure",
    "load_patient",
    "FS_HZ",
    "PRE_ICTAL_LENGTH_S",
    "POST_ICTAL_LENGTH_S",
    "MIN_ELECTRODES",
    "MAX_ELECTRODES",
    # connectivity
    "compute_fcm",
    "fcm_to_graph",
]
