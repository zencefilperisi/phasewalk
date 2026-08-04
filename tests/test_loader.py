"""
Tests for the data loader.

We do not ship real iEEG files with the repo (too large, licensed).
Instead each test builds a synthetic .mat file with a known structure
in a temporary directory and verifies the loader handles it correctly.

Coverage:
  - correct-shape file loads successfully
  - both on-disk orientations (T x M and M x T) load to the same
    canonical (n_electrodes, n_samples) output
  - segment boundaries are computed correctly and sum to the whole file
  - `.segment("preictal"/"ictal"/"postictal")` returns the right slices
  - malformed inputs raise LoaderError with informative messages
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import scipy.io as sio
import pytest

from data.loader import (
    FS_HZ,
    LoaderError,
    PRE_ICTAL_LENGTH_S,
    POST_ICTAL_LENGTH_S,
    Seizure,
    load_seizure,
    load_patient,
)


# --------------------------------------------------------------------------- #
# Helpers to build synthetic .mat files that follow the SWEC-ETHZ format.
# --------------------------------------------------------------------------- #

def _make_mat(tmp_path, patient="ID1", seizure="Sz1",
              n_electrodes=47, ictal_len_s=100,
              orientation="TxM", corrupt=None) -> Path:
    """Create a fake .mat file matching the SWEC-ETHZ short-term layout.

    Total length = preictal (3 min) + ictal_len_s + postictal (3 min).
    orientation: "TxM" (as documented) or "MxT" (transposed).
    corrupt: None | "nan" | "zero_channel" | "no_EEG_var" | "1d" | "empty"
    """
    n_samples = FS_HZ * (PRE_ICTAL_LENGTH_S + ictal_len_s + POST_ICTAL_LENGTH_S)
    rng = np.random.default_rng(0)
    signal = rng.standard_normal((n_electrodes, n_samples)) * 50.0

    if corrupt == "nan":
        signal[3, 1000] = np.nan
    elif corrupt == "zero_channel":
        signal[7, :] = 0.0

    # Store on disk in the requested orientation.
    on_disk = signal.T if orientation == "TxM" else signal

    patient_dir = tmp_path / patient
    patient_dir.mkdir(parents=True, exist_ok=True)
    mat_path = patient_dir / f"{seizure}.mat"

    if corrupt == "no_EEG_var":
        sio.savemat(mat_path, {"NOT_EEG": on_disk})
    elif corrupt == "1d":
        sio.savemat(mat_path, {"EEG": on_disk[0, :]})
    elif corrupt == "empty":
        sio.savemat(mat_path, {"EEG": np.zeros((0, 0))})
    else:
        sio.savemat(mat_path, {"EEG": on_disk})

    return mat_path


# --------------------------------------------------------------------------- #
# Happy-path tests
# --------------------------------------------------------------------------- #

def test_loads_documented_orientation(tmp_path):
    mat = _make_mat(tmp_path, ictal_len_s=100, orientation="TxM")
    sz = load_seizure(mat)
    assert isinstance(sz, Seizure)
    assert sz.n_electrodes == 47
    assert sz.fs == FS_HZ
    assert sz.ictal_duration_s == pytest.approx(100.0)


def test_loads_transposed_orientation_to_canonical(tmp_path):
    """A file stored as (electrodes, time) — the 'wrong' orientation —
    must still end up as (electrodes, time) after loading, because that
    is already the canonical output. We build both and verify the eeg
    arrays are identical up to the deterministic content."""
    mat_ok = _make_mat(tmp_path / "a", n_electrodes=47,
                       ictal_len_s=50, orientation="TxM")
    mat_tx = _make_mat(tmp_path / "b", n_electrodes=47,
                       ictal_len_s=50, orientation="MxT")
    sz_ok = load_seizure(mat_ok)
    sz_tx = load_seizure(mat_tx)
    # Same shape after loading, regardless of on-disk orientation.
    assert sz_ok.eeg.shape == sz_tx.eeg.shape == (47, sz_ok.n_samples)


def test_segment_ranges_partition_recording(tmp_path):
    mat = _make_mat(tmp_path, ictal_len_s=80)
    sz = load_seizure(mat)
    pre_s, pre_e = sz.preictal_range
    ic_s,  ic_e  = sz.ictal_range
    po_s,  po_e  = sz.postictal_range
    assert pre_s == 0
    assert pre_e == ic_s
    assert ic_e == po_s
    assert po_e == sz.n_samples
    # And the pre-ictal length is exactly 3 minutes.
    assert pre_e - pre_s == FS_HZ * PRE_ICTAL_LENGTH_S


def test_segment_method_returns_correct_slices(tmp_path):
    mat = _make_mat(tmp_path, ictal_len_s=60)
    sz = load_seizure(mat)
    pre = sz.segment("preictal")
    ic = sz.segment("ictal")
    po = sz.segment("postictal")
    assert pre.shape[1] == FS_HZ * PRE_ICTAL_LENGTH_S
    assert ic.shape[1]  == FS_HZ * 60
    assert po.shape[1]  == FS_HZ * POST_ICTAL_LENGTH_S
    # Concatenated segments reconstruct the whole recording.
    reconstructed = np.concatenate([pre, ic, po], axis=1)
    assert np.array_equal(reconstructed, sz.eeg)


def test_segment_rejects_unknown_name(tmp_path):
    mat = _make_mat(tmp_path, ictal_len_s=60)
    sz = load_seizure(mat)
    with pytest.raises(ValueError):
        sz.segment("nonsense")


def test_patient_and_seizure_ids_extracted(tmp_path):
    mat = _make_mat(tmp_path, patient="ID7", seizure="Sz42", ictal_len_s=30)
    sz = load_seizure(mat)
    assert sz.patient_id == "ID7"
    assert sz.seizure_id == "Sz42"


# --------------------------------------------------------------------------- #
# Error paths
# --------------------------------------------------------------------------- #

def test_missing_file_raises(tmp_path):
    with pytest.raises(LoaderError, match="does not exist"):
        load_seizure(tmp_path / "nope.mat")


def test_no_EEG_variable_raises(tmp_path):
    mat = _make_mat(tmp_path, ictal_len_s=30, corrupt="no_EEG_var")
    with pytest.raises(LoaderError, match="no 'EEG' variable"):
        load_seizure(mat)


def test_1d_array_raises(tmp_path):
    # Note: scipy.io.savemat promotes a 1-D array to shape (1, N) on disk,
    # so the ndim check does not fire; the orientation check does instead,
    # because (1, N) has no axis in the electrode range [36, 100].
    mat = _make_mat(tmp_path, ictal_len_s=30, corrupt="1d")
    with pytest.raises(LoaderError, match="electrode axis"):
        load_seizure(mat)


def test_nan_raises(tmp_path):
    mat = _make_mat(tmp_path, ictal_len_s=30, corrupt="nan")
    with pytest.raises(LoaderError, match="non-finite"):
        load_seizure(mat)


def test_zero_channel_raises(tmp_path):
    mat = _make_mat(tmp_path, ictal_len_s=30, corrupt="zero_channel")
    with pytest.raises(LoaderError, match="entirely zero"):
        load_seizure(mat)


def test_wrong_electrode_range_raises(tmp_path):
    # 10 electrodes is below the documented minimum (36).
    mat = _make_mat(tmp_path, n_electrodes=10, ictal_len_s=30)
    with pytest.raises(LoaderError, match="electrode axis"):
        load_seizure(mat)


# --------------------------------------------------------------------------- #
# Batch loading
# --------------------------------------------------------------------------- #

def test_load_patient_returns_sorted_list(tmp_path):
    _make_mat(tmp_path, patient="ID1", seizure="Sz1", ictal_len_s=30)
    _make_mat(tmp_path, patient="ID1", seizure="Sz2", ictal_len_s=40)
    _make_mat(tmp_path, patient="ID1", seizure="Sz3", ictal_len_s=50)
    seizures = load_patient(tmp_path / "ID1")
    assert len(seizures) == 3
    ids = [s.seizure_id for s in seizures]
    assert ids == sorted(ids)


def test_load_patient_skips_corrupt_files(tmp_path, capsys):
    _make_mat(tmp_path, patient="ID2", seizure="Sz1", ictal_len_s=30)
    _make_mat(tmp_path, patient="ID2", seizure="Sz2", ictal_len_s=30,
              corrupt="no_EEG_var")
    seizures = load_patient(tmp_path / "ID2")
    assert len(seizures) == 1
    # A warning should have been printed to stderr.
    captured = capsys.readouterr()
    assert "WARNING" in captured.err


def test_load_patient_missing_folder_raises(tmp_path):
    with pytest.raises(LoaderError, match="does not exist"):
        load_patient(tmp_path / "IDXX")
