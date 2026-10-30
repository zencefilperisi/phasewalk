"""
Unit tests for loader_bids -- the pure logic (id normalization, SOZ parsing,
onset/window rule, SOZ matching). These need no data and no mne; they encode the
real cases the ds003029 audits (steps 1-4b) surfaced, so a regression in the
parser is caught immediately.

Run:  py -m pytest tests/test_loader_bids.py -v
Requires: pandas, pytest  (mne/numpy not needed for these tests)

loader_bids.py lives in data/; this adds the repo root to sys.path so the
`from data.loader_bids import ...` works the same way as test_connectivity.py.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import pytest

from data.loader_bids import (
    normalize_id, expand_contacts, parse_seizure_window, match_soz,
)


def _events(rows):
    """rows: list of (onset_seconds, trial_type)."""
    return pd.DataFrame(
        [{"onset": t, "duration": 0.0, "trial_type": txt, "value": 0, "sample": 0}
         for t, txt in rows])


# --------------------------------------------------------------------------- #
# normalize_id (ADR 0011 decision 2)
# --------------------------------------------------------------------------- #
def test_normalize_id_strips_underscore():
    assert normalize_id("ummc_001") == "ummc001"
    assert normalize_id("ummc001") == "ummc001"

def test_normalize_id_case_and_punct():
    assert normalize_id("UMMC001") == "ummc001"
    assert normalize_id("pt1") == "pt1"
    assert normalize_id("jh101") == "jh101"


# --------------------------------------------------------------------------- #
# expand_contacts (ADR 0011 decision 1)
# --------------------------------------------------------------------------- #
def test_expand_semicolon_comma_ranges():
    # jh101's actual SOZ string
    assert expand_contacts(" LAT1-2, LAT6-7; LAH6") == \
        ["LAT1", "LAT2", "LAT6", "LAT7", "LAH6"]

def test_expand_colon_separator():
    # pt2's actual SOZ string -- the ':' that broke the step-1 parser
    assert expand_contacts("PST1-4: AST1-2, MST1-2") == \
        ["PST1", "PST2", "PST3", "PST4", "AST1", "AST2", "MST1", "MST2"]

def test_expand_keeps_case_and_singletons():
    assert expand_contacts("AnT1, AnT2, att1-3") == \
        ["AnT1", "AnT2", "att1", "att2", "att3"]

def test_expand_empty():
    assert expand_contacts(None) == []
    assert expand_contacts(float("nan")) == []
    assert expand_contacts("  ") == []


# --------------------------------------------------------------------------- #
# parse_seizure_window (ADR 0012)
# --------------------------------------------------------------------------- #
def test_window_short_seizure_capped_at_offset():
    # jh101/run-01: EEG-SZ onset, ELECTROGRAPHIC END 20.9 s later -> short, capped
    ev = _events([(0.0, "Segment: REC START BP ALL EEG"), (0.0, "A1+A2 OFF"),
                  (58.95, "SZ EVENT # (EEG SZ)"), (59.0, "Z SWC LAT1-2, 6-7, L"),
                  (79.875, "Z ELECTROGRAPHIC END"), (82.8, "Z POST ICTAL VOLT DE")])
    w = parse_seizure_window(ev, window=30.0)
    assert w is not None
    assert w.onset == pytest.approx(58.95)
    assert w.offset == pytest.approx(79.875)
    assert w.short is True
    assert w.win_end == pytest.approx(79.875)

def test_window_no_offset_uses_full_30s():
    # jh104/run-01: push-button onset, no end marker -> full 30 s, not short
    ev = _events([(0.0, "Segment: REC START REC EEG"), (54.855, "Schedule"),
                  (58.995, "SZ EVENT # (PB SZ)"), (61.755, "Z SEMIRHY SLOW PBT4"),
                  (142.830, "Z DEVOLUTION")])
    w = parse_seizure_window(ev, window=30.0)
    assert w.onset == pytest.approx(58.995)
    assert w.offset is None
    assert w.short is False
    assert w.win_end == pytest.approx(88.995)

def test_window_offset_keyword_seizure_over():
    # jh102/run-03: earliest onset is SZ EVENT; 'SEIZURE OVER' is an offset, not
    # an onset, and it is far beyond onset+30 -> full 30 s window
    ev = _events([(0.0, "Segment: REC START BP ALL EEG"),
                  (59.774, "SZ EVENT # (PB SZ)"),
                  (65.949, "Z ELECTROGRAPHIC ONS"),
                  (179.816, "Z SEIZURE OVER"),
                  (180.804, "Z 2 MINUTE SEIZURE")])
    w = parse_seizure_window(ev, window=30.0)
    assert w.onset == pytest.approx(59.774)
    assert w.offset == pytest.approx(179.816)
    assert w.short is False
    assert w.win_end == pytest.approx(89.774)

def test_window_bare_start_end_not_rec_start():
    # pt11/run-03: only 'start'/'end'; must NOT anchor on 'REC START' at t=0
    ev = _events([(0.0, "Segment: REC START Harmon Ref"), (0.0, "A1+A2 OFF"),
                  (67.620, "start"), (167.670, "end")])
    w = parse_seizure_window(ev, window=30.0)
    assert w.onset == pytest.approx(67.620)
    assert w.offset == pytest.approx(167.670)
    assert w.short is False

def test_window_no_onset_returns_none():
    # pt17/run-01: only administrative markers -> run must be dropped
    ev = _events([(0.0, "Segment: REC START EL-Ref EEG"), (0.0, "A1+A2 OFF")])
    assert parse_seizure_window(ev, window=30.0) is None

def test_window_very_short_pt13():
    # pt13/run-01: ~9 s seizure -> capped and flagged short
    ev = _events([(0.0, "Segment: REC START"), (68.08, "sz onset #18"),
                  (77.28, "offset")])
    w = parse_seizure_window(ev, window=30.0)
    assert w.duration == pytest.approx(9.2, abs=0.01)
    assert w.short is True

def test_window_ummc_clean_long():
    ev = _events([(0.0, "Segment: REC START"), (60.0, "sz onset"), (111.0, "sz offset")])
    w = parse_seizure_window(ev, window=30.0)
    assert w.onset == pytest.approx(60.0)
    assert w.short is False
    assert w.win_end == pytest.approx(90.0)


# --------------------------------------------------------------------------- #
# match_soz (ADR 0011 decisions 1, 5)
# --------------------------------------------------------------------------- #
def test_match_soz_good_bad_absent():
    good = ["LAT1", "LAT2", "LAH6"]
    allc = ["LAT1", "LAT2", "LAH6", "LAT6", "LAT7"]   # LAT6/LAT7 present but bad
    m = match_soz("LAT1-2, LAT6-7; LAH6", good, allc)
    assert m.matched == ["LAT1", "LAT2", "LAH6"]
    assert m.in_bad == ["LAT6", "LAT7"]
    assert m.unmatched == []

def test_match_soz_absent_contact():
    m = match_soz("att1-3", ["ANT1", "MEST1"], ["ANT1", "MEST1"])
    assert m.matched == []
    assert m.in_bad == []
    assert m.unmatched == ["att1", "att2", "att3"]

def test_match_soz_case_insensitive():
    m = match_soz("lat1", ["LAT1"], ["LAT1"])
    assert m.matched == ["lat1"]