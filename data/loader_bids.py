"""
ds003029 BIDS loader for the CTQW clinical-validation pipeline (Phase 3).

Implements the policy fixed in:
  - ADR 0011 (channel / SOZ inclusion): good ECOG/SEEG channels only; SOZ parsed
    from the clinical xlsx, matched case-insensitively to good channels; id
    normalization; pt01 excluded as a pt1 duplicate.
  - ADR 0012 (ictal window): anchor on the earliest seizure-onset marker; window
    [onset, onset+30s]; offset caps only short seizures; drop unmarked runs.

It produces, per ictal run, the good-channel ictal-window signal plus the
subject's SOZ reference, ready for the frozen FCM -> Layer P -> CTQW -> PR
pipeline (unchanged from Phase 2C / SWEC).

Public API:
    normalize_id(s)
    expand_contacts(s)                 -> list[str]
    parse_seizure_window(events_df, window=30.0) -> SeizureWindow | None
    load_clinical_table(data_root)     -> dict[norm_id -> dict]
    match_soz(soz_raw, good_names, all_names) -> SozMatch
    discover_ictal_runs(subject_dir)   -> list[RunPaths]
    load_ictal_window(run_paths, window=30.0) -> IctalRun | None

Requires: mne, pandas, numpy, openpyxl
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# Excluded / duplicate subjects per ADR 0011 decision 3.
EXCLUDE_SUBJECTS = {"pt01"}
WINDOW_SECONDS = 30.0

# --------------------------------------------------------------------------- #
# Marker patterns (ADR 0012). OFFSET takes precedence over ONSET so that
# "SEIZURE OVER" / "SZ OFFSET" are never mistaken for onsets.
# --------------------------------------------------------------------------- #
_ADMIN = re.compile(r"segment|rec start|a1\+a2|schedule|mark o", re.I)
_OFFSET = re.compile(
    r"offset|electrographic end|\bend\b|termination|post[-\s]?ictal|"
    r"seizure over|seizure off|\bsz end\b|\bover\b", re.I)
_ONSET = re.compile(
    r"onset|eeg sz start|sz event|ictal onset|definite onset|seizure|\bsz\b|\bstart\b",
    re.I)


# --------------------------------------------------------------------------- #
# Pure helpers (unit-testable without data)
# --------------------------------------------------------------------------- #
def normalize_id(s: str) -> str:
    """Lower-case and strip non-alphanumerics (ADR 0011 decision 2).
    'ummc_001' and 'ummc001' both -> 'ummc001'."""
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def expand_contacts(s) -> list[str]:
    """Parse a SOZ/bad contact string into individual contact names (ADR 0011).
    Separators ; , : ; ranges <prefix><a>-<b> expanded; original case kept.
    Anything unparsable is kept verbatim so matching flags it rather than
    silently dropping it."""
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return []
    out: list[str] = []
    for tok in re.split(r"[;,:]", str(s)):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"^([A-Za-z]+?)\s*(\d+)\s*-\s*(\d+)$", tok)
        if m:
            pre, a, b = m.group(1), int(m.group(2)), int(m.group(3))
            if a <= b and (b - a) < 200:
                out.extend(f"{pre}{k}" for k in range(a, b + 1))
            else:
                out.append(tok)
        else:
            out.append(tok)
    return out


def _is_timestamp_row(text: str) -> bool:
    return bool(re.match(r"^\s*\+?\d", str(text)))


@dataclass
class SeizureWindow:
    onset: float
    offset: float | None
    win_start: float
    win_end: float
    duration: float
    short: bool
    onset_text: str


def parse_seizure_window(events: pd.DataFrame, window: float = WINDOW_SECONDS
                         ) -> SeizureWindow | None:
    """Apply ADR 0012: earliest onset marker as anchor; window [onset,
    onset+window]; offset caps the window only if the seizure is shorter.
    Returns None if no onset marker exists (run should be dropped)."""
    tt = events["trial_type"].astype(str)
    on_t = on_txt = None
    off_t = None
    for _, r in events.iterrows():
        txt = str(r["trial_type"])
        if _ADMIN.search(txt) or _is_timestamp_row(txt):
            continue
        is_off = bool(_OFFSET.search(txt))
        is_on = bool(_ONSET.search(txt)) and not is_off
        t = float(r["onset"])
        if is_on and (on_t is None or t < on_t):
            on_t, on_txt = t, txt
    if on_t is None:
        return None
    # earliest offset at/after onset
    for _, r in events.iterrows():
        txt = str(r["trial_type"])
        if _ADMIN.search(txt) or _is_timestamp_row(txt):
            continue
        if _OFFSET.search(txt):
            t = float(r["onset"])
            if t >= on_t and (off_t is None or t < off_t):
                off_t = t
    win_start = on_t
    win_end = on_t + window
    short = False
    if off_t is not None and off_t < win_end:
        win_end = off_t
        short = True
    return SeizureWindow(on_t, off_t, win_start, win_end,
                         win_end - win_start, short, on_txt or "")


@dataclass
class SozMatch:
    matched: list[str] = field(default_factory=list)   # SOZ ∩ good channels
    in_bad: list[str] = field(default_factory=list)    # SOZ present but status=bad
    unmatched: list[str] = field(default_factory=list)  # SOZ not in montage


def match_soz(soz_raw, good_names, all_names) -> SozMatch:
    """Match parsed SOZ contacts to channel names (ADR 0011 decisions 1, 5).
    `good_names`: good iEEG channels (analysis set). `all_names`: every channel
    in the recording (to distinguish bad-but-present from truly absent)."""
    good_up = {n.upper() for n in good_names}
    all_up = {n.upper() for n in all_names}
    m = SozMatch()
    for c in expand_contacts(soz_raw):
        u = c.upper()
        if u in good_up:
            m.matched.append(c)
        elif u in all_up:
            m.in_bad.append(c)
        else:
            m.unmatched.append(c)
    return m


# --------------------------------------------------------------------------- #
# Clinical table
# --------------------------------------------------------------------------- #
def load_clinical_table(data_root) -> dict[str, dict]:
    """Read sourcedata/clinical_data_summary.xlsx, keyed by normalized id."""
    xlsx = Path(data_root) / "sourcedata" / "clinical_data_summary.xlsx"
    df = pd.read_excel(xlsx)
    out: dict[str, dict] = {}
    for _, r in df.iterrows():
        out[normalize_id(r["dataset_id"])] = r.to_dict()
    return out


# --------------------------------------------------------------------------- #
# Data access
# --------------------------------------------------------------------------- #
@dataclass
class RunPaths:
    subject: str
    run: str
    vhdr: Path
    channels_tsv: Path
    events_tsv: Path


def discover_ictal_runs(subject_dir) -> list[RunPaths]:
    """All ictal runs for a subject that have the FULL fileset: the BrainVision
    header (.vhdr), the data file (.eeg), channels.tsv and events.tsv. Runs
    whose .eeg data file is missing (an incomplete download) are skipped, so the
    pipeline runs on whatever is fully present."""
    subject_dir = Path(subject_dir)
    subj = subject_dir.name.replace("sub-", "")
    runs: list[RunPaths] = []
    for vhdr in sorted(subject_dir.glob("**/ieeg/*task-ictal*_ieeg.vhdr")):
        stem = vhdr.name.replace("_ieeg.vhdr", "")
        ch = vhdr.with_name(stem + "_channels.tsv")
        ev = vhdr.with_name(stem + "_events.tsv")
        eeg = vhdr.with_name(stem + "_ieeg.eeg")
        run = next((p for p in stem.split("_") if p.startswith("run-")), "run-?")
        if ch.exists() and ev.exists() and eeg.exists():
            runs.append(RunPaths(subj, run, vhdr, ch, ev))
    return runs


@dataclass
class IctalRun:
    subject: str
    run: str
    sfreq: float
    good_names: list[str]        # good iEEG channel names, in signal order
    all_names: list[str]         # every channel name (for SOZ bad/absent split)
    signal: np.ndarray           # (n_good_channels, n_window_samples)
    window: SeizureWindow
    dropped_reason: str | None = None


def _good_and_all(channels_tsv) -> tuple[list[str], list[str]]:
    c = pd.read_csv(channels_tsv, sep="\t")
    alln = c["name"].astype(str).tolist()
    good = c[(c["type"].isin(["ECOG", "SEEG"])) & (c["status"] == "good")]
    return good["name"].astype(str).tolist(), alln


def load_ictal_window(run: RunPaths, window: float = WINDOW_SECONDS
                      ) -> IctalRun | None:
    """Load one ictal run's good-channel signal over the ADR-0012 window.
    Returns None (with reason logged on the object) if the run must be dropped."""
    import mne
    mne.set_log_level("ERROR")

    events = pd.read_csv(run.events_tsv, sep="\t")
    win = parse_seizure_window(events, window)
    if win is None:
        return IctalRun(run.subject, run.run, 0.0, [], [],
                        np.empty((0, 0)), SeizureWindow(0, None, 0, 0, 0, False, ""),
                        dropped_reason="no onset marker")

    good_names, all_names = _good_and_all(run.channels_tsv)
    try:
        raw = mne.io.read_raw_brainvision(str(run.vhdr), preload=True, verbose="ERROR")
    except (FileNotFoundError, OSError) as e:
        return IctalRun(run.subject, run.run, 0.0, [], [], np.empty((0, 0)),
                        win, dropped_reason=f"signal unreadable ({e.__class__.__name__})")
    sfreq = float(raw.info["sfreq"])
    rec_dur = raw.n_times / sfreq
    if win.win_end > rec_dur:          # not enough signal after onset
        win.win_end = min(win.win_end, rec_dur)
        win.duration = win.win_end - win.win_start
    if win.duration < 1.0:
        return IctalRun(run.subject, run.run, sfreq, [], [],
                        np.empty((0, 0)), win, dropped_reason="window too short")

    # select good channels by name, preserving signal order
    name_to_idx = {n: i for i, n in enumerate(raw.ch_names)}
    good_idx = [name_to_idx[n] for n in good_names if n in name_to_idx]
    good_names = [raw.ch_names[i] for i in good_idx]
    data = raw.get_data()[good_idx]
    s0 = int(round(win.win_start * sfreq))
    s1 = int(round(win.win_end * sfreq))
    sig = data[:, s0:s1]
    return IctalRun(run.subject, run.run, sfreq, good_names, all_names, sig, win)


if __name__ == "__main__":
    import sys
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
        r"C:\Users\User\Desktop\ds003029")
    clinical = load_clinical_table(root)
    print(f"clinical rows: {len(clinical)}")
    subs = sorted(p for p in root.glob("sub-*")
                  if p.is_dir() and p.name.replace("sub-", "") not in EXCLUDE_SUBJECTS)
    print(f"subjects (excl. duplicates): {len(subs)}")