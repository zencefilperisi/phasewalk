"""
ds003029 — Step 1: end-to-end data-integrity sanity check on ONE subject.

Before building the full BIDS loader + clinical-validation pipeline on all 37
subjects, prove on a single subject that we can correctly:
  (1) read the BrainVision iEEG signal (MNE),
  (2) select the good iEEG channels from channels.tsv and align them to the
      signal,
  (3) locate the electrographic seizure onset/offset from events.tsv,
  (4) parse this subject's clinician SOZ contacts from the clinical summary
      xlsx, and
  (5) VERIFY the parsed SOZ contact names actually exist among the recording's
      good channels -- the critical name-matching step. Mismatches here would
      silently corrupt any clinical result, so we check it first.

Nothing is computed (no FCM / CTQW). This is a data-integrity check only.
Run it, paste the full output back, and we refine the parsers before scaling.

Requires:  py -m pip install mne pandas openpyxl numpy
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import mne

# --------------------------------------------------------------------------- #
# Config  (switch SUBJECT to "sub-pt1" if jh101 has no SOZ row -- see note)
# --------------------------------------------------------------------------- #
DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")
SUBJECT = "sub-jh101"
SESSION = "ses-presurgery"
TASK = "ictal"
RUN = "run-01"

mne.set_log_level("ERROR")

ieeg_dir = DATA_ROOT / SUBJECT / SESSION / "ieeg"
stem = f"{SUBJECT}_{SESSION}_task-{TASK}_acq-ecog_{RUN}"
vhdr = ieeg_dir / f"{stem}_ieeg.vhdr"
chan_tsv = ieeg_dir / f"{stem}_channels.tsv"
events_tsv = ieeg_dir / f"{stem}_events.tsv"
clinical_xlsx = DATA_ROOT / "sourcedata" / "clinical_data_summary.xlsx"

print("=" * 70)
print(f"ds003029 sanity check -- {SUBJECT} / {TASK} / {RUN}")
print("=" * 70)

# --------------------------------------------------------------------------- #
# (1) Read the BrainVision signal
# --------------------------------------------------------------------------- #
print("\n[1] Reading BrainVision signal ...")
raw = mne.io.read_raw_brainvision(str(vhdr), preload=True, verbose="ERROR")
data = raw.get_data()                       # (n_channels, n_samples)
sfreq = raw.info["sfreq"]
print(f"    signal shape : {data.shape}  (channels x samples)")
print(f"    sfreq        : {sfreq} Hz")
print(f"    duration     : {data.shape[1] / sfreq:.1f} s")
print(f"    n ch (raw)   : {len(raw.ch_names)}")

# --------------------------------------------------------------------------- #
# (2) Channel selection from channels.tsv
# --------------------------------------------------------------------------- #
print("\n[2] Channel selection from channels.tsv ...")
ch = pd.read_csv(chan_tsv, sep="\t")
print(f"    channels.tsv rows : {len(ch)}")
print(f"    types present     : {dict(ch['type'].value_counts())}")
print(f"    status counts     : {dict(ch['status'].value_counts())}")

is_ieeg = ch["type"].isin(["ECOG", "SEEG"])
is_good = ch["status"] == "good"
good = ch[is_ieeg & is_good].copy()
good_names = good["name"].tolist()
print(f"    good iEEG channels: {len(good_names)}")

# alignment: channels.tsv row order should match the signal's channel order
align_ok = (len(ch) == len(raw.ch_names)) and all(
    a == b for a, b in zip(ch["name"].tolist(), raw.ch_names)
)
print(f"    channels.tsv order matches raw signal order? {align_ok}")
missing_in_raw = [n for n in good_names if n not in raw.ch_names]
print(f"    good names NOT found in raw signal: {missing_in_raw if missing_in_raw else 'none'}")

# --------------------------------------------------------------------------- #
# (3) Seizure onset/offset from events.tsv
# --------------------------------------------------------------------------- #
print("\n[3] Events (electrographic onset/offset) ...")
ev = pd.read_csv(events_tsv, sep="\t")
print("    full events table:")
for _, r in ev.iterrows():
    print(f"      t={r['onset']:>9.3f}s  value={str(r['value']):>4}  {r['trial_type']}")

# tentative detection -- to be refined into a universal rule across subjects
tt = ev["trial_type"].astype(str)
onset_mask = tt.str.contains(r"\bSZ\b|EEG SZ|ONSET|LOW VOLTAGE FAST", case=False, regex=True)
offset_mask = tt.str.contains(r"ELECTROGRAPHIC END|OFFSET|POST ICTAL", case=False, regex=True)
onset_t = ev.loc[onset_mask, "onset"].min() if onset_mask.any() else None
offset_t = ev.loc[offset_mask, "onset"].min() if offset_mask.any() else None
print(f"\n    TENTATIVE electrographic onset : {onset_t} s")
print(f"    TENTATIVE electrographic offset: {offset_t} s")
if onset_t is not None and offset_t is not None:
    print(f"    => tentative ictal duration    : {offset_t - onset_t:.1f} s")

# --------------------------------------------------------------------------- #
# (4) SOZ contacts from the clinical summary xlsx
# --------------------------------------------------------------------------- #
print("\n[4] SOZ contacts from clinical_data_summary.xlsx ...")
clin = pd.read_excel(clinical_xlsx)
label = SUBJECT.replace("sub-", "")
clin["_did"] = clin["dataset_id"].astype(str).str.strip().str.lower()
row = clin[clin["_did"] == label.lower()]
soz_raw = None
if row.empty:
    print(f"    NO xlsx row for dataset_id == '{label}'.")
    print(f"    available dataset_ids: {sorted(clin['dataset_id'].astype(str).unique())}")
else:
    r = row.iloc[0]
    soz_raw = r.get("soz_contacts")
    print(f"    matched row    : dataset_id={r['dataset_id']}  center={r.get('clinical_center')}")
    print(f"    outcome        : {r.get('outcome')}  engel={r.get('engel_score')}  ilae={r.get('ilae_score')}")
    print(f"    surgery_type   : {r.get('surgery_type')}")
    print(f"    ez_hypo_region : {r.get('ez_hypo_brainregion')}")
    print(f"    soz_contacts   : {soz_raw!r}")
    print(f"    bad_contacts   : {r.get('bad_contacts')!r}")


def expand_contacts(s):
    """First-pass parser: 'PD1-4; AD1-4; ATT1-2' -> [PD1,PD2,PD3,PD4,AD1,...].
    Split on ; and , ; expand <prefix><start>-<end> ranges. Keep original case.
    Deliberately simple so parse failures stay visible in step (5)."""
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return []
    out = []
    for tok in re.split(r"[;,]", str(s)):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"^([A-Za-z]+?)\s*(\d+)\s*-\s*(\d+)$", tok)
        if m:
            pre, a, b = m.group(1), int(m.group(2)), int(m.group(3))
            for k in range(a, b + 1):
                out.append(f"{pre}{k}")
        else:
            out.append(tok)
    return out


soz_list = expand_contacts(soz_raw)
print(f"    parsed SOZ contacts ({len(soz_list)}): {soz_list}")

# --------------------------------------------------------------------------- #
# (5) CRITICAL: do parsed SOZ names match the recording's good channels?
# --------------------------------------------------------------------------- #
print("\n[5] SOZ <-> channel name matching (the critical check) ...")
good_upper = {n.upper(): n for n in good_names}
matched = [n for n in soz_list if n.upper() in good_upper]
unmatched = [n for n in soz_list if n.upper() not in good_upper]
print(f"    SOZ contacts matched to good channels: {len(matched)}/{len(soz_list)}")
print(f"      matched   : {matched}")
print(f"      UNMATCHED : {unmatched}")
if soz_list:
    print(f"    match rate: {len(matched) / len(soz_list):.0%}")
print("\n    (unmatched = parse error or cross-naming mismatch -> must fix "
      "before trusting any clinical result)")

print("\n" + "=" * 70)
print("Sanity check done. Paste the full output back.")
print("=" * 70)