"""
ds003029 - Step 4b: full event dump for the onset-rule stragglers.

The step-4 audit left two things to resolve before fixing the ictal-window /
onset rule:
  (1) 3 runs had NO detectable onset (pt11/run-03, pt17/run-01, pt17/run-03) --
      are they truly unmarked, or marked with different wording?
  (2) JHH uses "(EEG SZ)" = electrographic onset and "(PB SZ)" = push-button
      (patient/clinical) onset. For SOZ localization we want the electrographic
      onset. Do PB-only runs also carry an earlier electrographic marker we
      should prefer?

This prints the FULL event list for those runs so we can write a robust rule.
Read-only. Requires: py -m pip install pandas
"""
from pathlib import Path
import pandas as pd

DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")

# (subject, run) pairs to dump
TARGETS = [
    ("pt11", "run-03"),   # NO-ONSET
    ("pt17", "run-01"),   # NO-ONSET
    ("pt17", "run-03"),   # NO-ONSET
    ("jh104", "run-01"),  # JHH, PB-only in audit -> look for electrographic marker
    ("jh102", "run-03"),  # JHH, PB-only, NO-OFFSET -> inspect end markers too
]


def find_events(sub, run):
    hits = list(DATA_ROOT.glob(f"sub-{sub}/**/ieeg/*task-ictal*{run}*_events.tsv"))
    return hits[0] if hits else None


for sub, run in TARGETS:
    f = find_events(sub, run)
    print("=" * 70)
    print(f"{sub} / {run}")
    print("=" * 70)
    if f is None:
        print("   (no events file found)")
        continue
    ev = pd.read_csv(f, sep="\t")
    for _, r in ev.iterrows():
        print(f"   t={r['onset']:>9.3f}s  value={str(r['value']):>4}  {r['trial_type']}")
    print()