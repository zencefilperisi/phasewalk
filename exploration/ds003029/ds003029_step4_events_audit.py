"""
ds003029 - Step 4: electrographic onset/offset audit across ALL ictal runs.

Step 1 found a clean onset/offset for jh101, but the README warns seizure
markers are NOT consistently named across centers. Before we can define the
ictal segment, we must see every marker vocabulary and check that a plausible
onset AND offset can be found for every ictal run -- the timing analogue of the
step-2 SOZ audit.

Read-only (events.tsv). No signal, no CTQW.

For each ictal run we attempt a TENTATIVE onset (earliest onset-like marker) and
offset (earliest offset-like marker at/after onset), then flag runs where that
fails (no onset, no offset, or onset >= offset). We also dump the full distinct
vocabulary of seizure-relevant marker texts so we can design the final rule.

Requires:  py -m pip install pandas
"""
import re
from pathlib import Path
from collections import Counter

import pandas as pd

DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")

# Broad, tentative patterns -- deliberately permissive so we SEE what matches;
# the final rule gets fixed in an ADR from this audit's output.
ONSET_PAT = re.compile(
    r"ONSET|\bSZ\b|EEG SZ|SEIZURE|LOW VOLTAGE FAST|\bUEO\b|CLINICAL ONSET|SZ EVENT",
    re.IGNORECASE)
OFFSET_PAT = re.compile(
    r"OFFSET|ELECTROGRAPHIC END|\bEND\b|TERMINATION|POST[- ]?ICTAL|SZ END|SEIZURE END",
    re.IGNORECASE)
# marker vocabulary we care about (for the dump)
SZ_KEYWORDS = re.compile(
    r"ONSET|OFFSET|\bSZ\b|SEIZURE|ICTAL|\bEND\b|LOW VOLTAGE|\bUEO\b|TERMINATION|SWC|SPIK",
    re.IGNORECASE)


def ictal_event_files(root):
    return sorted(root.glob("sub-*/**/ieeg/*task-ictal*_events.tsv"))


files = ictal_event_files(DATA_ROOT)
print(f"found {len(files)} ictal event files\n")

rows = []
vocab = Counter()
for f in files:
    # subject/run short id
    sub = f.name.split("_")[0].replace("sub-", "")
    run = next((p for p in f.name.split("_") if p.startswith("run-")), "run-?")
    try:
        ev = pd.read_csv(f, sep="\t")
    except Exception as e:
        rows.append((sub, run, None, None, None, None, f"READ-ERROR: {e}"))
        continue
    tt = ev["trial_type"].astype(str)

    # collect vocabulary
    for s in tt[tt.str.contains(SZ_KEYWORDS)]:
        vocab[s.strip()] += 1

    on = ev[tt.str.contains(ONSET_PAT)]
    off = ev[tt.str.contains(OFFSET_PAT)]
    onset_t = float(on["onset"].min()) if len(on) else None
    onset_txt = on.loc[on["onset"].idxmin(), "trial_type"] if len(on) else ""
    # offset: earliest offset marker strictly after onset, else earliest offset
    offset_t = offset_txt = None
    if len(off):
        if onset_t is not None:
            after = off[off["onset"] >= onset_t]
            pick = after if len(after) else off
        else:
            pick = off
        offset_t = float(pick["onset"].min())
        offset_txt = pick.loc[pick["onset"].idxmin(), "trial_type"]

    flag = ""
    if onset_t is None:
        flag = "NO-ONSET"
    elif offset_t is None:
        flag = "NO-OFFSET"
    elif offset_t <= onset_t:
        flag = "OFFSET<=ONSET"
    rows.append((sub, run, onset_t, onset_txt, offset_t, offset_txt, flag))

# ---- per-run table -------------------------------------------------------- #
print(f"{'subject':<10}{'run':<8}{'onset':>9}{'offset':>9}{'dur':>8}  flag   onset_txt | offset_txt")
print("-" * 110)
for (sub, run, on_t, on_x, off_t, off_x, flag) in rows:
    dur = f"{off_t - on_t:7.1f}" if (on_t is not None and off_t is not None) else "      -"
    on_s = f"{on_t:9.2f}" if on_t is not None else "        -"
    off_s = f"{off_t:9.2f}" if off_t is not None else "        -"
    snip = f"{str(on_x)[:28]:<28} | {str(off_x)[:28]}"
    print(f"{sub:<10}{run:<8}{on_s}{off_s}{dur}  {flag:<13} {snip}")

# ---- summary -------------------------------------------------------------- #
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
ok = [r for r in rows if r[6] == ""]
print(f"ictal runs total        : {len(rows)}")
print(f"  onset+offset OK       : {len(ok)}")
for flg in ["NO-ONSET", "NO-OFFSET", "OFFSET<=ONSET"]:
    bad = [f"{r[0]}/{r[1]}" for r in rows if r[6] == flg]
    print(f"  {flg:<14}      : {len(bad)}  {bad}")

if ok:
    durs = [r[4] - r[2] for r in ok]
    import statistics
    print(f"\ntentative ictal durations (OK runs): "
          f"min={min(durs):.1f}s  median={statistics.median(durs):.1f}s  max={max(durs):.1f}s")

print("\n--- distinct seizure-relevant marker vocabulary (text : count) ---")
for txt, c in sorted(vocab.items(), key=lambda kv: -kv[1]):
    print(f"  {c:>3}  {txt}")

print("\n" + "=" * 70)
print("Audit done. Paste the full output back.")
print("=" * 70)