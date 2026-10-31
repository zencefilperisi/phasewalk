"""
ds003029 - Step 3: resolve the three straggler issues the Step-2 audit found,
so we can finalise the SOZ parser + inclusion policy before writing the loader.

(a) ummc002: SOZ tokens att1/att2/att3 matched NOTHING (not good, not bad).
    -> print ummc002's channel names so we can see what they should map to.
(b) pt01 vs pt1: pt1 has a clinical row, pt01 does not. Duplicate folder, or a
    distinct subject? -> compare their ictal runs and channel sets.
(c) 6 subjects have SOZ contacts marked status=bad. -> list exactly which SOZ
    contacts are bad and WHY (status_description), to decide include vs exclude.

Read-only (channels.tsv + clinical xlsx). No signal, no CTQW.
Requires:  py -m pip install pandas openpyxl
"""
import re
from pathlib import Path

import pandas as pd

DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")
CLINICAL_XLSX = DATA_ROOT / "sourcedata" / "clinical_data_summary.xlsx"


def expand_contacts(s):
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return []
    out = []
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


def ictal_tsvs(label):
    return sorted((DATA_ROOT / f"sub-{label}").glob("**/ieeg/*task-ictal*_channels.tsv"))


def soz_raw_for(clin, label):
    r = clin[clin["dataset_id"].astype(str).str.strip().str.lower() == label.lower()]
    return None if r.empty else r.iloc[0].get("soz_contacts")


clin = pd.read_excel(CLINICAL_XLSX)

# --------------------------------------------------------------------------- #
# (a) ummc002 channel names
# --------------------------------------------------------------------------- #
print("=" * 70)
print("(a) ummc002 -- what should att1/att2/att3 map to?")
print("=" * 70)
t = ictal_tsvs("ummc002")[0]
ch = pd.read_csv(t, sep="\t")
good = ch[(ch["type"].isin(["ECOG", "SEEG"])) & (ch["status"] == "good")]["name"].tolist()
print(f"SOZ raw: {soz_raw_for(clin, 'ummc002')!r}")
print(f"parsed : {expand_contacts(soz_raw_for(clin, 'ummc002'))}")
print(f"\nall {len(good)} good channel names:")
print(good)
print(f"\nchannels containing 'T' (candidate matches for att*):")
print([n for n in ch['name'].astype(str) if 'T' in n.upper()])

# --------------------------------------------------------------------------- #
# (b) pt01 vs pt1
# --------------------------------------------------------------------------- #
print("\n" + "=" * 70)
print("(b) pt01 vs pt1 -- duplicate or distinct?")
print("=" * 70)
for lab in ["pt1", "pt01"]:
    tsvs = ictal_tsvs(lab)
    print(f"\nsub-{lab}: {len(tsvs)} ictal runs")
    names_sets = []
    for tt in tsvs:
        c = pd.read_csv(tt, sep="\t")
        g = c[(c["type"].isin(["ECOG", "SEEG"])) & (c["status"] == "good")]["name"].tolist()
        names_sets.append(set(g))
        print(f"   {tt.name}  ({len(c)} ch, {len(g)} good)")
    if names_sets:
        print(f"   good channel example: {sorted(list(names_sets[0]))[:12]}")
# are pt1 and pt01 the same recordings?
try:
    g1 = set(pd.read_csv(ictal_tsvs('pt1')[0], sep='\t')['name'].astype(str))
    g01 = set(pd.read_csv(ictal_tsvs('pt01')[0], sep='\t')['name'].astype(str))
    inter = len(g1 & g01)
    print(f"\npt1 run-01 vs pt01 run-01: {inter} shared channel names "
          f"(pt1={len(g1)}, pt01={len(g01)}) -> "
          f"{'LIKELY SAME MONTAGE' if inter > 0.8*min(len(g1),len(g01)) else 'different'}")
except Exception as e:
    print("compare failed:", e)
print(f"\nany xlsx dataset_id containing '01': "
      f"{[x for x in clin['dataset_id'].astype(str).unique() if '01' in x]}")

# --------------------------------------------------------------------------- #
# (c) SOZ-in-bad conflicts: which contacts, and why bad?
# --------------------------------------------------------------------------- #
print("\n" + "=" * 70)
print("(c) SOZ contacts marked status=bad -- include or exclude?")
print("=" * 70)
for lab in ["jh103", "pt15", "pt8", "umf004", "ummc003", "ummc008"]:
    t = ictal_tsvs(lab)[0]
    ch = pd.read_csv(t, sep="\t")
    byname = {str(r["name"]).upper(): r for _, r in ch.iterrows()}
    soz = expand_contacts(soz_raw_for(clin, lab))
    print(f"\nsub-{lab}  (first ictal run: {t.name})")
    for c in soz:
        r = byname.get(c.upper())
        if r is None:
            continue
        if r["status"] == "bad":
            print(f"   BAD  {c:<8} type={r['type']:<5} desc={r['status_description']!r}")

print("\n" + "=" * 70)
print("Done. Paste the full output back.")
print("=" * 70)