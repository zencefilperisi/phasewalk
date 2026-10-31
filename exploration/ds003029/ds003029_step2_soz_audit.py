"""
ds003029 - Step 2: SOZ parse + channel-match audit across ALL subjects.

Step 1 passed cleanly on one easy subject (jh101, 5/5). A clean pass on one
easy case does NOT prove the parser generalises: SOZ formats vary across
centers (e.g. pt2 uses a ':' separator that step 1's parser could not handle).
This step stress-tests the SOZ parser and the SOZ<->channel name matching on
EVERY subject, to surface every clinical-annotation format we must handle
before trusting any clinical result.

Read-only: channels.tsv (good iEEG channel names) + clinical_data_summary.xlsx
(SOZ contacts). No signal is loaded, no FCM/CTQW. Runs in seconds.

For each subject, each parsed SOZ contact is classified as:
    matched : present among good iEEG channels            (usable)
    in_bad  : present in the recording but status=bad      (SOZ electrode
              excluded as bad -> a real conflict to decide on)
    unmatched: not present at all                          (parse error /
              cross-center naming mismatch -> must fix)

Output: a per-subject table + a summary (match-rate, every distinct unmatched
token, Engel-outcome distribution). Paste it back; we fix the parser from the
failures it surfaces, then write the final rule into an ADR.

Requires:  py -m pip install pandas openpyxl
"""
import re
from pathlib import Path
from collections import Counter

import pandas as pd

DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")
CLINICAL_XLSX = DATA_ROOT / "sourcedata" / "clinical_data_summary.xlsx"


def expand_contacts(s):
    """Parse e.g. 'PST1-4: AST1-2, MST1-2' -> [PST1..PST4, AST1, AST2, ...].
    Separators ; , : and ranges <prefix><a>-<b>. Original case preserved.
    Deliberately conservative: anything it cannot parse is kept verbatim so it
    shows up as 'unmatched' rather than being silently dropped."""
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
                out.append(tok)          # implausible range -> keep raw
        else:
            out.append(tok)
    return out


def channels_over_ictal_runs(subject_dir):
    """Union of good iEEG channel names and of ALL channel names across the
    subject's ictal runs. Returns (good_set, all_set, n_ictal_runs)."""
    good, allch = set(), set()
    tsvs = list(subject_dir.glob("**/ieeg/*task-ictal*_channels.tsv"))
    for t in tsvs:
        c = pd.read_csv(t, sep="\t")
        allch.update(c["name"].astype(str))
        m = c["type"].isin(["ECOG", "SEEG"]) & (c["status"] == "good")
        good.update(c.loc[m, "name"].astype(str))
    return good, allch, len(tsvs)


clin = pd.read_excel(CLINICAL_XLSX)
clin["_did"] = clin["dataset_id"].astype(str).str.strip().str.lower()

subjects = sorted(p for p in DATA_ROOT.glob("sub-*") if p.is_dir())
print(f"found {len(subjects)} subject folders\n")

recs = []
all_unmatched = Counter()
for sdir in subjects:
    label = sdir.name.replace("sub-", "")
    good, allch, n_ict = channels_over_ictal_runs(sdir)
    good_up = {c.upper() for c in good}
    all_up = {c.upper() for c in allch}
    r = clin[clin["_did"] == label.lower()]

    rec = dict(label=label, n_ict=n_ict, ngood=len(good), center="-",
               outcome="-", engel="-", nsoz=0, matched=0, in_bad=0,
               unmatched=0, unm_tokens=[], flag="")

    if r.empty:
        rec["flag"] = "NO-XLSX-ROW"
        recs.append(rec)
        continue
    r = r.iloc[0]
    rec["center"] = r.get("clinical_center")
    rec["outcome"] = r.get("outcome")
    rec["engel"] = r.get("engel_score")
    soz = expand_contacts(r.get("soz_contacts"))
    rec["nsoz"] = len(soz)

    if n_ict == 0 or not good:
        rec["flag"] = "NO-ICTAL"          # SOZ not matched (no ictal channels)
        recs.append(rec)
        continue

    mg = [c for c in soz if c.upper() in good_up]
    ib = [c for c in soz if c.upper() not in good_up and c.upper() in all_up]
    un = [c for c in soz if c.upper() not in all_up]
    rec.update(matched=len(mg), in_bad=len(ib), unmatched=len(un), unm_tokens=un)
    for c in un:
        all_unmatched[c] += 1
    recs.append(rec)

# ---- per-subject table ---------------------------------------------------- #
print(f"{'subject':<11}{'ict':>4}{'good':>5}{'center':>7}{'out':>5}{'eng':>5}"
      f"{'nSOZ':>5}{'mtch':>5}{'bad':>4}{'unm':>4}  flag / unmatched")
print("-" * 100)
for x in recs:
    extra = x["flag"] if x["flag"] else (str(x["unm_tokens"]) if x["unm_tokens"] else "")
    print(f"{x['label']:<11}{x['n_ict']:>4}{x['ngood']:>5}{str(x['center']):>7}"
          f"{str(x['outcome']):>5}{str(x['engel']):>5}{x['nsoz']:>5}{x['matched']:>5}"
          f"{x['in_bad']:>4}{x['unmatched']:>4}  {extra}")

# ---- summary -------------------------------------------------------------- #
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
usable = [x for x in recs if not x["flag"] and x["nsoz"] > 0]
print(f"subjects total                     : {len(recs)}")
print(f"  no xlsx row                      : {sum(1 for x in recs if x['flag']=='NO-XLSX-ROW')}")
print(f"  no ictal runs                    : {sum(1 for x in recs if x['flag']=='NO-ICTAL')}")
print(f"  xlsx row but empty SOZ           : {sum(1 for x in recs if not x['flag'] and x['nsoz']==0)}")
print(f"  usable (ictal + parsed SOZ)      : {len(usable)}")
full = [x for x in usable if x["unmatched"] == 0 and x["in_bad"] == 0]
print(f"    SOZ 100% matched to good ch    : {len(full)}")
partial = [x for x in usable if x["unmatched"] > 0]
print(f"    with >=1 UNMATCHED SOZ token   : {len(partial)}  -> {[x['label'] for x in partial]}")
inbad = [x for x in usable if x["in_bad"] > 0]
print(f"    with SOZ-marked-bad conflicts  : {len(inbad)}  -> {[x['label'] for x in inbad]}")
print(f"\nevery distinct UNMATCHED token (count): {dict(all_unmatched)}")

eng_dist = Counter(str(x["engel"]) for x in usable)
print(f"\nEngel distribution (usable)        : {dict(eng_dist)}")
print("  (Engel I = seizure-free; needed to stratify 'is the clinician SOZ trustworthy')")

print("\n" + "=" * 70)
print("Audit done. Paste the full output back.")
print("=" * 70)