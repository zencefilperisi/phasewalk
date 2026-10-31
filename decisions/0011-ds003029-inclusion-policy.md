# ADR 0011 — ds003029 channel and SOZ inclusion policy

**Date:** 2026-10-30
**Status:** Accepted
**Phase:** Phase 3 (clinical validation) — move from SWEC-ETHZ to the
Epilepsy-iEEG-Multicenter dataset (OpenNeuro ds003029).

## Context

SWEC-ETHZ short-term has no electrode-level seizure-onset-zone (SOZ) marks, so
clinical validation was impossible there (ADR 0010). ds003029 (Li et al.,
neural-fragility study; 37 published subjects across JHH / NIH / UMF / UMMC;
BrainVision iEEG; clinician SOZ contacts and surgical outcome in
`sourcedata/clinical_data_summary.xlsx`) provides both a larger cohort and the
SOZ labels needed for the central clinical test:

> do low-CTQW-PR (localized) electrodes correspond to the clinician SOZ?

Three read-only data-integrity audits (`ds003029_step1/2/3`) surfaced the
parsing and inclusion issues this ADR settles, so that no clinical result rests
on a silent name-matching or double-counting error. The methodology itself
(FCM → Layer P → CTQW → PR ranking) is frozen from Phase 2C and unchanged; this
ADR governs only how ds003029 data enters it.

## Decisions

**1. SOZ contact parsing.**
`soz_contacts` strings are split on `;`, `,` and `:`; tokens of the form
`<prefix><a>-<b>` are expanded to individual contacts (`<prefix>a` … `<prefix>b`);
matching to channel names is case-insensitive. Tokens that cannot be parsed, or
that correspond to no channel in the subject's montage, are DROPPED and logged
— never guessed. Example: ummc002 `att1-3` has no `ATT` channel in its montage
(ANT/MEST/SUBT only) — a clinical-table artifact — so those 3 contacts are
dropped; the subject keeps its 7 valid SOZ contacts.

**2. Subject-id mapping (folder ↔ clinical xlsx).**
The xlsx `dataset_id` column is inconsistently formatted (e.g. `ummc_001` with
an underscore vs folder `ummc001`). Matching normalizes both sides by
lower-casing and removing non-alphanumeric characters. This recovers subjects
(ummc001, ummc007) that an exact match missed.

**3. Duplicate subjects.**
`sub-pt01` and `sub-pt1` have identical montages (98 channels, 84 good, same
names, same 4 ictal runs) and pt01 has no distinct clinical row; pt01 is treated
as a duplicate of pt1 and EXCLUDED, to avoid double-counting one patient.
Inferred from montage identity, not byte-identity (a checksum/event-time check
can confirm later); excluding a possible duplicate is the conservative choice.

**4. Channel inclusion for analysis.**
Only channels with `type ∈ {ECOG, SEEG}` AND `status == good` (per the BIDS
`channels.tsv`) enter the FCM/CTQW pipeline. `status == bad` channels — whether
non-brain (WM/ventricle/CSF/out) or quality-flagged (`excel-file-labeled`) — are
excluded entirely, because a bad-quality signal cannot yield a trustworthy FCM
or CTQW value.

**5. SOZ reference over good channels only.**
A subject's SOZ reference set is (parsed SOZ contacts) ∩ (good channels). SOZ
contacts marked `status=bad` are dropped from the reference and logged per
subject (observed: jh103:1, pt15:2, pt8:3, umf004:6, ummc003:1, ummc008:1; all
were ECOG/SEEG flagged `excel-file-labeled`, i.e. quality-flagged brain
contacts, not non-brain). Rationale: CTQW can only be scored on channels we
actually analyze; a SOZ contact we cannot analyze cannot be scored. The
alternative — overriding `status` to force bad SOZ contacts into the analysis —
would feed corrupted signal into the FCM and is rejected.

**6. Usable cohort.**
A subject is usable for the clinical test iff it has ≥1 ictal run, a matched
clinical row, and ≥1 SOZ contact matched to a good channel. The audit gives ~32
usable; id-normalization (decision 2) is expected to recover ummc001/ummc007
(→ ~34), and pt01 is excluded (decision 3). Of the usable subjects, 20 are
Engel I (seizure-free).

## Consequences / honest caveats

- Dropped SOZ contacts (unparseable, non-existent, or bad-marked) slightly
  shrink each subject's SOZ reference. Per-subject counts are logged so the loss
  is visible, not silent.
- The clinician SOZ is a **hypothesis**. In Engel > I (non-seizure-free)
  patients that hypothesis was, in hindsight, at least partly wrong — the surgery
  did not stop the seizures. The clinical analysis will therefore be stratified
  by outcome: the Engel-I subgroup (n ≈ 20) is where the SOZ is the most
  trustworthy reference. Reporting correspondence pooled across all outcomes
  without stratifying would mis-state the result.
- pt01's exclusion rests on montage identity, not proven byte-identity.

## Deferred

- The **ictal window definition** (which samples constitute the ictal segment,
  derived from the inconsistent electrographic onset/offset markers in
  `events.tsv`) is a separate decision, pending the events audit (ds003029
  step 4), to be recorded in its own ADR.

# Note added 2026-10-30 to ADR 0011 (umf002–005 exclusion)

Four UMF subjects (umf002, umf003, umf004, umf005) are excluded from the
ds003029 cohort because **the published BrainVision signal is shorter than the
annotated seizure**, so the ADR-0012 ictal window [onset, onset+30 s] cannot be
formed.

## Verification: dataset property, not an incomplete download

This was checked against the OpenNeuro S3 source directly:

    s3://openneuro.org/.../sub-umf004_..._ictal_acq-seeg_run-01_ieeg.eeg
      = 355,472 bytes on S3  ==  355,472 bytes locally

i.e. the local file is a complete copy of what the server publishes — not a
truncated transfer. The `events.tsv` for these runs marks the seizure onset far
beyond the end of the published signal (e.g. umf004: onset at 135 s, while the
~355 KB SEEG file holds only a few seconds). The UMF SEEG files carry a later
processing date (Apr 2021) than the rest of the dataset (Dec 2020); the most
likely cause is that the signal was trimmed during the EDF→BrainVision
conversion while `events.tsv` kept the original EDF timeline.

umf001 (99 MB signal, onset 60 s, offset 91 s — intact) is **retained**.

## Effect on the cohort

- usable cohort: 33 → **29**
- Engel I: 19 → **17** (umf004, umf005 removed)

umf004 and umf005 have SOZ fractions of 0.66 and 0.50, where per-subject AUC is
near-undefined and contributes essentially noise, so their removal does not
change the primary result (Engel-I median AUC 0.39, one-sided Wilcoxon
p = 0.96). The cohort-build step now validates, per run, that the ictal window
falls within the published signal, so the frozen cohort equals the analysed
cohort (no silent NaNs downstream).