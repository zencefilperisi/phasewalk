# ADR 0003 — Switch from the long-term to the short-term SWEC-ETHZ dataset

**Status:** accepted (supersedes the dataset choice in ADR 0001)

**Date:** 2026-07-23

**Author:** Hatice Nalçacı

---

## Context

ADR 0001 chose the **long-term** SWEC-ETHZ iEEG dataset. That decision
was made before we had examined the actual file structure of the
database or accounted for the project owner's available disk space.
Two facts have since come to light:

1. **Disk constraint.** The project owner currently has ~44 GB of free
   disk space. A single long-term patient spans up to 294 hourly
   `.mat` files (one hour each, each up to ~1 GB), so a single
   long-term patient can approach or exceed the entire available disk.
   Downloading 2-3 long-term patients is not feasible right now.

2. **Better structural fit of the short-term set.** Inspecting the
   database page revealed that the **short-term** dataset is organized
   exactly around the object of our study:
   - Each patient is a single `.zip` (< 1 GB), containing one `.mat`
     file per recorded seizure (`Sz1.mat`, `Sz2.mat`, ...).
   - Each seizure file contains **3 minutes of pre-ictal**, the
     **ictal segment** (10 s to ~1000 s), and **3 minutes of
     post-ictal** signal.
   - The seizure onset is at a known, fixed sample (512 x 3 x 60 =
     92160), so pre-ictal / ictal / post-ictal segments are trivially
     separable.

Our analysis question is about seizure **propagation** — how activity
spreads across the iEEG network before, during, and after a seizure.
The short-term set hands us exactly those epochs, pre-segmented, with
the onset already marked by an expert epileptologist. The long-term
set would require first *finding* the seizures within hundreds of
hours of recording.

## Decision

Switch to the **short-term SWEC-ETHZ dataset** for Phase 2.

- Source: <http://ieeg-swez.ethz.ch/> (Short-term Dataset section)
- 16 patients, 100 seizures total
- 36-100 electrodes per patient, 512 Hz sampling
- Each seizure: 3 min pre-ictal + ictal + 3 min post-ictal
- Per-`.mat` `EEG` variable is a **T x M** array where, for the
  short-term set, **T = number of sampling points (time)** and
  **M = number of electrodes**. (Note: this is transposed relative to
  the long-term set, where T = electrodes and M = time. We must handle
  the orientation explicitly in the loader.)

Start with patients **ID1, ID2, ID3** (a few GB total), expand later
if needed.

## Alternatives considered

| Option | Why not (now) |
|--------|---------------|
| **Stay with long-term** (ADR 0001) | Disk cannot hold even one patient comfortably; seizures must be located within hundreds of hours; overkill for a propagation analysis that only needs peri-ictal epochs. May revisit later for resting-state/inter-ictal connectivity if disk allows. |
| **Zenodo curated subset** (3 short-term patients + 5 h of long-term ID18) | Reasonable, but the canonical short-term set is the standard reference and gives us the full 16-patient cohort to expand into. |
| **A different database entirely** | No reason to abandon SWEC-ETHZ; it remains a recognized standard. Only the long-vs-short choice changes. |

## Consequences

**Easier:**
- Fits the disk budget (each patient < 1 GB).
- Pre-segmented peri-ictal epochs — no seizure-hunting in long recordings.
- Expert-marked onset times come with the data.
- 16-patient cohort still allows cross-patient analysis.

**Harder / watch out for:**
- **Array orientation is transposed** relative to the long-term set.
  The loader must not assume a fixed orientation; it must check which
  axis is electrodes (M, expected 36-100) and which is time
  (T, expected ~512 x duration).
- Short recordings (a few minutes) mean fewer independent time windows
  for connectivity estimation than long-term would give. We will need
  to estimate functional connectivity on short epochs, which is
  statistically noisier; this must be acknowledged as a limitation.
- No inter-ictal "baseline far from any seizure" — only 3 min pre-ictal.
  If we later need true interictal baseline, we may have to add
  long-term data (revisit then).

**Commitments resulting from this decision:**
- The first piece of Phase 2 code is a robust `.mat` loader that
  detects array orientation and exposes (electrodes x time), the
  sampling frequency, and the known onset/offset sample indices.
- Until the loader reliably reads ID1-ID3 and reports sane shapes,
  no connectivity or walk code is written on real data.
- ADR 0001's "open data, no IRB" reasoning still holds and is unchanged.

## Open questions to revisit later

- Are 3-minute pre-ictal epochs long enough for stable functional
  connectivity estimates? If not, we either shorten the connectivity
  window (noisier) or bring in long-term inter-ictal data (more disk).
  Defer until we see real connectivity matrices.
- How to define functional connectivity (correlation, coherence,
  phase-locking value, ...)? This is the next methodological decision
  and will get its own ADR once we have data in hand.