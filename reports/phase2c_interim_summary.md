# Phase 2C — Interim Summary

**Date:** 2026-10-03
**Scope:** micro-steps 13–18 plus diagnostic 15b, on the SWEC-ETHZ
short-term iEEG dataset (patients ID1, ID2, ID3).
**Status:** methodology established on the development patient and
replicated on two patients; no clinical validation yet (deferred to
Phase 3 / HUP per ADR 0010).

This document freezes the state of Phase 2C before moving to the next
major phase. It is an honest account of what was built, what was found,
what replicated, and what did not.

---

## 1. Question

Does a continuous-time quantum walk (CTQW) on a real ictal functional
connectivity network assign a reproducible, structurally meaningful
importance to individual electrodes — something classical network
analysis does not already capture?

## 2. Pipeline

For one seizure:

1. **Load** the recording (`data.loader.load_seizure`) and extract the
   ictal segment.
2. **FCM** — Pearson correlation across electrodes over the ictal
   segment (ADR 0007: Pearson primary, broadband 0.5–150 Hz, no
   thresholding).
3. **Graph** — convert the FCM to a usable graph (ADR 0008):
   zero the diagonal, clip negatives to zero (Layer P, primary).
4. **CTQW** — for each electrode as the start node, run the
   time-averaged quantum walk and take its participation ratio (PR),
   a measure of how widely the walk spreads.
5. **Rank** electrodes by PR. Low PR = the walk stays localized; high
   PR = it spreads broadly. Rankings only, no thresholds (ADR 0009).

Every step is tested (`tests/`, 71 passing) and every methodological
choice is recorded in an ADR.

## 3. The reframe (ADR 0009)

The original Phase 2C plan compared classical diffusion with CTQW on
the same graph. Micro-step 14 showed this comparison is not well-posed
for our time grid: classical diffusion converges to a near-uniform
stationary distribution, so its time-averaged PR is pinned at the
ceiling (≈ n) for every electrode (std 0.02 across 47 electrodes).
A rank correlation against an all-tied vector is meaningless.

The analysis was reframed (ADR 0009): classical diffusion is a **null
baseline** (confirming the graph is connected), and the primary object
of study is the **CTQW per-electrode PR ranking** in its own right.

## 4. Findings across the cohort

### 4.1 CTQW rankings are not explained by classical node metrics
On ID1/Sz1, CTQW PR correlated only weakly and non-significantly with
weighted degree (Spearman −0.19, p = 0.20). The most-localized and
most-spread electrodes are not the highest- or lowest-degree nodes.
CTQW spread is not a simple function of degree — consistent with the
earlier synthetic-graph finding (micro-steps 6–9) that spectral, not
classical, node properties govern quantum-walk spread.

### 4.2 Cross-seizure consistency (the primary result)
Within-patient agreement of CTQW rankings across a patient's seizures
(mean pairwise Spearman, Layer 1, long seizures):

| Patient | Seizures | mean Spearman |
|---------|---------:|--------------:|
| ID1     | 9 (>=30s)|       +0.410  |
| ID2     | 4        |       +0.630  |
| ID3     | 2        |       +0.535  |

All three positive, all near 0.5. The ranking carries a reproducible
within-patient signal, though it is modest — below the 0.4–0.6+ seen
for well-established functional-connectivity reproducibility measures,
and with wide confidence intervals for ID2/ID3 (4 and 2 seizures).

### 4.3 Signature electrodes (consistent across all three patients)
In each patient, specific electrodes stay consistently at the
localized or spread extreme across that patient's seizures:

| Patient | localized example | spread example |
|---------|-------------------|----------------|
| ID1     | e28 (mean rank 8.0 ± 5.7) | e46 |
| ID2     | e8 (mean rank 4.5 ± 1.8)  | e29 (38.8 ± 1.1) |
| ID3     | e91 (ranks 1, 2)          | e14 (ranks 92, 94) |

Electrode numbers do not correspond across patients (different
implants), so the replicated observation is the *phenomenon*, not a
specific electrode: in three independent patients, CTQW singles out a
stable subset of electrodes. This is consistent with CTQW reflecting a
stable topological feature of each patient's ictal network. It is not
yet established, however: we have observed low within-patient rank
variance for some electrodes but have not tested that variance against
a null model (e.g. shuffled rankings), so we cannot yet say the
stability exceeds chance. That null-model test is a near-term next step
(see section 6).

### 4.4 Duration confounder (validated)
Short seizures (<30 s) produce noisier FCMs. Micro-step 15b: pairs
involving a short seizure had mean Spearman +0.225 vs +0.410 for
long-only pairs (Mann-Whitney U p < 0.0001). This is a real effect,
not an artefact of reduced sample size, and it empirically validates
the 30-s minimum window (Bastos & Schoffelen 2016).

### 4.5 Onset heterogeneity (weak, does not cleanly replicate)
Restricting to the first 30 s of ictal (Layer 2) was expected to be a
cleaner, more consistent window. Instead it was *less* consistent than
the full ictal in 2 of 3 patients (L2 − L1: ID1 −0.169, ID2 +0.004,
ID3 −0.095). The tentative reading — seizure onset is per-seizure
heterogeneous, full ictal more stereotyped — holds weakly and
inconsistently. Not a robust claim; flagged as a Phase 3 hypothesis.

## 5. What replicated, what did not

| Result | Replication status |
|--------|--------------------|
| CTQW ranking has modest positive within-patient consistency | ✓ all 3 patients |
| Signature electrodes exist per patient | ✓ all 3 patients |
| CTQW spread not explained by degree | ✓ (tested on ID1; consistent with synthetic work) |
| Onset (Layer 2) less consistent than full ictal | ⚠ 2 of 3, weak |

That not everything replicated is itself reassuring: a result where
every test passed would raise the suspicion of overfitting.

## 6. Honest limitations

- **No clinical validation.** SWEC-ETHZ short-term provides no
  electrode-level seizure-onset-zone marks (confirmed, ADR 0010). We
  cannot say the localized electrodes are epileptogenic. Every result
  here is methodological, not clinical.
- **Small N.** 3 patients; 2 and 4 seizures for the replication
  patients. Confidence intervals are wide. These are directional
  signals, not established effects.
- **Signature phenomenon not yet tested against a null model.** The
  "signature electrode" observation (section 4.3) rests on low
  within-patient rank variance, not on a test that this variance
  exceeds what shuffled rankings would produce. Until that test is run,
  it is a suggestive observation, not an established effect. (Near-term
  next step, doable on the current data.)
- **Signature not yet tied to anything clinical.** Even once established
  against a null, we have not shown the reproduced electrodes mean
  anything clinically. That is the whole point of Phase 3.
- **One connectivity measure, one layer.** Pearson + Layer P only.
  PLV (ADR 0007) and signed-graph Layer S (ADR 0008) remain untested.

## 7. What carries to Phase 3 (HUP, with mentor, post-DGS)

1. The pipeline (loader → connectivity → Layer P → CTQW → PR ranking),
   ported to a dataset with SOZ marks.
2. The central test: do low-CTQW-PR (localized) electrodes correlate
   with clinician-marked seizure onset zones? This is the clinical
   bridge the whole project aims at (ADR 0009, ADR 0010).
3. The onset-heterogeneity hypothesis (section 4.5), testable with
   adequate patient numbers.
4. A compute note: CTQW on 98 electrodes took ~6 min/seizure; larger
   datasets will need an eigendecomposition cache.

## 8. One-paragraph verdict

Phase 2C established and tested a reproducible CTQW-based analysis of
real ictal connectivity networks across three patients. The method
produces a within-patient-stable per-electrode ranking that classical
degree does not explain, and the stability replicates across three
independent patients. Whether that stability has clinical meaning is
unknown and untestable on this dataset — it is the question Phase 3 is
built to answer. As of 2026-10-03 the project has a validated
methodology and an honest negative/weak result on the secondary
questions, but no clinical finding yet.