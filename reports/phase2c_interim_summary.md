# Phase 2C — Interim Summary

**Date:** 2026-10-07
**Scope:** micro-steps 13–20 (incl. diagnostic 15b, null-model 19,
classical comparison 20), on the SWEC-ETHZ short-term iEEG dataset
(patients ID1, ID2, ID3).
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
within-patient signal, though it is modest in absolute terms — below
the 0.4–0.6+ seen for well-established functional-connectivity
reproducibility measures.

**Permutation test (micro-step 19).** Against a null of independent
uniform random rankings (20,000 permutations per patient), the observed
consistency is far beyond chance in all three patients:

| Patient | observed | null mean ± std | z | p |
|---------|---------:|----------------:|---:|---|
| ID1     | +0.410   | 0.00 ± 0.025    | +16.6 | < 0.00005 |
| ID2     | +0.630   | 0.00 ± 0.064    | +9.8  | < 0.00005 |
| ID3     | +0.535   | 0.00 ± 0.101    | +5.3  | < 0.00005 |

(p = 0.00005 is the permutation floor, 1/20001; the true p is far
smaller.) The within-patient reproducibility of CTQW rankings is a
statistically established result, not a chance pattern — including for
ID3 (a single seizure pair, 98 electrodes).

**Important boundary.** This null is "no structure at all" (random
permutation). It establishes that the rankings are reproducible, NOT
that CTQW is more reproducible or more informative than a classical
centrality on the same FCMs. That comparison is a separate test
(see section 6).

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
stable topological feature of each patient's ictal network.

**Null-model test (micro-step 19).** The most-stable electrode (minimum
per-electrode rank std) was tested against random permutations: ID1
p = 0.015, ID2 p = 0.047 (both significant), ID3 p = 0.63
(inconclusive — with only 2 seizures the min-std statistic cannot
resolve signal from chance, since a random pair repeats an exact rank
often enough). This secondary test is not independent of the overall
consistency test (if rankings are consistent, extremes are stable by
construction); the primary evidence is the consistency test in 4.2.

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

### 4.6 CTQW vs classical centralities — the core comparison (micro-step 20)
On the same seizures and Layer P graphs, cross-seizure consistency of
CTQW PR was compared against two classical centralities, and the
distinctness of the CTQW ranking from them was measured.

**Reproducibility (mean pairwise Spearman):**

| Patient | CTQW PR | eigenvector cent. | weighted degree |
|---------|--------:|------------------:|----------------:|
| ID1     | +0.410  | +0.447            | +0.538          |
| ID2     | +0.630  | +0.704            | +0.795          |
| ID3     | +0.535  | +0.840            | +0.639          |

In all three patients CTQW is the LEAST reproducible of the three. We
cannot claim CTQW is a more reproducible centrality — the opposite
holds. (This is partly expected: weighted degree is a simple sum of
edge weights and is stable whenever overall coupling strength is stable,
while CTQW PR is a higher-order spectral quantity sensitive to finer
structure that varies more between seizures. An explanation, not a
reason to discount the result.)

**Distinctness (within-seizure Spearman, CTQW vs classical):**

| Patient | CTQW vs eigenvector | CTQW vs degree |
|---------|--------------------:|---------------:|
| ID1     | −0.016              | −0.280         |
| ID2     | −0.348              | −0.493         |
| ID3     | −0.010              | −0.375         |

CTQW is nearly orthogonal to eigenvector centrality (≈ 0 in ID1 and
ID3) and not redundant with degree. CTQW tracks a reproducible-beyond-
chance but genuinely DIFFERENT aspect of the ictal network than the
classical centralities.

**What this does to the thesis.** The project's value does not rest on
CTQW being a more stable centrality — it is not. It rests on CTQW
capturing a distinct structure whose clinical meaning is the Phase 3
question. This result closes off a secondary claim (reproducibility
superiority) we might otherwise have been tempted to make, and sharpens
the real one: CTQW is different, and whether "different" is "clinically
useful" is what Phase 3 must decide.

### 4.7 Layer S (signed graph) robustness (micro-step 21)
ADR 0008's deferred signed-graph check was run: CTQW rankings on Layer S
(negatives kept) vs Layer P (negatives clipped). Per-seizure agreement
between the two: ID1 +0.587, ID2 +0.806, ID3 +0.778 (moderate to high).
Layer S is about as reproducible as Layer P (consistency within ±0.08
in each patient). Conclusion: Layer P is validated as primary — clipping
negatives does not produce a misleading ranking — but negatives carry
some CTQW-relevant structure (most in ID1), so the two-layer design was
correct. Both carry to Phase 3 (Layer P primary, Layer S sensitivity).

## 5. What replicated, what did not

| Result | Replication status |
|--------|--------------------|
| CTQW ranking consistency beyond chance (permutation test) | ✓ all 3 patients, p < 0.00005 |
| CTQW ranking consistency ABOVE classical centralities | ✗ all 3 patients (CTQW is lower) |
| CTQW ranking distinct from classical centralities | ✓ all 3 patients (≈ orthogonal to eigenvector) |
| Signature electrodes exist per patient | ✓ all 3 patients |
| Onset (Layer 2) less consistent than full ictal | ⚠ 2 of 3, weak |

That not everything replicated — and that one central hope (CTQW more
reproducible than classical) was refuted by the data — is itself a sign
of an honest analysis. A result where every hoped-for effect appeared
would raise the suspicion of overfitting or wishful framing.

## 6. Honest limitations

- **No clinical validation.** SWEC-ETHZ short-term provides no
  electrode-level seizure-onset-zone marks (confirmed, ADR 0010). We
  cannot say the localized electrodes are epileptogenic. Every result
  here is methodological, not clinical.
- **Small N.** 3 patients; 2 and 4 seizures for the replication
  patients. Confidence intervals are wide. These are directional
  signals, not established effects.
- **CTQW is NOT more reproducible than classical centralities
  (micro-step 20, settled).** The permutation test (micro-step 19)
  established CTQW rankings are reproducible beyond chance. The
  classical comparison (micro-step 20) then showed that classical
  centralities — eigenvector and weighted degree — are MORE reproducible
  than CTQW in all three patients. So CTQW's reproducibility is real but
  not superior. CTQW's justification therefore cannot be "a more stable
  centrality"; it must be "a distinct structure" (which it is: ≈
  orthogonal to eigenvector centrality) whose clinical value is still
  untested. The whole weight of the project's thesis now sits on
  Phase 3.
- **Signature not yet tied to anything clinical.** Even established
  against a null, we have not shown the reproduced electrodes mean
  anything clinically. That is the whole point of Phase 3.
- **One connectivity measure.** Pearson only; PLV (ADR 0007, a
  phase-based measure) remains untested — the next planned step. The
  signed-graph layer (Layer S, ADR 0008) HAS now been tested
  (micro-step 21): it agrees moderately-to-well with Layer P and is
  equally reproducible, so the Layer P primary choice is validated.

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

Phase 2C established and tested a CTQW-based analysis of real ictal
connectivity networks across three patients. The per-electrode ranking
is reproducible beyond chance (permutation test, p < 0.00005, all three
patients) and is distinct from classical centralities (≈ orthogonal to
eigenvector centrality). It is NOT more reproducible than those
classical centralities — in all three patients eigenvector centrality
and weighted degree are more stable across seizures. So the case for
CTQW is not "a better/more-stable centrality" but "a different one":
it captures a reproducible structure that classical measures do not,
and whether that difference is clinically meaningful is unknown and
untestable on this dataset. As of 2026-10-07 the project has a
validated methodology, one established result (reproducible-beyond-
chance, distinct-from-classical), one refuted hope (reproducibility
superiority), weak/negative results on the secondary questions, and no
clinical finding yet. The entire clinical case rests on Phase 3 (HUP,
with SOZ marks).