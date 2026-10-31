# ADR 0013 — ds003029 clinical-validation analysis plan (pre-registration)

**Date:** 2026-10-30
**Status:** Accepted (pre-registered BEFORE any cohort-level result was computed)
**Phase:** Phase 3 (ds003029 clinical validation)

## Why this is written now

This ADR fixes the entire analysis — metric, population, hypothesis, statistics,
and decision criteria — **before** any cohort-level number is looked at. The
only CTQW-vs-SOZ number produced so far is a single-subject integration check
(jh105, mean-PR AUC ≈ 0.51), used to verify the pipeline runs, not as evidence.

Pre-registering addresses the specific risk raised in Phase 2C: that we might,
consciously or not, steer the analysis toward the result we hope for. By
committing to the plan here and reporting whatever comes out, a null result
stays as informative and publishable as a positive one. We will NOT change the
primary metric or primary population after seeing results; post-hoc choices, if
any, will be labelled exploratory.

## Cohort (frozen before computing AUCs)

Per ADR 0011/0012, a subject is **usable** iff it has: a matched clinical row, at
least one COMPLETE ictal run (header + .eeg + channels + events), and ≥1 SOZ
contact matched to a good channel. Excluded: pt01 (duplicate of pt1, ADR 0011);
jh107 (no clinical row); pt1 (signal metadata-only on the server — if its signal
remains unavailable it is dropped, not reconstructed from pt01). The usable list
and each subject's SOZ set are enumerated and frozen by a cohort-build step
before any AUC is computed. Expected ≈ 32 usable, ≈ 20 Engel I (final N recorded
with the results once the download is complete).

## Per-subject quantity

For each complete ictal run of a subject:
1. ictal window [onset, onset+30 s] (ADR 0012),
2. Pearson FCM (ADR 0007), Layer P graph (ADR 0008),
3. CTQW time-averaged participation ratio (PR) per electrode (H = A; t-grid
   linspace(0.1, 100, 300); the batched eigendecomposition routine, equivalent
   to Phase 2C).

Per electrode, PR is averaged across the subject's runs (common good channels).
The **per-subject statistic** is

    AUC_i = AUC of score(-PR) predicting SOZ membership over subject i's good channels,

i.e. the probability that a random SOZ electrode is more localized (lower PR)
than a random non-SOZ electrode. AUC is used because it is invariant to SOZ
prevalence, which varies widely across subjects (jh105 has 31 % SOZ). AUC > 0.5
means SOZ electrodes tend to be more localized under the CTQW.

## Primary hypothesis and test

- **Population:** Engel-I subjects only. Rationale (ADR 0011): the clinician SOZ
  is a hypothesis that was, in hindsight, borne out only where surgery made the
  patient seizure-free; in Engel > I patients the SOZ is not a trustworthy
  localization target.
- **H1:** the per-subject AUCs of the Engel-I group are greater than 0.5.
- **Test:** one-sample, one-sided Wilcoxon signed-rank on {AUC_i − 0.5},
  α = 0.05. Effect reported as median AUC (IQR) and mean AUC.
- This is the single confirmatory test. Everything below is secondary and will
  not be treated as confirmatory even if significant.

## Secondary analyses (exploratory, labelled as such)

1. **All-outcomes pooled:** same test on every usable subject.
2. **Engel > I subgroup:** interpretive contrast — a weaker/absent effect here
   would be consistent with the SOZ-trustworthiness argument (not a formal
   control).
3. **Classical-centrality baseline:** recompute AUC_i using eigenvector
   centrality, weighted degree and betweenness in place of CTQW PR, on the same
   graphs. Paired Wilcoxon (CTQW vs each classical) over Engel-I subjects.
   Question: does CTQW localize the SOZ *better than* a classical centrality?
   **Pre-specified score directions** (fixed here to remove a degree-of-freedom):
   CTQW uses low PR (score = −PR, the localization hypothesis); the three
   classical centralities use high centrality (score = +value, the "epileptic
   hub" hypothesis). AUC is read in these fixed directions for every subject;
   we do not pick, per method, whichever direction scores higher.
   **Pre-committed expectation:** based on the Phase 2C result (CTQW was distinct
   from but not superior to classical centralities), we do NOT expect CTQW to
   beat classical, and will report the comparison either way.
4. **Label-permutation null:** per subject, shuffle SOZ labels (preserving SOZ
   count) many times, recompute AUC, to confirm the observed AUC distribution is
   or is not beyond chance (the clinical analogue of micro_step 19).
5. **Per-electrode pooled view:** AUC pooling all electrodes with subject as a
   grouping factor — exploratory only.

## Sensitivity analyses (pre-specified)

- **Short seizures:** repeat the primary test excluding runs with < 30 s ictal
  (ADR 0012 duration flag), since short windows give noisier FCMs
  (micro_step_15b).
- **SOZ fraction (reported, never filtered):** AUC is prevalence-robust by
  construction, so NO subject is excluded on SOZ size. Excluding large-SOZ
  subjects would drop the hardest cases (where AUC is naturally near 0.5, as
  jh105 showed) and could inflate the result — precisely the researcher
  degree-of-freedom bias this pre-registration guards against. Instead, AUC is
  reported against SOZ fraction descriptively so any dependence is visible. The
  only subjects that cannot yield a defined AUC are those with no non-SOZ good
  channels (SOZ = 100 %); such a subject, if any exists, is noted and omitted of
  necessity, not by choice.
- **Layer S (ADR 0008) and PLV (ADR 0007):** rerun the primary test under the
  signed layer and the phase-based measure as robustness, if and only if the
  primary 30 s / Layer P / Pearson result warrants it (these multiply compute).

## Decision criteria (pre-committed)

- **Support for the clinical hypothesis** = Engel-I median AUC meaningfully
  above 0.5 AND Wilcoxon p < 0.05 AND the effect survives the short-seizure
  sensitivity.
- **No support** = the Engel-I AUC distribution is ≈ 0.5 or the test is
  non-significant. This is reported, plainly, as a negative result — a
  legitimate and publishable outcome, not a failure to be worked around.
- We will not switch the primary metric (AUC), primary population (Engel I), or
  primary test (one-sided Wilcoxon) after seeing results.

## What is explicitly NOT claimed

- No causal or clinical-utility claim. A positive result is a statistical
  association, in this cohort, between CTQW-PR localization and the clinician
  SOZ — not evidence that CTQW should guide surgery.
- No claim that CTQW is clinically deployable, nor (unless the baseline analysis
  shows it) that CTQW outperforms classical centralities.

## Recording

Results will be written to a Phase 3 summary citing this ADR as the
pre-registration. The cohort list and SOZ sets are frozen by the cohort-build
step before AUC computation; that frozen list is saved alongside the results.