# ADR 0001 — Dataset choice for Phase 2

**Status:** accepted

**Date:** 2026-05-16

**Author:** project owner (with literature search support)

---

## Context

Phase 2 of the project requires real iEEG recordings from epilepsy
patients. The aims of Phase 2 are:

- to build functional connectivity matrices from the iEEG channels of
  each patient,
- to run classical random walks and continuous-time quantum walks
  (CTQW) on those matrices,
- to compare how the two propagation models behave across different
  brain states (interictal, pre-ictal, ictal).

This requires a dataset with:

1. **Long enough recordings** to extract stable functional connectivity
   estimates across multiple time windows.
2. **Multiple patients** with verified seizure events, to allow
   cross-patient statistical analysis.
3. **Reasonably high sampling rate** (>= 256 Hz) so that we can study
   frequencies up to the gamma band.
4. **Clinically annotated seizure onset and offset times**, so we can
   slice the recordings into interictal / pre-ictal / ictal segments.
5. **Open licence**, no IRB requirement at the user's institution,
   compatible with the user not yet having a formal academic affiliation.

## Decision

We use the **SWEC-ETHZ long-term iEEG database**, hosted by the
Sleep-Wake-Epilepsy-Center of Inselspital Bern and the Integrated
Systems Laboratory of ETH Zurich.

Source: <http://ieeg-swez.ethz.ch/>

Key properties:
- 18 patients with pharmacoresistant epilepsy
- 116 seizures across 2656 hours of continuous recording
- 24-128 channels per patient (strip, grid, and depth electrodes)
- Sampling at 512 or 1024 Hz
- Anonymised, freely available for research and education
- Band-pass filtered 0.5-150 Hz; common-median referenced
- Format: MATLAB `.mat` files, one per hour of recording

## Alternatives considered

| Candidate | Why rejected |
|-----------|--------------|
| **SWEC-ETHZ short-term subset** | Each recording is only ~10 minutes (pre-ictal + ictal + post-ictal). Too short to estimate stable functional connectivity over multiple non-overlapping windows. Good for benchmarking but not for the propagation analysis we need. |
| **OpenNeuro ds003876 (RNS Epilepsy iBIDS)** | Stimulation-paired iEEG, only 8 patients, single-day recordings. The stimulation design is interesting clinically but introduces a confound we cannot easily handle. |
| **OpenNeuro ds003844 (Clinical Epilepsy iEEG to BIDS)** | 12 patients, BIDS-compatible, but seizure counts and recording durations are not consistently large. BIDS is a plus but not enough to compensate for the smaller cohort. |
| **OpenNeuro ds005602 (IDEAS)** | Imaging only (T1w, FLAIR). No iEEG. Outside scope. |
| **OpenNeuro ds005873 (SeizeIT2)** | Wearable EEG (behind-the-ear), not iEEG. Spatial resolution and channel count insufficient for the network analysis we plan. |
| **Human Connectome Project** | Healthy population, no seizures. We initially considered using it for "healthy baseline" but Phase 2 needs within-patient ictal/interictal contrast, which HCP cannot provide. |

## Consequences

**Easier:**
- Standard dataset, used in ~30% of recent epilepsy iEEG papers — a
  reviewer recognises it instantly.
- Long enough for stable connectivity estimates in many time windows.
- Variety of patients enables a basic statistical analysis (within-
  patient ictal vs interictal comparisons; possibly across-patient
  pattern aggregation).

**Harder:**
- Format is `.mat` rather than BIDS — we need to write our own loader.
  Not a serious obstacle but adds setup time.
- Electrode counts and locations vary widely across patients. We
  cannot directly average patients in electrode space; analyses must
  either be within-patient or in a derived space (e.g. graph metrics
  rather than per-channel quantities).
- Total dataset size is large (~500-800 GB). We will start with 2-3
  patients (~30-50 GB) and expand only if needed.
- The dataset has no formal "control" group of healthy controls. All
  contrasts are within-patient (state contrasts: ictal vs interictal).

**Commitments resulting from this decision:**
- All analyses in Phase 2 will be designed at the within-patient
  level first, with cross-patient aggregation as a secondary step.
- The data loader becomes the first piece of Phase 2 code; until it
  is reliable and tested, no analysis code is written.
- We must contact the SWEC-ETHZ group (`ieeg@iis.ee.ethz.ch`) early
  to acknowledge usage and ask whether they have any guidance for
  our intended analysis. This is courtesy, not legal requirement,
  but it lays the ground for collaboration.

## Open questions to revisit later

- Should we also include a structural connectome (DTI/dMRI) for any
  of these patients? SWEC-ETHZ does not include structural imaging.
  If our analysis later requires it, we will need a separate dataset
  and a way to bridge the two — a non-trivial decision deferred to
  a future ADR.
- How many patients do we actually need for adequate statistical
  power? This depends on effect sizes we have not yet measured.
  Revisit once first results are available.
