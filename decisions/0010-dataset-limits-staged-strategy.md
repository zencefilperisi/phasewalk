# ADR 0010 — Dataset limits for clinical validation and staged strategy

**Status:** accepted

**Date:** 2026-09-19

**Author:** project owner

---

## Context

Micro-step 14b produced the first CTQW-based per-electrode
localization pattern on a real ictal network (ID1/Sz1, Layer P).
The natural next question — "do the electrodes CTQW localizes
correspond to the clinician-identified seizure onset zone (SOZ)?" —
is what ADR 0009's Open Questions section already flagged as the
first genuine clinical bridge for the project.

To answer that question we need, for each patient, an electrode-level
label of which channels are marked as SOZ (and ideally the resection
zone and surgical outcome for validation). A literature search on
2026-09-19 confirmed:

1. **SWEC-ETHZ short-term dataset** provides `Sz*.mat` files with a
   single `EEG` variable and nothing else. No electrode-level SOZ
   marks, no channel labels, no electrode locations, no resection
   zone, no surgical outcome. This is documented at
   ieeg-swez.ethz.ch: "In each .mat file, the recording is stored
   in a variable called EEG as a TxM array". That is the complete
   per-seizure metadata.

2. **Burrello 2019 (TBME)** — the reference paper for this dataset —
   *does* report identification of "ictogenic brain regions", but
   explicitly "at two levels of granularity: cerebral hemispheres
   and lobes". This is coarser than the electrode-level labels we
   would need. A statement like "electrode 3 is inside the SOZ" is
   not derivable from what Burrello 2019 reports.

3. **HUP (Hospital of the University of Pennsylvania) iEEG dataset**
   does include electrode locations, clinician-identified SOZ marks,
   resection/ablation zones, and Engel outcome classifications. This
   is the kind of metadata a clinical-bridge analysis needs.

This means the plan implied by ADR 0009's Open Questions — CTQW-based
localization validated against SOZ marks on this dataset — is
**not executable on SWEC-ETHZ short-term**. The dataset choice
(fixed by ADR 0003) constrains what claims Phase 2C can make.

## Options considered

**A. Migrate to the HUP dataset now.**
Rebuild the loader, connectivity, and downstream pipeline on HUP.
Gains: real clinical validation becomes possible. Costs: 2–3 weeks
of re-work; several new steps (electrode-location handling, brain-
atlas alignment, montage schemes) require clinical judgment; the
solo working setup (no mentor until Feb–Mar 2027) makes those steps
risky to get right without oversight.

**B. Stay on SWEC-ETHZ and drop the clinical target.**
Keep the current dataset. Reframe Phase 2C as a purely methodological
result: "CTQW captures a per-electrode structure on real ictal
networks that classical graph metrics do not explain." No SOZ claim.
Gains: no re-work; the reframed micro-step 14b result already stands
as an honest observation. Costs: the paper narrative loses its
clinical hook, weakening publishability.

**C. Staged strategy — SWEC-ETHZ now, HUP after mentor.**
Keep SWEC-ETHZ as the methodology-development dataset. Build the
CTQW-on-ictal-network pipeline, robustness checks (Layer S per
ADR 0008), cross-seizure and cross-patient stability, all on
SWEC-ETHZ. Defer the SOZ-validation extension to HUP until a mentor
is available (post Feb–Mar 2027 DGS). Present current results as
methodological observation; present HUP validation as future work.

## Decision

Adopt **Option C** (staged strategy).

**Phase 2C on SWEC-ETHZ (now through DGS):**
- Central quantity: CTQW time-averaged participation ratio per
  electrode on the Layer P ictal FCM (per ADR 0009).
- Robustness checks: Layer S (ADR 0008), cross-seizure stability
  within a patient (ID1's 13 seizures), cross-patient replication
  (ID2 and ID3 per ADR 0006).
- Claim ceiling: purely methodological. Statements are of the form
  "CTQW localizes at electrodes X, Y, Z", never "electrodes X, Y, Z
  are epileptogenic." The clinical interpretation is explicitly out
  of scope for the SWEC-ETHZ portion of the work.

**Phase 3 on HUP (after DGS, with mentor):**
- Port the pipeline to HUP.
- Add SOZ-marked electrodes and resection zones as ground truth.
- Test whether CTQW localization patterns correlate with SOZ marks
  (rank-based analysis, per ADR 0009's no-thresholds rule).
- If the correlation holds, this is the clinical bridge and the
  scientific contribution the project was aiming at from the start.
- If it does not, the SWEC-ETHZ methodological result still stands
  and can be published on its own as a graph-theoretic observation
  about CTQW on real neural connectivity.

## Implications for other ADRs

- **ADR 0003** (dataset choice) stays as accepted. It correctly
  identified SWEC-ETHZ as a workable methodology-development dataset
  under our disk and structural constraints; it did not commit to
  SWEC-ETHZ being the *validation* dataset, and this ADR clarifies
  that split.
- **ADR 0006** (ID1 development, ID2/ID3 replication) stays in force
  for the SWEC-ETHZ phase.
- **ADR 0007** (FCM methodology) and **ADR 0008** (Layer P / Layer S)
  stay in force; they will need to be re-examined when the pipeline
  is ported to HUP, because HUP's sampling rate, filtering, and
  montage conventions differ.
- **ADR 0009** (Phase 2C analysis reframing) — its Open Question
  about SOZ correlation is now formally deferred to Phase 3 rather
  than being a live open question of Phase 2C.

## Consequences

**Easier:**
- No pipeline re-work now. Every piece of infrastructure already
  built (loader, connectivity, walks, tests) stays.
- The current results have a clear, defensible ceiling on what they
  claim, which removes the temptation to overreach in the write-up.
- Future mentor onboarding is cleaner: SWEC-ETHZ work is a self-
  contained methodological chapter with well-documented ADRs; HUP
  extension is a separate, additive step.

**Harder:**
- The scientific pay-off for the SWEC-ETHZ portion is limited by
  design. A reviewer will legitimately ask "so what?" of any purely
  methodological result. The write-up must address this by
  explicitly positioning the SWEC-ETHZ work as pipeline validation
  and by stating the HUP extension as the intended follow-up.
- The staged plan depends on the DGS outcome and mentor availability
  in Feb–Mar 2027. If either slips, Phase 3 slips with it.
- Some Phase 2C results (like "electrode 3 is the most localized in
  ID1/Sz1") will feel unfinished until the HUP validation confirms
  or refutes their clinical meaning. This is a real cost.

**Commitments:**
- No SOZ-related claims are made using SWEC-ETHZ. Any language in
  micro-step scripts, notebooks, or preprints that could be read as
  a clinical claim on SWEC-ETHZ data is either rewritten or removed.
- The HUP extension is a distinct future phase (Phase 3), scheduled
  post-mentor. It is not attempted solo.
- If circumstances change (e.g. a mentor becomes available earlier,
  or the SWEC-ETHZ metadata expands), this ADR is revisited.

## Open questions to revisit

- What exactly does the HUP dataset look like in terms of preprocessing,
  sampling rate, montage convention, and SOZ-labeling format? This
  determines how much of the current loader / connectivity code
  survives the port. Deferred until mentor conversations begin.
- Is there a smaller, faster clinical-validation option than a full
  HUP port — e.g. a SWEC-ETHZ patient whose SOZ was published in a
  separate clinical paper we could match by hand? Low-probability
  lead; noted for future search.
- Should the "methodological observation" from SWEC-ETHZ be prepared
  as a standalone preprint independent of the HUP work? To be
  decided once cross-seizure and cross-patient replication (ID2,
  ID3) results are in.