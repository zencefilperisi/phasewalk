# ADR 0005 — Patient-level analysis, not group-level

**Status:** accepted

**Date:** 2026-07-27

**Author:** Hatice Nalçacı

---

## Context

The three patients we begin Phase 2 with (ID1, ID2, ID3) differ
substantially:

| patient | electrodes | seizures | epilepsy profile |
|---------|-----------:|---------:|------------------|
| ID1     | 47         | 13       | many seizures, highly variable duration |
| ID2     | 42         | 4        | few but long seizures                   |
| ID3     | 98         | 2        | very high electrode density, few seizures |

These differences are not noise — they reflect real clinical
heterogeneity in the underlying epilepsy and in the pre-surgical
implantation plan tailored to each patient.

A conventional "group analysis" (compute a metric per patient, then
average or compare across patients) is tempting but does not work here:

- **Matrix-size confound.** Functional-connectivity matrices and their
  spectral properties depend on their size. A 42x42 matrix and a
  98x98 matrix produce different eigenvalue distributions purely
  because of dimension. Comparing raw numbers across such matrices
  conflates a size effect with a biological effect.
- **Statistical power.** Three patients cannot support group-level
  statistics in any meaningful sense.
- **Anatomical incomparability.** Electrode locations are
  patient-specific; there is no shared "channel 12" across patients
  the way there might be a shared "region V1" in a group fMRI study.

## Decision

Phase 2 analyses are conducted at the **within-patient level**. Each
patient is analysed independently; results are compared across patients
**qualitatively**, not statistically.

Concretely:

1. All within-patient statistics (e.g., "in ID1, seizures Sz1 through
   Sz13 show a consistent propagation pattern with correlation X")
   are permitted and are the primary output.
2. Across-patient statements are limited to qualitative comparisons
   ("the same pattern is also observed in ID2 and ID3") or
   presence/absence claims ("the effect that appears in ID1 does not
   appear in ID3").
3. No number computed by averaging across patients will be reported
   as a headline finding. We may compute such numbers privately for
   sanity checks but will not present them as evidence.
4. Any dimensionless, size-invariant metric (e.g., normalized
   participation ratio expressed as a fraction of the maximum
   possible) may be *compared numerically* across patients; anything
   size-dependent (raw eigenvalue spectra, raw connectivity strengths)
   may not.

## Alternatives considered

| Option | Why rejected |
|--------|--------------|
| Group-level averaging across patients | Combines incompatible matrix sizes; three patients is too few for meaningful group statistics; anatomical channels are not comparable. |
| Reduce all patients to a common electrode subset | Would throw away most of ID3's data (from 98 down to ~42 channels) with no principled basis for choosing which electrodes to keep. |
| Register all patients to a common brain atlas via electrode coordinates | Correct in principle but requires anatomical MRI data we do not have and is outside the project's scope. Deferred. |

## Consequences

**Easier:**
- Every patient stands on its own; a null result in ID3 does not
  invalidate a positive result in ID1.
- No pretence of statistical rigor at group level.
- Analyses can be tuned to each patient's actual data without
  worrying about cross-patient consistency of parameters.

**Harder:**
- Publication story becomes "we show this pattern in ID1, replicated
  in ID2, and observed differently in ID3" — this is a weaker claim
  than "we show this pattern across a cohort".
- Reviewers may ask why we did not expand to more patients. The
  answer must be prepared: disk constraint, methodological
  development stage, cross-patient normalization is a separate
  research problem.
- All metrics must be checked for size-dependence before being
  compared across patients.

**Commitments:**
- Every analysis notebook produces per-patient output; cross-patient
  summaries appear only in a final synthesis document.
- The paper's Methods section will state explicitly that this is a
  within-patient study and justify the choice.
- If Phase 2 expands beyond three patients later, this ADR may be
  revised, but the within-patient primacy is preserved.

## Open questions to revisit

- If a Phase 2 result turns out to depend on matrix size in a way
  we can quantify and correct for, cross-patient comparison could
  become meaningful. Reconsider then.
- If we later obtain a larger cohort (e.g., all 16 short-term
  patients), some group-level statistics may become possible.
  Revisit then.