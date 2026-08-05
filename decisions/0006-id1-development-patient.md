# ADR 0006 — ID1 as the development patient, with parameter independence

**Status:** accepted

**Date:** 2026-07-30

**Author:** Hatice Nalçacı

---

## Context

We have three patients (ID1: 13 seizures, ID2: 4, ID3: 2). The
methodology of Phase 2 — how to build functional connectivity matrices,
which spectral / walk metrics to compute, which time windows to use —
is not yet fixed and will be developed iteratively as we look at real
data.

Two roles a patient can play in this process:

- **Development role.** The data on which the analyst inspects results
  intermediate outputs, tries variants, decides parameter values.
- **Test / replication role.** The data on which a *fixed* analysis
  pipeline is applied to see whether the findings from the development
  role hold up.

If the same data plays both roles, we risk the well-known
**development-set contamination** trap of machine learning: parameters
end up subtly tuned to one dataset, and "replication" on another
dataset is not really replication, it is applying an already-fitted
pipeline.

## Decision

ID1 is the **development patient**. Methodology, code, visual style,
and analysis flow are developed and iterated on ID1's 13 seizures.
ID2 and ID3 are **replication patients**: once a method is finalized,
it is applied to them without further tuning.

**Critical additional constraint (to avoid contamination):**
Parameter values are chosen from **patient-agnostic principles** —
literature-standard defaults, physiologically justified bounds, or
theoretical considerations — **not** by inspecting ID1's data.

Examples of allowed parameter choices:
- "Minimum ictal duration = 60 s (matches Bastos & Schoffelen 2016
  recommendation for spectral connectivity estimation)."
- "Frequency band = 4-30 Hz (covers theta, alpha, beta — canonical
  epilepsy analysis bands)."
- "Connectivity metric = Pearson correlation of band-limited signals
  (simplest robust choice, comparable across sites)."

Examples of NOT allowed parameter choices:
- "We tried three thresholds on ID1 and 0.3 worked best."
- "The window length was chosen because it maximizes the effect in ID1."
- "We dropped channels below a variance threshold picked by looking
  at ID1's histograms."

The distinction is: the parameter must be defensible **before seeing
the data**. If we cannot write down a rationale that does not depend
on inspecting ID1, we cannot use that parameter.

## Alternatives considered

| Option | Why rejected |
|--------|--------------|
| Cross-validation across the three patients | Only three patients — not enough for meaningful cross-validation folds. |
| Random split of ID1's seizures into dev and test subsets | Possible in principle but wastes half of ID1's statistical power on hold-out data. Since ID2/ID3 already provide external replication, splitting ID1 is redundant. |
| Develop on ID2 and test on ID1 | ID1 has by far the most seizures, so it is the most informative development set. Reversing wastes ID1's richness. |
| No formal roles; iterate on all three together | This is the contamination trap. Rejected explicitly. |

## Consequences

**Easier:**
- Clear division of labor across the three patients.
- Every methodological decision now has a "why is this choice
  patient-agnostic?" question attached to it — forces us to justify
  choices explicitly.
- Replication on ID2 and ID3 becomes a genuine test of generalization.

**Harder:**
- Some parameters may look sub-optimal for ID1 because we picked
  them from literature rather than fitting them. We must resist the
  temptation to "improve" them on ID1.
- If ID2 or ID3 results are weak, we cannot rescue them by tuning
  parameters; we must report the weakness honestly.
- Slower development in early Phase 2, because we spend time
  justifying choices we could otherwise have picked by looking.

**Commitments:**
- Every parameter used in Phase 2 code will have an inline comment
  citing its source (literature, ADR, or theoretical reason).
- No parameter value is chosen by inspecting ID1's numerical output.
  If we ever need to (e.g., a truly patient-specific parameter),
  we will document it as an exception in its own ADR.
- ID2 and ID3 analyses will be run only after the ID1 pipeline is
  frozen. Any change to the pipeline after ID2/ID3 have been touched
  requires re-running everything and noting the revision.

## Open questions to revisit

- What counts as "inspecting ID1's data"? Looking at a raw EEG trace
  is fine (we already did in micro-step 10); looking at a distribution
  of a metric across ID1's seizures and picking a threshold is not.
  The line will need concrete examples as they arise.
- If a parameter genuinely has no patient-agnostic justification and
  must be tuned, we will document that as a limitation and use
  ID1's median as a defensible tuning rule, rather than
  cherry-picking.