# ADR 0009 — Reframing the Phase 2C central analysis

**Status:** accepted (supersedes the implicit analysis plan in ADR 0007)

**Date:** 2026-09-16

**Author:** project owner

---

## Context

ADR 0007 fixed the FCM-extraction methodology, and ADR 0008 fixed how
the FCM is used as a graph. Together they implied a Phase 2C central
analysis of the form: "Compute spread (time-averaged participation
ratio) per electrode for classical diffusion and CTQW on the Layer P
graph, then compare the two rankings."

Micro-step 14 implemented that analysis. On ID1/Sz1 the result showed
a clear methodological problem, not a scientific one:

    Classical PR: min = 46.86, max = 46.96, mean = 46.93   (n = 47)
    Quantum   PR: min = 7.88,  max = 41.58, mean = 25.30
    Spearman rho(classical, quantum) = +0.195  (p = 0.189)

The classical PR is **saturated at the ceiling** across every start
node. Classical diffusion on a connected graph converges to its
stationary distribution (approximately uniform for the Layer P
graph), so the time-averaged distribution over a long window
(t_grid = 0.1 to 100 was used) is essentially uniform regardless of
start node, giving PR ~= n for every electrode.

CTQW does *not* converge to a stationary distribution (verified in
micro-step 5 and in the earlier synthetic-graph work), so its
time-averaged PR retains information about the start node, and does
vary across electrodes (7.9 to 41.6).

The Spearman correlation between the two rankings is therefore
statistically meaningless: a rank correlation against an all-tied
vector is dominated by tie-breaking noise, not by any relationship.
No conclusion about "classical vs quantum agreement" can be drawn
from this number. The comparison as originally framed is
**not well-posed for our time grid**.

## Options considered

Two ways forward were on the table.

**Option 1 — Shorten the time grid.**
Restrict the time average to short times where classical diffusion is
still spreading (e.g. t in [0.1, 5]) so that classical PR also varies
across start nodes and a genuine ranking comparison becomes possible.

This was rejected because choosing "which short window" would either
have to be tuned to make ID1's classical output look interesting —
a direct violation of ADR 0006's parameter-independence rule — or
justified from a literature standard we do not yet have for this
specific setting.

**Option 2 — Reframe the question.**
Accept that classical diffusion carries no per-electrode information
in the long-time limit and use it as a **null baseline** rather than
as a peer comparator. Study the CTQW distribution across electrodes
in its own right, asking which electrodes have low time-averaged
spread (are "localized" by quantum interference) and how that pattern
relates to graph-topological properties (degree, and later betweenness,
clustering, eigenvector).

## Decision

Adopt Option 2. The Phase 2C central analysis is reframed as:

1. **Primary output:** the distribution of time-averaged quantum
   participation ratio across electrodes for a given seizure on its
   Layer P graph.
2. **Interpretation of extremes:** electrodes with the lowest quantum
   PR are those most localized by interference on this ictal network;
   electrodes with the highest quantum PR are those from which the
   walk spreads most broadly. These are the electrodes of primary
   interest, ranked rather than thresholded.
3. **Baseline reference:** classical diffusion is included as a
   sanity check that the graph is connected (classical PR should be
   near n for every start node, as expected). It is not treated as a
   peer comparator; the classical vs quantum comparison via Spearman
   rank correlation is dropped.
4. **Follow-up structural question:** relate the per-electrode
   quantum PR to graph-topological node properties (degree first,
   later betweenness / eigenvector / clustering — using the same
   pipeline as micro-steps 6-9 on synthetic graphs).

Rankings only, no thresholds. Statements are of the form "electrode X
has the k-th lowest / highest quantum PR", never "electrode X has
PR below threshold T".

## Implications for other ADRs

- **ADR 0007** stays in force as written; only the downstream
  interpretation of the Layer 1 / Layer 2 outputs is reframed here.
  Its "Open questions to revisit" section is updated with a pointer:
  "The classical-vs-quantum rank correlation question is dropped;
  see ADR 0009 for the reframed analysis."
- **ADR 0008** stays in force. Layer P remains the primary graph.
  Layer S (signed graph, CTQW only) is still on the roadmap for a
  later micro-step and will be compared against Layer P using the
  reframed analysis of this ADR (rank agreement of quantum PR
  between the two layers).
- **ADR 0006** parameter-independence rule is reinforced: this
  reframing was the direct consequence of refusing to tune the time
  grid to ID1's classical output.

## Consequences

**Easier:**
- The analysis is now well-posed: it asks something CTQW can answer
  and does not depend on classical diffusion having per-node variation.
- Rankings are dimensionless and permit clean cross-seizure and
  cross-layer comparisons even where absolute PR values are not
  directly comparable.
- The parameter-independence rule (ADR 0006) is satisfied by
  construction: no tunable threshold, no ID1-tuned window choice.

**Harder:**
- The paper narrative loses its most natural framing ("classical vs
  quantum on the same graph — here is what changes"). The new
  framing ("what does CTQW localize on the ictal network") is
  more subtle and needs a clearer motivation section.
- Interpreting "low PR = interesting electrode" requires a downstream
  argument for why this quantity has clinical meaning. That argument
  is not made in this ADR; it is a future step (relating quantum PR
  to seizure onset zone marks, if such marks are available for the
  SWEC-ETHZ patients — currently unclear).

**Commitments:**
- The Phase 2C code path uses the reframed analysis from micro-step
  14b onwards. Micro-step 14 stays in `exploration/` as the record
  of the earlier attempt.
- Any published finding will state explicitly that the classical
  arm is a null baseline and that the primary object of study is the
  CTQW-only quantity.

## Open questions to revisit

- Does the ranking of quantum PR change substantially across
  different seizures of the same patient (ID1's 13 seizures)? If
  yes, the interpretation "electrode X is localized" is
  seizure-specific rather than patient-specific.
- Does the same pattern hold at the Layer S (signed graph) level?
  Deferred to the micro-step that implements Layer S.
- If seizure-onset-zone electrode marks are available for ID1
  (needs checking), does low quantum PR correlate with them? This
  would be the first genuine clinical bridge for the project.