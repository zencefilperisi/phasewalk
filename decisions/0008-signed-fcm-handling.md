# ADR 0008 — Handling negative values and the diagonal in the FCM-as-graph

**Status:** accepted

**Date:** 2026-09-18

**Author:** Hatice Nalçacı

---

## Context

ADR 0007 fixed how we *compute* the functional connectivity matrix
(FCM) from ictal iEEG. But that ADR does not say how the resulting
matrix is *used as a graph*. In practice, taking a Pearson-correlation
FCM and feeding it directly into the walk code raises two concrete
problems:

**Problem 1 — Negative entries.** Pearson correlations lie in [-1, 1].
Real ictal FCMs contain both positive and negative off-diagonal
entries (in ID1/Sz1 the observed range was roughly [-0.51, +0.77]).
The two walk models we use handle this differently:

- **Continuous-time quantum walk (CTQW):** the Hamiltonian must be
  Hermitian; a real symmetric matrix satisfies this regardless of
  sign. So `expm(-i H t)` is well-defined for signed H, evolution
  stays unitary, and probability is conserved. The walk simply
  responds to negative "hopping amplitudes" with different
  interference patterns.
- **Classical diffusion:** driven by the graph Laplacian `L = D - A`.
  If A contains negative entries, `D` (degrees defined as row sums)
  is no longer non-negative, `L` loses the properties that make
  `expm(-L t) p_0` a valid probability distribution, and outputs
  can become negative. Classical diffusion is not defined on a
  signed graph in the usual sense.

Because Phase 2C's central question is *classical vs quantum on the
same graph*, we need a graph on which both models are defined.

**Problem 2 — The diagonal.** The FCM has 1s on the diagonal by
construction (a channel is perfectly correlated with itself). When we
use the FCM as a graph Hamiltonian, diagonal entries are on-site
energies (the "cost of staying at a node"). A uniform diagonal of 1
adds a constant offset to every eigenvalue, changing global phase but
not the observable dynamics. It is nevertheless standard practice in
network-analysis pipelines to zero the diagonal before treating the
matrix as an adjacency, because we want the graph to describe
between-node relations, not self-loops.

## Literature note

The iEEG / brain-connectivity literature does not offer a single
standard for negative values (see the searches recorded on 2026-09-16):

1. **Discard negatives** (most common): keep only the positive part
   of the correlation matrix, zero the negatives. Used by van den
   Heuvel and many others.
2. **Analyse positive and negative networks separately** as two graphs
   with different structural properties.
3. **Keep the signed graph** and adapt the analysis to signed weights.
   This is a newer direction (see e.g. the "adaptive signed random
   walks" preprint, 2023).

On the quantum-walk side, CTQW on signed and complex-weighted graphs
is mathematically defined and has been studied (Kirkland & von Bommel;
Coutinho et al.; the "generalized Laplacian" line of work, 2025), but
it is not yet an established standard tool in applied neuroimaging.

## Decision

Two sub-decisions.

### 1. Diagonal: always zero it before the matrix is used as a graph

Before passing an FCM to any graph function (walk code, degree
statistics, spectral analysis), set `np.fill_diagonal(fcm, 0)`. This
is not a tunable parameter; it is a definitional step encoded inside
the FCM-to-graph conversion.

Rationale:
- Standard practice throughout network analysis on correlation matrices.
- A uniform diagonal contributes only a global phase to CTQW and does
  not affect classical diffusion structure; removing it does not
  change observable dynamics but does keep the graph description clean.
- Simpler code: the FCM-to-graph conversion is a single, deterministic
  transformation.

### 2. Negative values: two-layer approach (positive-only primary, signed robustness check)

**Primary analysis (Layer P):**
Set all negative FCM entries to zero. The resulting matrix has entries
in [0, 1], is symmetric, and has zero diagonal (from sub-decision 1).
Both classical diffusion and CTQW are well-defined and directly
comparable on this graph.

This is Option A from the design discussion — the most common choice
in the literature and the one that lets us keep the *classical vs
quantum* comparison meaningful.

**Robustness check (Layer S), deferred to a later micro-step:**
Repeat the CTQW analysis with the signed FCM (diagonal zeroed but
negatives kept). Compare node rankings and spread metrics against
Layer P.

- **If Layer P and Layer S agree** (same nodes central, same qualitative
  spread pattern): discarding negatives did not affect the finding.
  Report Layer P as primary.
- **If they disagree**: negatives carry information that changes the
  interpretation. Report both, and treat the disagreement as a
  finding worth discussing.

Layer S applies to **CTQW only**, because classical diffusion is not
well-defined on a signed graph. This asymmetry is itself a fact worth
noting in any paper we write: quantum dynamics extends naturally to
signed connectivity in a way classical diffusion does not.

Layer S is deferred (not implemented in the first pass) because the
primary comparison (classical vs quantum) requires a graph on which
both are defined; running Layer S first would leave the central
comparison uncomputed. Layer S enters at micro-step 15 or 16, once
Layer P results are in hand and there is something to compare against.

## Alternatives considered

| Option | Why rejected (as primary) |
|--------|---------------------------|
| **Absolute value** `|FCM|` (Option B in the design discussion) | Preserves magnitudes but collapses "in-phase" and "anti-phase" pairs into one, which misrepresents the underlying biology. If we ever want a robustness check on magnitude only, this is the natural third layer. |
| **Signed graph as primary** (Option C in the design discussion) | Most honest but breaks classical diffusion, which the central *classical vs quantum* comparison depends on. Also loses the ADR-0006 "literature-standard" justification. Kept as a deferred robustness check (Layer S) rather than dropped. |
| **Threshold negatives at some level** (keep only strong negatives) | Reintroduces the free-parameter problem ADR 0007 explicitly avoided. |
| **Analyse positive and negative networks as two separate graphs** | Doubles the analysis and forces us to invent a way to combine or contrast results across the two graphs — a research question in itself, not a methodological detail. Deferred indefinitely. |

## Consequences

**Easier:**
- Layer P is a clean, well-defined graph on which both classical
  diffusion and CTQW give interpretable outputs.
- The central comparison of Phase 2C (classical vs quantum on the
  same graph) is now well-posed.
- Every step from FCM to graph is a single deterministic function
  with no hidden parameters.

**Harder:**
- Layer P discards information (about half of the off-diagonal
  entries in a typical Pearson FCM, in expectation). A reviewer can
  and will ask why. The answer must point to this ADR and to Layer S.
- Layer S is a commitment: we cannot ship a paper based only on
  Layer P without also running Layer S, because we said we would.
  Skipping it silently would violate the ADR.
- Layer S introduces an asymmetry — CTQW gets a robustness check on
  signed graphs, classical does not (because it cannot). We must
  present this honestly rather than gloss over it.

**Commitments:**
- Any FCM-to-graph function in the codebase applies both steps
  (`fill_diagonal(0)` + `clip negatives to 0` for Layer P), and
  exposes the layer as an explicit argument (`layer="positive"` or
  `layer="signed"`).
- Every Phase 2C finding derived from Layer P must be accompanied by
  its Layer S counterpart (for CTQW), or by an explicit note that
  Layer S has not yet been computed for this specific finding.
- The absolute-value option (Layer M for "magnitude") is available
  as a possible third robustness check if Layer P and Layer S
  disagree in ways we cannot resolve. Not committed at this point.

## Open questions to revisit

- If Layer P and Layer S give substantially different node rankings,
  we need a principled way to report both — probably a per-node
  concordance metric (Spearman rank correlation between the two
  rankings) plus a discussion of the discrepancies.
- The magnitude option (Layer M = |FCM|) may be worth adding as a
  third layer if Layers P and S give conflicting messages, to
  disambiguate "sign matters" from "magnitude matters".
- Whether the signed-graph CTQW literature has settled on a preferred
  spreading metric that differs from the participation ratio we use
  now — this needs checking before Layer S is implemented.