# ADR 0002 — Quantum-walk conventions and the hub/periphery open question

**Status:** accepted (with one open question flagged for Phase 2)

**Date:** 2026-07-19

**Author:** Hatice Nalçacı

---

## Context

Before applying continuous-time quantum walks (CTQW) to real brain
networks in Phase 2, we built and validated the core walk machinery on
a small synthetic Watts-Strogatz graph (n=20, k=4, p=0.3, seed=42).
This was done step by step in the `micro_step_*` exploration scripts
and then consolidated into the tested `quantum_walks` module.

Two conventions had to be fixed, and one unexpected observation needs
to be recorded so we do not forget it.

## Decision

1. **Hamiltonian convention.** The CTQW Hamiltonian defaults to the
   adjacency matrix, `H = A` (Mulken & Blumen 2011). The Laplacian
   convention `H = L` (Childs 2010) is available via a parameter. Both
   are kept because the literature uses both and we may need to compare.

2. **Spreading metric.** We use the participation ratio
   `PR = 1 / sum(p_i^2)` to quantify how delocalized a distribution is
   (PR=1 fully localized, PR=n fully uniform).

3. **Return probability** is the probability mass at the start node as
   a function of time. It is our primary signature distinguishing
   quantum (oscillatory) from classical (monotone) dynamics.

These were validated by 14 unit tests encoding the relevant physical
laws: probability conservation, correct t=0 localization, classical
convergence to uniform, quantum non-convergence to uniform, oscillatory
vs monotone return probability, and PR bounds.

## Validated observations (synthetic graph)

- **Probability is conserved** to 6+ decimal places for the quantum
  walk at all tested times (t up to 100). Confirms unitarity is
  correctly implemented.
- **Classical diffusion converges to the uniform distribution**; the
  quantum walk does not. Confirmed numerically and by test.
- **Quantum return probability oscillates** (rebounds to 0.4-0.7 at
  later times after the initial drop from 1.0), a direct visual and
  numerical signature of interference. Classical return probability
  does not rebound.

## Open question flagged for Phase 2

When comparing a quantum walk started at a **hub** node (high degree)
versus a **peripheral** node (low degree), the *instantaneous*
participation ratio did NOT show a clean separation: both rise quickly
to a similar band (~8-14) and then oscillate strongly, overlapping.

This is unexpected. Two competing explanations:

1. **Substantive:** in a quantum walk the starting node's degree may
   matter less for asymptotic spreading than it does classically,
   because interference reorganizes the distribution. If real, this is
   an interesting finding.

2. **Methodological artifact:** the instantaneous PR oscillates so
   strongly that any hub/periphery difference is buried in the
   fluctuations. The fix would be to compare the **time-averaged**
   distribution / PR rather than instantaneous values.

We implemented `time_averaged_distribution` precisely to test
explanation (2) in Phase 2. **Do not draw any conclusion about
hub vs periphery from the instantaneous curves alone.** This must be
resolved with time-averaged quantities on the real connectivity
matrices, not the toy graph.

## Consequences

- The `quantum_walks` module is now the validated foundation for
  Phase 2. Real brain connectivity matrices will be fed into the same
  functions; only the input graph changes.
- The hub/periphery question is the first concrete analysis to run on
  real data, using time-averaged metrics.
- We have NOT yet shown any advantage of CTQW over classical walks for
  any seizure-relevant question. Everything so far is method-building
  and replication of known CTQW properties. The scientific contribution,
  if any, still lies entirely ahead in Phase 2.
