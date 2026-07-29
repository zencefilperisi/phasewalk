# ADR 0004 — Handling seizure duration heterogeneity

**Status:** accepted

**Date:** 2026-05-27

**Author:** project owner (with methodological review)

---

## Context

The bulk scan of ID1, ID2 and ID3 (micro-step 11) revealed that
seizure durations vary enormously, both within and across patients:

| patient | n seizures | shortest (s) | longest (s) |
|---------|------------|--------------|-------------|
| ID1     | 13         | 10           | 252         |
| ID2     | 4          | 96           | 301         |
| ID3     | 2          | 73           | 125         |

In ID1 specifically, four of the thirteen seizures are 20 seconds or
shorter (Sz4, Sz9, Sz10, Sz11) — roughly 30% of that patient's data.

This creates a methodological problem for the planned analyses.
Functional connectivity matrices and quantum-walk statistics estimated
from ~10 seconds of data are substantially noisier than those from
~100 seconds of data. Naively pooling short and long seizures could
mix reliable and unreliable estimates, contaminating any group-level
statistic.

Two obvious extreme approaches:

- **Discard short seizures** below some duration threshold. Clean
  statistically but throws away potentially informative clinical data
  and risks selection bias (short seizures may have a distinct
  pathophysiology).
- **Keep everything unconditionally.** Preserves all data but risks
  duration acting as a hidden confounder.

## Decision

Neither extreme. We adopt a **duration-stratified reporting** policy:

1. **Include all seizures in the primary analysis** — no seizure is
   discarded on the basis of duration alone.
2. **For every reported result, also report it restricted to a
   long-seizure subset** (threshold to be documented per analysis;
   typical choice: seizures with an ictal segment of at least
   60 seconds, matching common functional-connectivity minimums in
   the literature).
3. **If the two versions agree** (same qualitative pattern, quantitative
   values within reasonable noise): duration is not a confounder,
   report the full-data result as primary.
4. **If they disagree**: report both explicitly and treat the
   discrepancy as a finding, not a nuisance. Discrepancy tells us
   short-seizure dynamics differ from long-seizure dynamics — which
   is itself worth writing about.

This is a stricter form of A2 in the design discussion: not just
"keep everything and check", but "keep everything, always check, and
treat disagreement as information".

## Alternatives considered

| Option | Why rejected |
|--------|--------------|
| Hard cut-off at 30 or 60 s | Loses ~30% of ID1's data. Also, the threshold choice would look arbitrary to a reviewer; no clean physiological justification for any specific number. |
| Weight seizures by duration | Sounds principled but hides the short-seizure signal inside a weighted average, making it hard to see whether short seizures behave differently. |
| Analyse short and long seizures separately as the primary analysis | Doubles the reporting burden and reduces power in each subset. Better as a robustness check (which is what this ADR institutes) than as the main design. |

## Consequences

**Easier:**
- No data is thrown away.
- Duration effects become visible instead of hidden.
- The threshold used per analysis is a documented choice, not a
  hidden pipeline parameter.

**Harder:**
- Every result requires a paired subset analysis. Extra work per
  finding.
- Reporting becomes more verbose (two numbers instead of one).
- If short and long seizures disagree systematically, we must
  interpret that, not just report it.

**Commitments:**
- The functional-connectivity code will accept an `min_ictal_duration_s`
  parameter and expose it in every analysis notebook.
- Any published finding will be accompanied by the corresponding
  long-seizure-only replication (or a clear statement of which
  seizures were included).
- The choice of threshold per analysis will be justified inline in
  the analysis code, not buried in a config file.

## Open questions to revisit

- Once we have real connectivity matrices in hand, we may find that
  the useful duration threshold is well below 60 s (e.g., stable
  estimates from 20 s of ictal signal). If so, revise the threshold
  and note it in this ADR.
- If a specific short-seizure pattern turns out to be systematically
  different from long ones, we may need a separate ADR on
  interpretation.