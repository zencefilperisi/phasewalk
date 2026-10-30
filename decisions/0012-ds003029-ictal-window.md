# ADR 0012 — ds003029 ictal window and seizure-onset definition

**Date:** 2026-10-30
**Status:** Accepted
**Phase:** Phase 3 (ds003029 clinical validation)

## Context

ds003029 marks seizures with free-text clinician annotations in `events.tsv`.
The step-4 / 4b audits (all 106 ictal runs) established:

- **Onset is reliably markable** (103/106 runs).
- **Offset is NOT reliable** — 20/106 runs (mostly JHH) have no clean end marker;
  the seizure is described by evolution notes that trail off.
- **Durations vary enormously** (8.5 s to 857 s, median 102 s), so a full
  onset-to-offset segment would produce wildly unequal, incomparable FCMs.
- **Marker wording is highly inconsistent** across centers, and JHH distinguishes
  push-button (`(PB SZ)`) from electrographic (`(EEG SZ)` / `Z ELECTROGRAPHIC
  ONSET`) onset, present inconsistently (e.g. jh102/run-03 has both ~6 s apart;
  jh104/run-01 has only push-button).

## Decision

**1. Onset anchor = the EARLIEST seizure-onset marker of any type**, matched
case-insensitively, with offset/administrative markers excluded. Patterns:
`onset`, `sz onset`, `eeg onset`, `eeg sz start`, `SZ EVENT`, `seizure #N`,
`ictal onset`, `definite onset`, and a bare `start`. Uniform across centers
(Option A). Rationale: in a 30 s window the few-second gap between push-button
and electrographic onset is immaterial, and one uniform definition avoids the
run-to-run inconsistency a "prefer-electrographic-else-fall-back" rule would
introduce (that inconsistency would itself be a confound).

**2. Ictal window = [onset, onset + 30 s].** 30 s is carried over from Phase 2C /
SWEC (Bastos & Schoffelen 2016; minimum length for a stable FCM), giving
cross-phase methodological consistency and equal-length, comparable windows
across all subjects regardless of seizure duration.

**3. Offset is used ONLY as a short-seizure cap.** If an offset marker exists and
`offset - onset < 30 s`, the window is capped at offset and the run flagged
`short`. If no offset marker exists, the full 30 s is used — safe because every
NO-OFFSET run is a long JHH seizure whose last annotation is well beyond
onset + 30 s (verified in the audit).

**4. Runs with NO onset marker are dropped** (pt17/run-01, pt17/run-03).
pt11/run-03 is recovered via the bare `start` marker. No subject is lost.

**5. Short seizures (< 30 s of ictal) are flagged** and handled as a duration
sensitivity, consistent with the Phase 2C duration-confounder finding
(micro_step_15b): short windows give noisier FCMs, so pooled results are
reported with and without short-seizure runs.

## Consequences / caveats

- The JHH anchor is the annotated seizure event (push-button or EEG, whichever
  is earliest), which can differ from strict electrographic onset by a few
  seconds; immaterial at the 30 s scale but recorded here.
- A fixed 30 s window deliberately ignores seizure dynamics after the first 30 s;
  this trades late-ictal information for cross-subject comparability — the right
  trade for an onset-zone localization question.
- The loader asserts that ≥ (window length) of signal exists after onset and
  drops any run that fails; given pre-ictal padding places onset ~30–150 s into
  ≥140 s recordings, this effectively never triggers.

## Deferred / future

- A window-length sensitivity (10 s / 20 s) can be run later, mirroring the SWEC
  Layer-1 / Layer-2 comparison, if the 30 s result warrants it.