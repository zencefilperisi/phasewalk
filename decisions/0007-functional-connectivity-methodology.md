# ADR 0007 — Functional connectivity extraction methodology (Phase 2C)

**Status:** accepted

**Date:** 2026-08-04

**Author:** Hatice Nalçacı

---

## Context

Phase 2C requires turning each seizure recording (an
electrodes-by-time matrix) into a **functional connectivity matrix
(FCM)** — a symmetric electrodes-by-electrodes matrix whose entry
`(i, j)` measures how strongly channels `i` and `j` co-vary. The FCM
is the graph on which classical and quantum walks will subsequently
be run.

The iEEG functional-connectivity literature offers no single "gold
standard" method. A brief review (May 2026) surfaced several widely
used measures:

- **Pearson correlation** — amplitude-based, simplest, fastest.
- **Coherence** — frequency-domain amplitude coupling.
- **Phase Locking Value (PLV)** — phase-based synchrony.
- **Amplitude Envelope Correlation (AEC)** — envelope of band-limited
  signal correlation.
- **Lagged coherence / Phase Lag Index (PLI)** — insensitive to
  zero-lag interactions (volume conduction robust).

A recent review (Front. Neurol. / IOP J. Neural Eng., 2025) points
out that **the choice of measure materially changes the reported
findings**: coherence and PLV report ictal *hyper*-synchronization,
while lagged coherence and PLI report ictal *desynchronization* of
the epileptogenic zone. This is not a small caveat — it says the
"finding" is partly a property of the method.

Because Phase 2C is exploratory and the goal is to compare quantum
against classical walks on the same graphs (not to argue for one
connectivity measure over another), we make the following minimal,
defensible choices, following ADR 0006's parameter-independence rule.
No parameter below was chosen by inspecting any patient's data.

## Decision

Five sub-decisions, each independently revisable:

### 1. Primary connectivity measure: Pearson correlation

We compute the FCM as the Pearson correlation matrix of the ictal
signal across channels:

```
FCM[i, j] = corr(EEG[i, :], EEG[j, :])   # over the chosen segment
```

Rationale:
- Simplest and fastest measure, minimising the number of hidden
  parameters (window length, taper, band edges) that a reviewer
  could question.
- Widely used as a baseline in the iEEG functional-connectivity
  literature.
- The purpose of Phase 2C is to compare walk models on a graph, not
  to advocate a particular connectivity measure. The simpler the
  graph-construction step, the cleaner that comparison.

### 2. Robustness measure: Phase Locking Value (PLV)

Every analysis that uses the Pearson FCM as its input will be
repeated with a PLV-based FCM. The two results are then compared:

- **Agreement** → the finding is not an artefact of the connectivity
  measure. Report it with confidence.
- **Disagreement** → the finding depends on the measure. Report both
  outcomes explicitly; the disagreement itself is a finding worth
  discussing.

Rationale: Pearson and PLV are the two most-used measures in the
epilepsy iEEG literature, and they disagree in known systematic ways
(amplitude- vs phase-driven). Using them as a pair is the standard
way to check whether an ictal result is method-dependent.

### 3. Frequency band: broadband (as-recorded)

The FCM is computed on the signal as delivered by the loader — no
band-pass filtering beyond the 0.5-150 Hz filter applied by the
SWEC-ETHZ acquisition pipeline.

Rationale:
- Starting broadband minimises the number of methodological choices
  under discussion in the first paper. Frequency-band analyses are
  a natural extension (a future ADR).
- Broadband correlation is the simplest possible input for the graph
  analysis and mirrors what the underlying signal actually contains.
- Restricting to a single band this early risks tuning the band
  choice to whichever produces the "cleanest" figure in ID1 — the
  exact ADR-0006 trap.

Frequency-resolved (delta / theta / alpha / beta / gamma / high-gamma)
analyses are deferred to a later ADR once the broadband pipeline is
stable and interpretable.

### 4. Segmentation: two-layer analysis

We produce **two** FCMs per seizure and compare their downstream
outputs:

- **Layer 1 (primary):** the FCM is computed from the *entire* ictal
  segment, whatever its duration. A 10-second seizure contributes a
  matrix estimated from 10 seconds; a 252-second seizure contributes
  a matrix estimated from 252 seconds.
- **Layer 2 (robustness):** the FCM is computed from the *first 30
  seconds of the ictal segment* only. Seizures whose ictal segment
  is shorter than 30 seconds are excluded from Layer 2 and reported
  as excluded.

The 30-second threshold is taken from the practical minimum used in
the connectivity-estimation literature (see e.g. Bastos & Schoffelen,
Front. Syst. Neurosci. 2016). It was not chosen by inspecting ID1.

Layer 1 asks "what does the seizure as a whole look like as a graph?"
Layer 2 asks two questions at once: "does the seizure-duration
difference drive our results?" and "when we restrict every seizure
to the same clinical phase (onset + early propagation, before any
plateau or termination phase in long seizures), does the pattern
change?" A disagreement between Layer 1 and Layer 2 is itself a
finding — it tells us the effect depends on either duration or on
which phase of the seizure is included.

**Comparison protocol:** When comparing Layer 1 with Layer 2, the
comparison is restricted to the subset of seizures that qualify for
Layer 2 (ictal duration >= 30 s). Layer 1 results on this subset are
compared against Layer 2 results on the same subset. This isolates
the effect of "full ictal vs first 30 s" from the effect of "different
sample of seizures". Short seizures are NOT excluded from the primary
Layer 1 analysis — they contribute to the full-cohort Layer 1 result,
which is reported alongside the matched subset comparison. This
satisfies ADR 0004's rule that no seizure is discarded on the basis
of duration alone.

This structure is a direct application of ADR 0004's duration-stratified
reporting policy: keep all seizures in the primary analysis, always
report a matched-duration replication, and treat any disagreement as
information, not as a nuisance to hide.

### 5. Thresholding: none — the FCM stays fully weighted

The FCM is passed to the walk code as a fully connected weighted
matrix. No edges are thresholded out.

Rationale:
- Any threshold is a free parameter that can be silently tuned to
  ID1's data. Refusing to threshold removes the temptation entirely.
- The quantum-walk machinery in `quantum_walks/walks.py` already
  handles weighted matrices — its Hamiltonian is the input matrix
  directly, no binarization required.
- Classical diffusion likewise operates on weighted Laplacians without
  requiring a threshold.
- If thresholding turns out to be necessary for a specific downstream
  analysis, it will be introduced by its own ADR with a documented
  patient-agnostic rationale.

## Alternatives considered (per sub-decision)

| Sub-decision | Alternative | Why rejected (now) |
|--------------|-------------|--------------------|
| 1 (measure)  | AEC as primary | Similar simplicity to Pearson but less standard as a baseline. Would still need PLV as robustness check. Adds nothing over Pearson at this stage. |
| 1 (measure)  | Lagged coherence / PLI | Solves the volume-conduction concern, but SWEC-ETHZ is intracranial (volume conduction less severe than scalp EEG), and lagged measures introduce their own parameter (which lags). Deferred. |
| 2 (robustness) | Coherence instead of PLV | Coherence and Pearson share the amplitude-driven character; the pair is less informative than a Pearson+PLV pair, which spans amplitude vs phase. |
| 3 (band)     | Six-band decomposition from the start | Multiplies analyses by six and forces early choices we cannot yet justify. Deferred to a follow-up ADR. |
| 4 (segments) | Duration-matched: use the shortest seizure's length across all seizures | The shortest ictal segment in ID1 is 10 s; restricting all seizures to 10 s throws away 96% of Sz5's ictal signal. Layer 2's 30-second threshold discards short seizures instead, which is less destructive to the informative long ones. Also, "shortest = 10 s" would be an ID1-tuned choice; "30 s from literature" is not. |
| 4 (segments) | Sliding windows within each seizure | Introduces window-length and overlap parameters, which are another ID1-tuning risk. Deferred to a later ADR once the static two-layer analysis is understood. |
| 5 (threshold) | Percentile threshold (e.g. keep top 20% of edges) | The percentile itself is a free parameter and is a well-known route to unconsciously tuning a result to one dataset. |

## Consequences

**Easier:**
- One line of code per seizure produces its FCM.
- The walk code can be applied unchanged; the graph is just no longer
  synthetic.
- Two-layer analysis is a small addition (compute two FCMs, run the
  downstream code twice) and gives us built-in robustness reporting.
- Every parameter has a documented, patient-agnostic origin.

**Harder:**
- Pearson is not the most sophisticated measure. Reviewers may ask
  "why not X?" — we must be ready to point to this ADR and the
  robustness plan.
- Every finding requires the paired PLV computation. Roughly doubles
  compute time for connectivity estimation (still very cheap).
- If Layer 1 and Layer 2 disagree systematically, we owe the paper a
  discussion of why, not a quiet averaging.

**Commitments:**
- The Phase 2C connectivity code will expose `measure` (`"pearson"` or
  `"plv"`) and `layer` (`"full"` or `"first30s"`) as explicit
  parameters; no hidden defaults.
- Every reported finding will be accompanied by its cross-measure
  (Pearson vs PLV) and cross-layer (full vs first-30s) counterpart,
  or an explicit statement of why the counterpart was not computed.
- Choices in this ADR will be revisited if the paired analyses
  reveal systematic method- or segment-dependence that the current
  minimalism cannot explain.

## Open questions to revisit

- If Pearson and PLV disagree on a large fraction of seizures, we
  will need to decide whether to add a third measure (e.g. AEC or
  PLI) as a tiebreaker, and how to report the three-way comparison.
- If Layer 1 and Layer 2 disagree, the natural next step is a proper
  time-resolved analysis (sliding windows within each seizure). This
  gets its own ADR at that point.
- Frequency-band decomposition is deferred; once the broadband
  pipeline is stable, a follow-up ADR will introduce it.
- Whether to add anatomy-informed constraints (e.g., grouping
  electrodes by lobe) remains open and depends on whether SWEC-ETHZ
  provides usable electrode coordinates for our patients.
- If Layer 1 and Layer 2 produce systematically different results, we
  will need a procedure to distinguish (a) duration acting as a
  confounder, (b) different seizure phases producing different graphs,
  and (c) both. This will require a sliding-window analysis, which
  is deferred to a follow-up ADR (candidate 0008), to be written when
  we actually observe the disagreement — not before.
- The classical-vs-quantum rank correlation question raised implicitly
  here has been dropped; see ADR 0009 for the reframed analysis.

  # Note added 2026-09-30 (from micro-steps 15, 16, 15b diagnostic)

**Status of the ADR itself:** unchanged. This note extends the ADR
with empirical findings; it does not revise any of the original
decisions.

## Empirical background

Cross-seizure consistency of CTQW-based per-electrode participation
ratio rankings was tested on ID1 (per ADR 0006, ID1 is the development
patient). Three measurements:

| Setting                                                | Seizures | Mean pairwise Spearman rho |
|--------------------------------------------------------|---------:|---------------------------:|
| Layer 1, all seizures (micro-step 15)                  |       13 |                    +0.310  |
| Layer 1, subset ictal >= 30 s (micro-step 16)          |        9 |                    +0.410  |
| Layer 2, first 30 s of same subset (micro-step 16)     |        9 |                    +0.241  |

The +0.100 shift from all-seizures Layer 1 to long-only subset Layer 1
raised a legitimate concern: is that shift a duration effect, or just
an artefact of dropping n? Micro-step 15b decomposed the 78 pairs into
those involving at least one short seizure (Group A) and those among
long seizures only (Group B):

| Group                                | Pairs | Mean rho |
|--------------------------------------|------:|---------:|
| A (>= 1 short seizure, <30 s)        |    42 |   +0.225 |
| B (both seizures long, >=30 s)       |    36 |   +0.410 |

Mann-Whitney U p < 0.0001. Group B mean matches the Layer 1 subset
mean exactly — sanity check passed.

## Two findings

**Finding 1 — short seizures (<30 s) genuinely add noise to the FCM,
not merely reduce n.**
The individual pair rhos are systematically lower when a short seizure
is involved. This is a real signal, not a sampling artefact. The 30-s
minimum window from Bastos & Schoffelen (2016) is validated empirically
for this dataset: with fewer than about 15,000 samples per channel, the
Pearson FCM is too variable to give a stable per-electrode CTQW
ranking.

**Finding 2 — Layer 2 (first 30 s) is LESS consistent than Layer 1
(full ictal) on the same seizures.**
Contrary to the assumption that a fixed window would reduce noise and
raise cross-seizure consistency, restricting to the first 30 s DROPS
mean pairwise Spearman by -0.169 (from +0.410 to +0.241). The best
reading of this is biological:

- The seizure **onset** phase is per-seizure heterogeneous — different
  onset zones, different propagation patterns.
- The **full ictal** recording averages over onset + plateau +
  resolution. The shared plateau structure carries a more stable
  ranking signature that shows up across seizures.

## What this changes in interpretation

The original motivation for Layer 2 (ADR 0007 sub-decision 4) was
"robustness check for FCM sampling bias." That framing assumed a
fixed-window FCM would give more comparable rankings across seizures.
The data shows this assumption is wrong for CTQW ranking consistency
in ID1.

Layer 2 is **retained** but its interpretive role is reframed:

- **Layer 1 (full ictal, primary):** characterises the stereotyped
  ictal network signature of the patient. Higher cross-seizure
  consistency. This is what CTQW ranking is best at.
- **Layer 2 (first 30 s, secondary):** now framed as characterising
  the seizure **onset** network — a different biological question, not
  a check on Layer 1. Lower cross-seizure consistency by design,
  because onset is heterogeneous.

The two layers therefore answer different questions:
- "What is this patient's ictal-network signature?" → Layer 1
- "What does the onset network look like across seizures?" → Layer 2

## What does NOT change

- ADR 0007 sub-decisions 1 (Pearson primary), 2 (PLV robustness),
  3 (broadband 0.5-150 Hz), and 5 (no thresholding) are unchanged.
- The Comparison Protocol (when comparing Layer 1 with Layer 2, use
  the matched Layer-2-qualifying subset for both) is unchanged.
- ADR 0009's primary quantity (Layer 1 CTQW per-electrode PR ranking)
  is unchanged.
- ADR 0004's rule that no seizure is discarded on the basis of duration
  alone in the PRIMARY analysis is unchanged: short seizures still
  contribute to the full-cohort Layer 1 result, and both the
  full-cohort (+0.310) and long-only-subset (+0.410) numbers are
  reported alongside each other.

## Resolved open question

The earlier "Open questions to revisit" item — "If Layer 1 and Layer 2
produce systematically different results, we will need a procedure to
distinguish (a) duration acting as a confounder, (b) different seizure
phases producing different graphs, (c) both" — is now resolved:

- (a) Duration IS a confounder for short seizures — validated by
  micro-step 15b (Mann-Whitney U p < 0.0001).
- (b) Different seizure phases DO produce different graphs — validated
  by micro-step 16 (Layer 1 vs Layer 2 mean per-seizure rho = +0.478,
  and within-cohort Layer 2 mean is -0.169 lower than Layer 1).
- (c) Both effects operate; they were separated by the sequence
  15 → 16 → 15b.

## New open question

Does the "onset heterogeneous, full-ictal stereotyped" pattern hold
in ID2 and ID3? If yes, it is a patient-general property of ictal
networks and part of the methodological contribution. If no, ID1 has
an unusual dynamic and the framing needs revision. To be tested in
the ID2 replication (planned micro-step 17).