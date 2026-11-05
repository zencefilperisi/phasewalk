# ADR 0014 — Phase 4: peripheral-SOZ hypothesis (pre-registration)

**Date:** 2026-11-01
**Status:** Accepted (pre-registered BEFORE any held-out dataset is examined)
**Phase:** Phase 4 (classical functional-network localization of the SOZ)

## Why this is written now

Phase 3 tested, and rejected, a pre-registered quantum-walk hypothesis for
seizure-onset-zone (SOZ) localization (ADR 0013; CTQW participation ratio does
not localize the SOZ above chance). In the course of that analysis, an
**exploratory, post-hoc** observation emerged: in the ds003029 discovery data,
SOZ electrodes tend to have **lower** weighted degree and eigenvector centrality
(i.e. a more peripheral, weakly-connected position in the ictal functional
network) than non-SOZ electrodes — the reverse of the pre-registered "epileptic
hub" direction, with per-subject AUC ≈ 0.78 in Engel-I patients.

Because that direction was chosen **after** seeing the hub direction fail, it is
exploratory and cannot be claimed from ds003029. This ADR locks the hypothesis,
direction, metric, population and decision rule **before** any held-out dataset
is examined, so that a held-out test is genuinely confirmatory. ds003029 is the
discovery set and is **not** eligible as the confirmatory dataset.

## What the discovery data already shows (exploratory, for context only)

- Effect present in Engel-I (n = 17): median AUC(−degree → SOZ) = 0.777,
  one-sided Wilcoxon p = 1e-4; eigenvector similar.
- **Not an electrode-type artefact:** all 17 Engel-I subjects are pure-ECoG, so
  there is no depth-vs-surface confound within subjects by construction
  (`ds003029_confound_electrode_type.py`).
- One clear exception (pt3, SOZ fraction 0.38) runs the other way, suggesting
  the effect may be specific to focal (minority) SOZs.
- Plausible mechanism: the ictal "core" that generates the seizure sustains its
  own rhythm and decouples from the propagating activity, lowering its ictal
  functional degree (consistent with the isolated-ictal-core literature). This
  is an ICTAL mechanism; it is not expected to hold interictally.

## Pre-registered confirmatory hypothesis

**H1.** In the ictal functional connectivity network, SOZ electrodes have
**lower weighted degree** than non-SOZ electrodes, in seizure-free (Engel I)
patients.

- **Score direction (LOCKED):** low degree → SOZ (score = −degree). This is the
  single pre-committed direction; we will **not** read whichever direction
  scores higher.
- **Primary metric:** per-subject AUC of (−degree) predicting SOZ membership over
  the subject's good iEEG channels (prevalence-robust).
- **Primary population:** Engel I (seizure-free), where the clinician SOZ is the
  most trustworthy reference (as in ADR 0013).
- **Primary test:** one-sample, one-sided Wilcoxon signed-rank on {AUC_i − 0.5},
  α = 0.05; effect reported as median AUC (IQR).

## Confirmatory dataset

The confirmatory test is run on the **first held-out iEEG dataset with (a) ictal
recordings, (b) electrode-level SOZ annotations and (c) surgical outcome, that is
not ds003029.** Target: HUP ictal iEEG via the IEEG Portal (ieeg.org), the source
used by the neural-fragility studies. The dataset is named and its
inclusion/windowing fixed in a follow-up ADR **before** any AUC is computed on
it, mirroring ADR 0011/0012. Pipeline (ictal window → Pearson FCM → Layer P →
weighted degree) is carried over unchanged from Phases 2C/3.

If no ictal + SOZ held-out dataset can be obtained with the available (solo,
pre-mentor) access, the confirmatory test is deferred; the interictal
generalization below is then run and reported **as a different question**, never
as confirmation of H1.

## Secondary / pre-specified analyses

1. **Eigenvector centrality** in place of degree (same direction).
2. **All-outcomes** population (not only Engel I).
3. **SOZ-fraction moderator:** AUC vs SOZ fraction; pre-specified expectation that
   the effect is stronger for focal (small) SOZs and may vanish or reverse for
   diffuse SOZs (the pt3 pattern). Reported, never used to filter subjects.
4. **Electrode-type stratification** (ECoG vs SEEG) on the held-out data, since
   it may contain SEEG unlike the all-ECoG discovery cohort.
5. **Label-permutation null** (as ADR 0013).
6. **Interictal generalization (exploratory, different question):** test whether
   the peripheral-SOZ structure also appears in interictal functional networks
   (e.g. HUP / ds004100). A null here does not bear on H1, because the proposed
   mechanism is ictal-specific.

## Confounds explicitly acknowledged

- **Grid-edge / spatial degree:** edge electrodes have fewer close neighbours and
  may show lower functional degree. If electrode coordinates are available on the
  held-out data, this is controlled (e.g. adjust for or stratify by edge
  position); otherwise it is flagged as a residual confound.
- **Signal-quality low degree:** a channel weakly correlated with others could be
  noisy rather than an isolated core; only good-status channels are used, and the
  interictal-vs-ictal contrast (iii, below) distinguishes an ictal mechanism from
  a static channel property.

## Mechanism check on the discovery data (exploratory, not confirmatory)

Separately, on ds003029 (discovery), test whether the SOZ's low degree is
**ictal-specific**: compare each SOZ electrode's functional degree in the ictal
window against a pre-ictal baseline window from the same recording. If SOZ degree
drops specifically during the ictal window, that supports the decoupling-core
mechanism; if it is already low interictally, the effect is a static network
property. This is exploratory evidence about the mechanism, not a test of H1.

## Decision criteria (pre-committed)

- **Support for H1** = on the held-out data, Engel-I median AUC meaningfully > 0.5
  AND one-sided Wilcoxon p < 0.05 AND the effect survives the electrode-type and
  (where testable) spatial-edge checks.
- **No support** = AUC ≈ 0.5 or non-significant on held-out data → reported as a
  negative result; the ds003029 lead is then recorded as a discovery that did not
  replicate.
- The primary metric (AUC), direction (low degree → SOZ), population (Engel I) and
  test (one-sided Wilcoxon) will not be changed after seeing held-out results.

## What is explicitly NOT claimed

- No causal or clinical-utility claim.
- The ds003029 result is exploratory; only the held-out test is confirmatory.
- This is a **classical** network hypothesis; it does not revive the quantum-walk
  thesis, which Phase 3 rejected.