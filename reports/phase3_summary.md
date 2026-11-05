# Phase 3 — Summary (ds003029 clinical validation)

**Date:** 2026-11-05
**Scope:** The pre-registered CTQW test for seizure-onset-zone (SOZ)
localization on the Epilepsy-iEEG-Multicenter dataset (ds003029), its
negative result, and the exploratory classical lead that followed, with
its confound and mechanism characterization.
**Status:** The pre-registered quantum-walk hypothesis is **rejected**. An
exploratory classical hypothesis (peripheral SOZ) emerged, is robust to the
main confound, has a plausible ictal mechanism, and is itself pre-registered
(ADR 0014) for held-out confirmation. No confirmed clinical finding yet.

This document is the honest written record of Phase 3, and the skeleton the
Phase 4 held-out result will slot into.

---

## 1. Question (pre-registered)

Does the continuous-time quantum walk (CTQW) per-electrode participation
ratio (PR) identify the clinician-marked SOZ — specifically, are low-PR
(localized) electrodes the SOZ? The full analysis plan (metric, population,
statistic, decision rule) was fixed before any cohort result was seen
(ADR 0013).

## 2. Data and pipeline

- **Dataset:** ds003029 (Li et al. 2021, neural fragility), multi-centre
  ECoG/SEEG, BrainVision, with clinician SOZ annotations and surgical
  outcome (Engel). Inclusion/channel policy in ADR 0011; ictal window
  [onset, onset+30 s] in ADR 0012.
- **Cohort (frozen before analysis):** 29 usable subjects, **17 Engel I**
  (seizure-free; the primary population, where the SOZ reference is most
  trustworthy). umf002–005 excluded (signal shorter than the annotated
  seizure, verified against source).
- **Pipeline:** ictal window → Pearson functional connectivity → Layer P
  graph → CTQW time-averaged PR per electrode. Per-subject AUC of
  (−PR) predicting SOZ membership (prevalence-robust). All tested and
  reproducible; the batched CTQW routine is equivalence-tested against the
  per-node implementation.

## 3. Primary result — negative

In the Engel-I group (n = 17):

| quantity | value |
|---|---|
| median AUC(−PR → SOZ) | **0.39** |
| one-sided Wilcoxon (AUC > 0.5) | **p = 0.96** |
| label-permutation null | observed 0.44 vs null 0.50 (below chance) |

CTQW-PR does **not** localize the SOZ above chance; if anything the effect is
in the opposite direction. The result holds for the all-outcomes and Engel > I
populations and under the short-seizure sensitivity. Because the analysis was
pre-registered, this is a credible negative result, not a framing artefact.

**Classical baseline (pre-specified "hub" direction, +centrality):** weighted
degree AUC 0.22, eigenvector 0.24, betweenness 0.47; CTQW 0.39. CTQW is
statistically higher than degree/eigenvector (paired p = 0.003) but all are
below 0.5 — CTQW is "less wrong in the hub direction," not useful. The
pre-registered hub hypothesis for classical centralities is also rejected.

## 4. The exploratory lead — peripheral SOZ (NOT pre-registered here)

That degree/eigenvector score ≈ 0.22–0.24 in the hub direction means the
**reverse** direction (low centrality → SOZ) scores ≈ 0.76–0.78. I.e.:

> **SOZ electrodes occupy a peripheral, low-degree position in the ictal
> functional connectivity network.**

This direction was chosen after the pre-registered directions failed, so it is
exploratory and cannot be claimed from ds003029. Two characterizations were run
to decide whether it is worth a formal held-out test.

### 4.1 Confound check — electrode type (not an artefact)
Depth (SEEG) electrodes can be intrinsically less connected; if the SOZ sits on
depths, "peripheral = SOZ" could be a placement artefact. But **all 17 Engel-I
subjects are pure-ECoG**, so within each subject the AUC compares SOZ vs non-SOZ
among same-type electrodes — no depth-vs-surface confound exists by
construction. The effect in pure-ECoG subjects: median AUC(−degree → SOZ)
= **0.777**, one-sided Wilcoxon **p = 1e-4**. The lead is a genuine
within-ECoG network property, not an electrode-type artefact.

### 4.2 Mechanism check — ictal decoupling on a static baseline
Comparing each electrode's degree in a pre-ictal window ([onset−35 s,
onset−5 s]) vs the ictal window:

| window | median AUC(−degree → SOZ) |
|---|---|
| pre-ictal | 0.667 |
| ictal | 0.754 |

paired Wilcoxon (ictal > pre-ictal) **p = 0.013**. So the SOZ becomes
**significantly more peripheral during the seizure** (consistent with the
isolated-ictal-core idea: the seizure-generating core sustains its own rhythm
and decouples). But the pre-ictal AUC is already 0.67 — a **static peripheral
tendency exists before the seizure too**; the picture is baseline tendency +
ictal strengthening, not a purely ictal-emergent effect. The effect is
heterogeneous across subjects (jh105: pre-ictal 0.33 → ictal 0.80, a dramatic
decoupling; ummc006 and pt11 reverse).

### 4.3 A pre-specified caveat
pt3 (SOZ fraction 0.38, a diffuse SOZ) runs the other way (AUC 0.35),
suggesting the effect may be specific to focal (minority) SOZs — a moderator
carried into the Phase 4 pre-registration.

## 5. Honest status of the lead

The peripheral-SOZ lead is **exploratory**: a post-hoc direction on the
discovery data. The confound and mechanism checks make it worth pursuing — it
is not an electrode-type artefact and it has a plausible, partly ictal-specific
mechanism — but they do **not** confirm it. Confirmation requires a held-out
ictal + SOZ dataset, with direction, metric, population and decision rule fixed
in advance. That pre-registration is ADR 0014.

## 6. Limitations

- The negative CTQW result is solid (pre-registered, multi-centre); the
  classical lead is only a discovery.
- Small N (17 Engel I); all-ECoG, so no evidence yet that the lead generalises
  to SEEG.
- Residual confounds for the lead not testable without electrode coordinates
  (e.g. grid-edge electrodes have fewer neighbours → lower degree).
- The pre-ictal window is close to onset, not true interictal; it suggests but
  does not establish that an interictal validation would succeed.
- Held-out ictal + SOZ data access is pending (IEEG Portal / ieeg.org),
  solo and pre-mentor.

## 7. What carries to Phase 4

1. **ADR 0014** — the pre-registered peripheral-SOZ hypothesis (low degree →
   SOZ, Engel-I primary, AUC metric, one-sided Wilcoxon; direction locked).
2. **Held-out confirmatory dataset** — target: HUP ictal iEEG via ieeg.org;
   inclusion/windowing to be fixed in a follow-up ADR before any AUC is
   computed on it.
3. **Interictal generalisation** (HUP / ds004100) as a secondary, explicitly
   different question (the mechanism is partly ictal).
4. The frozen pipeline (loader → connectivity → Layer P → degree) and the
   SOZ-fraction moderator analysis.

## 8. One-paragraph verdict

Phase 3 pre-registered and then rejected the project's central quantum
hypothesis: on a multi-centre open iEEG dataset, the CTQW participation-ratio
ranking does not localize the clinician SOZ above chance (Engel I, n = 17,
median AUC 0.39, one-sided Wilcoxon p = 0.96), and is not superior to classical
centralities. In its wake, an exploratory **classical** signal emerged — SOZ
electrodes sit in a peripheral, low-degree position in the ictal network
(reverse-direction AUC ≈ 0.78) — which is robust to the electrode-type confound
(all-ECoG cohort), strengthens significantly during the seizure (ictal vs
pre-ictal p = 0.013) on top of a pre-ictal baseline tendency, and is heterogeneous
across patients. The quantum thesis is rejected; the project's live, pre-registered
question (ADR 0014) is now classical, awaiting a held-out ictal + SOZ dataset for
confirmation.