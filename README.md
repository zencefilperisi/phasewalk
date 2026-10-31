# Phasewalk

Continuous-time classical and quantum walks applied to intracranial
EEG recordings of epileptic seizures, with a stochastic neural-mass
model as the underlying theoretical scaffold.

The project frames a seizure as a bifurcation-induced transition in
excitatory-inhibitory dynamics, uses stochastic simulation to identify
early-warning signals of that transition, and applies quantum walks
on functional-connectivity graphs to ask whether they carry information
about where and how ictal activity is organised across brain regions.

## Scope and honest positioning

Phases 1 and 2A are pure replications of prior work and contribute no
new science; they build the tested infrastructure used later and serve
as replication controls:

- **Phase 1** reproduces the classical Wilson-Cowan analysis (Borisyuk
  & Kirillov 1992) — fixed points, subcritical Hopf bifurcation,
  bistability window.
- **Phase 2A** reproduces the stochastic Wilson-Cowan analysis of
  Negahbani et al. 2015 — noise-induced variance and lag-1
  autocorrelation growth as early-warning indicators, and noise-induced
  escape from the low-activity basin.

The quantum-walk line of inquiry (Phases 2B, 2C, 3) is where the
project's potential contribution lay. Its results, stated plainly, are
**negative or null**:

- **Phase 2B — quantum walks on graphs (complete).** The core walk
  machinery (module + tests) reproduces known properties (probability
  conservation, oscillatory return, non-convergence to uniform). One
  negative result on synthetic graphs: no simple node-level centrality
  reliably predicts continuous-time quantum-walk (CTQW) spread
  (`decisions/0002`).
- **Phase 2C — CTQW on real iEEG, SWEC-ETHZ (complete).** On the
  short-term dataset, the CTQW per-electrode participation-ratio (PR)
  ranking is reproducible across a patient's seizures beyond chance
  (permutation test, p < 1e-4) and is distinct from classical
  centralities (near-orthogonal to eigenvector centrality) — but it is
  **not** more reproducible than them, and the dataset carries no
  electrode-level seizure-onset-zone (SOZ) marks, so no clinical claim
  is possible. See `decisions/0007`–`0010` and
  `phase2c_interim_summary.md`.
- **Phase 3 — CTQW clinical validation, ds003029 (complete).** A
  pre-registered test on a multi-centre open iEEG dataset of whether
  low-PR (localized) electrodes identify the clinician SOZ. The analysis
  plan — metric, population, statistic and decision rule — was fixed
  **before** any cohort result was seen (`decisions/0013`). **Result:
  negative.** In the primary population (seizure-free / Engel I,
  n = 17), CTQW-PR does not localize the SOZ above chance (median
  AUC 0.39; one-sided Wilcoxon p = 0.96; label-permutation null below
  chance), and it is not superior to classical centralities. An
  exploratory, *not* pre-registered observation: low-degree / low-
  eigenvector (peripheral) electrodes track the SOZ in the reverse
  direction (AUC ≈ 0.78) — a hypothesis for future pre-registered work,
  not a claim here.

**Bottom line.** Across synthetic graphs and two independent real iEEG
datasets, continuous-time quantum walks have not provided a clinically
useful, or classically superior, signal for seizure-onset localization.
The value of this repository is a rigorously documented, pre-registered,
honestly-reported negative line of inquiry, with a reproducible and
tested pipeline that the next (likely classical) question can build on.

## Roadmap

| Phase | Topic                                                                 | Status             |
|-------|-----------------------------------------------------------------------|--------------------|
| 1     | Classical Wilson-Cowan: phase plane, bifurcation, fixed points        | complete           |
| 2A    | Stochastic Wilson-Cowan: critical slowing down, escape                | complete           |
| 2B    | Quantum walks on graphs: core machinery + synthetic exploration       | complete           |
| 2C    | CTQW on real iEEG (SWEC-ETHZ): reproducibility + characterization      | complete           |
| 3     | CTQW clinical validation (ds003029, multi-centre, pre-registered)     | complete — negative|
| 4     | (earlier QAOA / variational plans dropped — see decision log)         | not planned        |

## Installation

```bash
git clone https://github.com/<user>/phasewalk.git
cd phasewalk

python -m venv phasewalk-env
source phasewalk-env/bin/activate    # Windows: phasewalk-env\Scripts\activate

pip install -r requirements.txt
```

## Usage

```bash
# Phase 1 — deterministic Wilson-Cowan
python notebooks/simulate_oscillation.py
python notebooks/visualize_phase_space.py
python notebooks/analyze_bifurcation.py

# Phase 2A — stochastic Wilson-Cowan
python notebooks/explore_stochastic_dynamics.py
python notebooks/detect_critical_slowing.py
python notebooks/noise_induced_transitions.py

# Phase 2C — CTQW on SWEC-ETHZ iEEG (step-by-step in exploration/)
#   see the micro-step scripts and phase2c_interim_summary.md

# Phase 3 — CTQW clinical validation on ds003029
python exploration/ds003029/ds003029_cohort_build.py      # freeze the cohort
python exploration/ds003029/ds003029_cohort_analysis.py   # pre-registered analysis
python exploration/ds003029/ds003029_plot_engelI_auc.py   # primary-result figure

# Tests (all phases)
pytest tests/
```

The tested underlying functions live in `classical/`, `quantum_walks/`
and `data/`; scripts under `exploration/` are the development record and
are not part of the tested code path.

## Layout

```
phasewalk/
├── classical/                  # Phase 1 and 2A: (stochastic) Wilson-Cowan
│   ├── model.py                # equations and deterministic integrator
│   ├── analysis.py             # fixed points, eigenvalues, Hopf threshold
│   ├── stochastic.py           # Euler-Maruyama SDE integrator, ensemble,
│   │                           # critical-slowing-down metrics
│   └── visualize.py            # reusable plotting helpers
├── quantum_walks/              # Phase 2B: continuous-time walks on graphs
│   └── walks.py                # classical diffusion, CTQW, participation
│                               # ratio, return probability, time-average,
│                               # and a batched eigendecomposition routine
├── data/                       # Phase 2C / 3: iEEG connectivity + loaders
│   ├── connectivity.py         # functional-connectivity matrices (Pearson,
│   │                           # PLV) and FCM -> graph (Layer P / Layer S)
│   └── loader_bids.py          # ds003029 BIDS loader: channel + SOZ
│                               # inclusion (ADR 0011) and ictal window
│                               # (ADR 0012)
├── notebooks/                  # Phase 1/2A exploratory scripts
├── exploration/                # step-by-step development scripts, kept for
│   └── ds003029/               # transparency; NOT the tested code path
│                               #   Phase 3 cohort build + analysis + figure,
│                               #   frozen cohort and per-subject results
├── decisions/                  # Architecture Decision Records (ADRs)
├── theory/                     # literature review, derivations, notes
├── tests/                      # pytest suite (see note below)
├── requirements.txt
└── README.md
```

Not committed to the repository (excluded via `.gitignore`):

```
├── data/raw, ds003029/, *.eeg, *.mat   # raw iEEG recordings (download from source)
└── phasewalk-env/                      # local Python virtual environment
```

## Data

Two iEEG datasets are used; neither is redistributed here — download each
from its source under its own licence.

- **Phase 2C — SWEC-ETHZ short-term iEEG dataset**
  (<http://ieeg-swez.ethz.ch/>): 16 patients, 100 seizures, 512 Hz,
  research-and-education licence. Used for methodology and
  reproducibility; it carries no electrode-level SOZ annotations.
- **Phase 3 — Epilepsy-iEEG-Multicenter dataset (ds003029)** on
  OpenNeuro (Li et al. 2021, "neural fragility"): multi-centre ECoG/SEEG
  in BrainVision format, with clinician SOZ annotations and surgical
  outcome (Engel). Public; the raw signal is not stored in this repo.
  Inclusion, windowing and the pre-registered analysis plan are in
  `decisions/0011`–`0013`.

## Decision log

Every non-trivial methodological choice is a standalone ADR in
`decisions/`.

| ID   | Title                                                      | Status             |
|------|------------------------------------------------------------|--------------------|
| 0001 | Dataset choice for Phase 2                                 | superseded by 0003 |
| 0002 | Quantum-walk conventions                                   | accepted           |
| 0003 | Switch to short-term dataset                               | accepted           |
| 0004 | Seizure duration policy                                    | accepted           |
| 0005 | Patient-level analysis                                     | accepted           |
| 0006 | ID1 as development patient (parameter independence)        | accepted           |
| 0007 | Connectivity measure: Pearson primary, PLV robustness      | accepted           |
| 0008 | Graph construction: Layer P (positive) vs Layer S (signed) | accepted           |
| 0009 | Reframe: CTQW participation-ratio ranking as primary object| accepted           |
| 0010 | SOZ unavailable on SWEC-ETHZ: staged clinical strategy     | accepted           |
| 0011 | ds003029 channel and SOZ inclusion policy                  | accepted           |
| 0012 | ds003029 ictal window and seizure-onset definition         | accepted           |
| 0013 | ds003029 analysis plan (pre-registration)                  | accepted           |

## Key references

- Wilson, H. R., & Cowan, J. D. (1972). Excitatory and inhibitory
  interactions in localized populations of model neurons.
  *Biophysical Journal*, 12(1), 1-24.
- Borisyuk, R. M., & Kirillov, A. B. (1992). Bifurcation analysis of a
  neural network model. *Biological Cybernetics*, 66(4), 319-325.
- Jirsa, V. K., et al. (2014). On the nature of seizure dynamics.
  *Brain*, 137(8), 2210-2230.
- Negahbani, E., et al. (2015). Noise-induced precursors of state
  transitions in the stochastic Wilson-Cowan model. *The Journal of
  Mathematical Neuroscience*, 5:9.
- Mulken, O., & Blumen, A. (2011). Continuous-time quantum walks: Models
  for coherent transport on complex networks. *Physics Reports*,
  502(2-3), 37-87.
- Li, A., et al. (2021). Neural fragility as an EEG marker of the seizure
  onset zone. *Nature Neuroscience*, 24, 1465-1474. (ds003029)
- Burrello, A., et al. (2019). Laelaps: an energy-efficient seizure
  detection algorithm from long-term human iEEG recordings without false
  alarms. *DATE 2019*.

Full literature positioning is in `theory/02_literature_review.md`.

## Development notes

- Testing: `pytest tests/` — all changes keep the suite green (96 tests as of Phase 3).
- Reproducibility: deterministic components use fixed seeds; stochastic
  components accept a `seed`; the Phase 3 cohort is frozen
  (`exploration/ds003029/cohort_frozen.json`) before any result is
  computed, so the analysed cohort cannot drift after the fact.
- Style: seaborn `deep` palette + `cividis` colormap; see
  `notebooks/_style.py`.

## License

MIT.