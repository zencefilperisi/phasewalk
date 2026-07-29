# Phasewalk

Continuous-time classical and quantum walks applied to intracranial
EEG recordings of epileptic seizures, with a stochastic neural-mass
model as the underlying theoretical scaffold.

The project frames a seizure as a bifurcation-induced transition in
excitatory-inhibitory dynamics, uses stochastic simulation to identify
early-warning signals of that transition, and applies quantum walks
on functional-connectivity graphs to model how ictal activity spreads
across brain regions.

## Scope and honest positioning

Two phases are complete and are pure replications of prior work:

- **Phase 1** reproduces the classical Wilson-Cowan analysis (Borisyuk
  & Kirillov 1992) — fixed points, subcritical Hopf bifurcation,
  bistability window.
- **Phase 2A** reproduces the stochastic Wilson-Cowan analysis of
  Negahbani et al. 2015 — noise-induced variance and lag-1
  autocorrelation growth as early-warning indicators of the
  bifurcation, and noise-induced escape from the low-activity basin.

Neither phase contributes new science; both build the tested
infrastructure needed for the phases below and serve as replication
controls. This is documented explicitly in
`theory/02_literature_review.md`.

Two remaining phases are where the project's potential contribution
lies:

- **Phase 2B** — continuous-time quantum walks on synthetic and real
  brain graphs. Core walk machinery (module + tests) is complete.
  Synthetic-graph explorations confirmed several known properties
  (probability conservation, oscillatory return, non-convergence to
  uniform) and one negative result: no simple node-level centrality
  measure reliably predicts quantum-walk spread across random graphs
  (see `decisions/0002-quantum-walk-conventions.md`).
- **Phase 2C** — application to real intracranial EEG (SWEC-ETHZ
  short-term dataset). In progress: data loaded and inspected;
  functional-connectivity extraction is the next step.

Phases 3 (variational quantum Wilson-Cowan) and 4 (QAOA feature
selection) from earlier plans have been dropped to concentrate on
depth in Phase 2. See `decisions/` for the reasoning.

## Roadmap

| Phase | Topic                                                            | Status              |
|-------|------------------------------------------------------------------|---------------------|
| 1     | Classical Wilson-Cowan: phase plane, bifurcation, fixed points   | complete            |
| 2A    | Stochastic Wilson-Cowan: critical slowing down, escape           | complete            |
| 2B    | Quantum walks on graphs: core machinery + synthetic exploration  | complete            |
| 2C    | Quantum walks on real iEEG functional-connectivity graphs        | in progress         |
| 3, 4  | (dropped in favour of deepening Phase 2 — see decision log)      | not planned         |

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

# Tests (all phases)
pytest tests/
```

Phase 2B/2C analysis scripts are still being developed; the tested
underlying functions live in `quantum_walks/`.

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
│                               # ratio, return probability, time-average
├── notebooks/                  # exploratory scripts (Jupyter-compatible)
│   ├── _style.py               # shared visual style (deep + cividis)
│   ├── simulate_oscillation.py
│   ├── visualize_phase_space.py
│   ├── analyze_bifurcation.py
│   ├── explore_stochastic_dynamics.py
│   ├── detect_critical_slowing.py
│   └── noise_induced_transitions.py
├── exploration/                # step-by-step learning scripts (micro-steps),
│                               # kept for transparency of the development
│                               # process; NOT part of the tested code path
├── decisions/                  # Architecture Decision Records (ADRs)
│                               # documenting every methodological choice
├── theory/                     # literature review, derivations, notes
├── tests/                      # pytest suite; 39 tests as of Phase 2B
├── requirements.txt
└── README.md
```

Not committed to the repository (excluded via `.gitignore`):

```
├── data/                       # raw iEEG data (SWEC-ETHZ .mat files)
└── phasewalk-env/              # local Python virtual environment
```

## Data

Phase 2C uses the **SWEC-ETHZ short-term iEEG dataset**
(<http://ieeg-swez.ethz.ch/>), 16 patients with 100 seizures,
sampled at 512 Hz, provided under a research-and-education licence.
Data is not redistributed with this repository; users must download
it directly from the source.

See `decisions/0003-switch-to-short-term-dataset.md` for the rationale
behind choosing the short-term over the long-term subset, and
`decisions/0004`, `0005`, `0006` for the methodological policies
governing how the data is handled.

## Decision log

Every non-trivial methodological choice is recorded in `decisions/`
as a standalone ADR (Architecture Decision Record). Current entries:

| ID   | Title                                        | Status              |
|------|----------------------------------------------|---------------------|
| 0001 | Dataset choice for Phase 2                   | superseded by 0003  |
| 0002 | Quantum-walk conventions                     | accepted            |
| 0003 | Switch to short-term dataset                 | accepted            |
| 0004 | Seizure duration policy                      | accepted            |
| 0005 | Patient-level analysis                       | accepted            |
| 0006 | ID1 as development patient                   | accepted            |

## Key references

- Wilson, H. R., & Cowan, J. D. (1972). Excitatory and inhibitory
  interactions in localized populations of model neurons.
  *Biophysical Journal*, 12(1), 1-24.
- Borisyuk, R. M., & Kirillov, A. B. (1992). Bifurcation analysis of
  a neural network model. *Biological Cybernetics*, 66(4), 319-325.
- Jirsa, V. K., Stacey, W. C., Quilichini, P. P., Ivanov, A. I., &
  Bernard, C. (2014). On the nature of seizure dynamics.
  *Brain*, 137(8), 2210-2230.
- Negahbani, E., Steyn-Ross, D. A., Steyn-Ross, M. L., Wilson, M. T.,
  & Sleigh, J. W. (2015). Noise-induced precursors of state
  transitions in the stochastic Wilson-Cowan model.
  *The Journal of Mathematical Neuroscience*, 5:9.
- Maturana, M. I., et al. (2020). Critical slowing down as a biomarker
  for seizure susceptibility. *Nature Communications*, 11:2172.
- Mulken, O., & Blumen, A. (2011). Continuous-time quantum walks:
  Models for coherent transport on complex networks.
  *Physics Reports*, 502(2-3), 37-87.
- Burrello, A., Cavigelli, L., Schindler, K., Benini, L., & Rahimi, A.
  (2019). Laelaps: An energy-efficient seizure detection algorithm
  from long-term human iEEG recordings without false alarms. In
  *DATE 2019* (best-paper nomination).
- Burrello, A., Schindler, K., Benini, L., & Rahimi, A. (2018).
  One-shot learning for iEEG seizure detection using end-to-end
  binary operations. In *BioCAS 2018* (best-paper award).

Full literature positioning is in `theory/02_literature_review.md`.

## Development notes

- Style: seaborn `deep` categorical palette + `cividis` sequential
  colormap; DPI 300 for saved figures. See `notebooks/_style.py`.
- Testing: `pytest tests/` — all changes must keep the test suite green.
- Reproducibility: deterministic components use fixed seeds
  (typically 42); stochastic components accept a `seed` argument and
  ensembles use independent child seeds derived from a master seed.

## License

MIT (to be added once the repository is public-ready).