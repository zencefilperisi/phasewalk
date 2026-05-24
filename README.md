# Phasewalk

A study of epileptic seizure dynamics through the lens of quantum
computation. The project frames a seizure as a **phase transition** in a
neural-mass model and explores three quantum tools for analysing it:

1. **Quantum simulation of the Wilson-Cowan model** via variational
   quantum circuits.
2. **Quantum walks** as a model of seizure propagation across a network.
3. **QAOA-based feature selection** for seizure prediction from EEG.

## Roadmap

| Phase | Topic                                                        | Status      |
|-------|--------------------------------------------------------------|-------------|
| 1     | Classical Wilson-Cowan: simulation, phase plane, bifurcation | complete    |
| 2     | Quantum walks — propagation on a small network               | pending     |
| 3     | Variational quantum formulation of Wilson-Cowan              | pending     |
| 4     | QAOA feature selection for seizure prediction (CHB-MIT/Bonn) | pending     |

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
# Time-series and trajectory in the oscillatory regime
python notebooks/simulate_oscillation.py

# Streamlines, nullclines, and trajectory overlays
python notebooks/visualize_phase_space.py

# Bifurcation diagram and eigenvalue tracking
python notebooks/analyze_bifurcation.py

# Tests
pytest tests/
```

## Layout

```
phasewalk/
├── classical/                  # Phase 1: classical Wilson-Cowan
│   ├── model.py                # equations and integrator
│   ├── analysis.py             # fixed-point search, eigenvalues, Hopf
│   └── visualize.py            # reusable plotting helpers
├── notebooks/                  # exploratory scripts (Jupyter-compatible)
│   ├── _style.py               # shared visual style
│   ├── simulate_oscillation.py
│   ├── visualize_phase_space.py
│   └── analyze_bifurcation.py
├── theory/                     # LaTeX/Markdown derivations and notes
├── tests/                      # pytest suite
├── requirements.txt
└── README.md
```

Modules added in later phases:

```
├── quantum_walks/     # Phase 2
├── quantum_wc/        # Phase 3
├── seizure_opt/       # Phase 4
└── data/              # EEG datasets (excluded by .gitignore)
```

## References

- Wilson, H. R., & Cowan, J. D. (1972). Excitatory and inhibitory
  interactions in localized populations of model neurons.
  *Biophysical Journal*, 12(1), 1-24.
- Jirsa, V. K., Stacey, W. C., Quilichini, P. P., Ivanov, A. I., &
  Bernard, C. (2014). On the nature of seizure dynamics.
  *Brain*, 137(8), 2210-2230.
- Strogatz, S. H. (2018). *Nonlinear Dynamics and Chaos* (2nd ed.).
  CRC Press.

## License

MIT (to be added once the repository is published).
