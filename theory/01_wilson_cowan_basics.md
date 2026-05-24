# Phase 1 — Foundations of the Wilson-Cowan Model

## 1. The Model

Mean-field dynamics for two coupled populations (excitatory E,
inhibitory I):

```
tau_E * dE/dt = -E + S_E(w_EE * E - w_EI * I + P)
tau_I * dI/dt = -I + S_I(w_IE * E - w_II * I + Q)
```

Sigmoid activation:

```
S(x) = 1 / (1 + exp(-a * (x - theta)))
```

Both ``E(t)`` and ``I(t)`` lie in ``[0, 1]``, representing the average
firing rate of each population. Sigmoid saturation guarantees that the
state remains in this interval.

## 2. Parameters and Typical Ranges

| Symbol  | Meaning                          | Typical range |
|---------|----------------------------------|---------------|
| w_EE    | E → E coupling                   | 10 – 16       |
| w_EI    | I → E coupling                   | 12 – 15       |
| w_IE    | E → I coupling                   | 10 – 15       |
| w_II    | I → I coupling                   | 3 – 8         |
| P       | external drive to E              | variable      |
| Q       | external drive to I              | usually 0     |
| tau_E   | E time constant                  | 1 (norm.)     |
| tau_I   | I time constant                  | 2 (norm.)     |
| a, theta| sigmoid slope and threshold      | 1.3, 4.0 etc. |

## 3. Fixed-Point Analysis

Fixed points ``(E*, I*)`` satisfy ``dE/dt = 0`` and ``dI/dt = 0``
simultaneously. Because the equations are nonlinear, we solve them
numerically:

- **Nullclines.** The E-nullcline is the locus of points where
  ``dE/dt = 0``; similarly for the I-nullcline.
- **Intersections.** Points where both nullclines cross are fixed
  points.
- **Numerical contours.** ``matplotlib.contour`` traces the zero-level
  set directly — no analytic inversion is needed.

## 4. Linearisation and Stability

Near a fixed point ``(E*, I*)``, write ``E = E* + e`` and
``I = I* + i``. To leading order in the perturbation:

```
d/dt [e, i]^T = J [e, i]^T,    J = Jacobian evaluated at (E*, I*)
```

Eigenvalues of J classify the fixed point:

| Eigenvalues                                  | Type           | Flow              |
|----------------------------------------------|----------------|-------------------|
| Re(λ₁), Re(λ₂) < 0, complex conjugate pair   | Stable focus   | Inward spiral     |
| Re(λ₁), Re(λ₂) < 0, real                     | Stable node    | Direct inflow     |
| Re(λ₁), Re(λ₂) > 0                           | Unstable       | Outward flow      |
| λ₁ > 0, λ₂ < 0                               | Saddle         | Mixed             |
| Re(λ₁) = Re(λ₂) = 0                          | Centre / Hopf  | Marginal oscill.  |

## 5. Bifurcations

**Hopf bifurcation.** As a control parameter is varied, a complex pair
of eigenvalues crosses the imaginary axis. The fixed point loses
stability and a limit cycle is born. Clinical interpretation: **seizure
onset**.

- *Supercritical Hopf*: the limit cycle grows continuously from zero
  amplitude. Matches seizures preceded by an aura.
- *Subcritical Hopf*: bistability between a stable fixed point and a
  large-amplitude limit cycle; transitions are abrupt. Matches
  seizures with sudden onset.

**Saddle-node.** Two fixed points collide and annihilate.

## 6. Clinical Correspondence

- *Healthy resting state*: a single globally attracting fixed point.
- *Pre-ictal*: the system drifts toward a bifurcation; statistical
  signatures appear in the EEG (rising variance, slowing autocorrelation
  — known as critical slowing down).
- *Ictal*: the system is on a limit cycle; periodic bursts are visible
  as synchronous rhythmic waves on EEG.
- *Post-ictal*: the system returns to a fixed point but parameters may
  remain near the bifurcation threshold.

## 7. Phase 1 Goals

- [x] Implement the Wilson-Cowan equations (`classical/model.py`)
- [x] Phase portrait and nullcline visualisation (`classical/visualize.py`)
- [ ] Fixed-point finder with Jacobian/eigenvalue classification
      (`classical/analysis.py`)
- [ ] One-parameter bifurcation diagram (sweeping P)
- [ ] Two-parameter bifurcation map (in the (P, w_EI) plane)
- [ ] Stochastic version (additive Wiener noise) — bridge to Phase 2

## References

1. Wilson, H. R., & Cowan, J. D. (1972). Excitatory and inhibitory
   interactions in localized populations of model neurons.
   *Biophysical Journal*, 12(1), 1-24.
2. Borisyuk, R. M., & Kirillov, A. B. (1992). Bifurcation analysis of a
   neural network model. *Biological Cybernetics*, 66(4), 319-325.
3. Jirsa, V. K., Stacey, W. C., Quilichini, P. P., Ivanov, A. I., &
   Bernard, C. (2014). On the nature of seizure dynamics.
   *Brain*, 137(8), 2210-2230.
4. Strogatz, S. H. (2018). *Nonlinear Dynamics and Chaos*. CRC Press.
