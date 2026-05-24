# Exploration scripts

These `micro_step_*.py` scripts are the step-by-step learning process
through which the quantum-walk machinery was first built and understood,
before being consolidated into the tested `quantum_walks` module.

They are kept for transparency and as a record of how the understanding
was developed. They are NOT part of the production code path and are not
covered by the test suite.

| Script | What it introduced |
|--------|--------------------|
| `micro_step_1.py` | Building and drawing a small-world graph (NetworkX) |
| `micro_step_2.py` | Extracting adjacency and Laplacian matrices |
| `micro_step_3.py` | Classical diffusion via the matrix exponential |
| `micro_step_4.py` | Quantum walk and side-by-side comparison |
| `micro_step_5.py` | Probability conservation, hub vs periphery, return probability |

The consolidated, tested versions of these dynamics live in
`../quantum_walks/walks.py`.
