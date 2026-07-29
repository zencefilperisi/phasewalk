# Exploration scripts

These `micro_step_*.py` scripts are the step-by-step learning process
through which the project's methodology was first developed and
understood, before being consolidated into the tested modules and
production code paths.

Scripts fall into two groups:

- Scripts 1-9 built the quantum-walk machinery on synthetic graphs.
  Their consolidated, tested output lives in `../quantum_walks/walks.py`.
- Scripts 10-11 are the first contact with real intracranial EEG data
  from the SWEC-ETHZ short-term dataset. Their output will be
  consolidated into a `../data/` module during Phase 2C.

They are kept in the repository for transparency and as a record of
how the understanding was developed. They are NOT part of the tested
production code path.

| Script | What it introduced |
|--------|--------------------|
| `micro_step_1.py`  | Building and drawing a small-world graph (NetworkX) |
| `micro_step_2.py`  | Extracting adjacency and Laplacian matrices |
| `micro_step_3.py`  | Classical diffusion via the matrix exponential |
| `micro_step_4.py`  | Quantum walk and side-by-side comparison with classical |
| `micro_step_5.py`  | Probability conservation, hub vs periphery, return probability |
| `micro_step_6.py`  | Time-averaged distribution; resolving the hub vs periphery question |
| `micro_step_7.py`  | Four node centralities (degree, betweenness, eigenvector, clustering) as predictors of spread |
| `micro_step_8.py`  | Robustness check across 100 graphs — clustering effect is not consistent |
| `micro_step_9.py`  | Spectral participation as a much stronger (but tautological — see decision log) predictor |
| `micro_step_10.py` | First inspection of a real iEEG seizure recording |
| `micro_step_11.py` | Bulk scan and shape check across all seizure files for ID1, ID2, ID3 |