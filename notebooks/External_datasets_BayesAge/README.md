# BayesAge 2.0 models for external datasets

This workflow can be run in two equivalent ways, and this is demonstrated using Xu et al data:

---

## Option A — Jupyter notebook
Use this if you prefer an interactive environment or want to visualize intermediate steps.

- Open the notebook in an interactive jupyter session.

## Option B — Script + Bash
Use this if you prefer command-line execution (e.g., on HPC or in a pipeline).

- In your chosen compute environment, run the following. We did this on a HPC. For example:

```sbatch root_dir/scripts/bash/run_bayes_brain.sh
```