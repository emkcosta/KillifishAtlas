#!/bin/bash
#SBATCH -J bayes_brain                   # job name
#SBATCH -o logs/bayes_brain_%A.out       # stdout log
#SBATCH -e logs/bayes_brain_%A.err       # stderr log
#SBATCH -p batch                    # partition (use what your lab uses)
#SBATCH --account=twc               # account
#SBATCH -t 24:00:00                 # walltime (hh:mm:ss)
#SBATCH -N 1                        # nodes
#SBATCH --cpus-per-task=16          # CPUs for sklearn n_jobs=-1
#SBATCH --mem=64G                   # RAM

set -euo pipefail

# --- Paths (edit if needed) ---
PROJECT_DIR="/labs/twc/Emma/atlas_gdrive_backup/revisions/CodeCheck_revisions"   # project root
PY_SCRIPT="${PROJECT_DIR}/code/BayesAge/bayes_xu_brain.py"  # your script path

# Create logs dir if missing
mkdir -p "${PROJECT_DIR}/Xu_model_outputs/BayesAge/logs"

# --- Reproducibility + sane threading ---
export PYTHONHASHSEED=42
# Route parallelism to sklearn/joblib (n_jobs), not BLAS:
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# If you prefer to give each worker a couple BLAS threads, set the above to 2
# and reduce --cpus-per-task proportionally.

# --- Conda env (same as you've used before) ---
# Adjust these if your Anaconda path/env name differs
source /scg/apps/software/anaconda/2023.03-0/etc/profile.d/conda.sh
conda activate bayes

# --- Diagnostics (helpful in logs) ---
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "CWD before cd: $PWD"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-not_set}"
python --version

# --- Run ---
cd "${PROJECT_DIR}/Xu_model_outputs/BayesAge"
echo "CWD after cd: $PWD"

# srun binds the task to the allocated resources
srun python "${PY_SCRIPT}"

# --- Done ---
echo "Finished at: $(date)"

conda deactivate
