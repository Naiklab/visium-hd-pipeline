#!/bin/bash
# =============================================================================
# LSF submit template — Visium HD pipeline (Minerva / Rocky9 / LSF)
# =============================================================================
# Copy this file, replace every <PLACEHOLDER>, then submit from your project
# directory (LSF writes the -oo/-eo logs relative to that directory):
#
#     mkdir -p logs
#     bsub < submit_<step>.sh
#
# Working resource profiles from real runs are in the README FAQ
# ("How much memory / walltime do I actually need?").
# =============================================================================

#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n <NCORES>                    # cores, e.g. 2
#BSUB -W <HH:MM>                     # walltime, e.g. 144:00 (premium max)
#BSUB -R "rusage[mem=<MEM_MB>]"      # MB PER CORE, e.g. 300000 = 300 GB/core
#BSUB -R "span[hosts=1]"             # keep all cores on one host
#BSUB -J "<job_name>"
#BSUB -oo logs/<job_name>.out
#BSUB -eo logs/<job_name>.err

# NOTE: we intentionally do NOT use `set -e` here — `module` is a shell
# function and a nonzero return from e.g. `module list` would abort the job
# before it starts. We capture and report the R exit code explicitly instead.

echo "========================================================================"
echo "Job $LSB_JOBID ($LSB_JOBNAME) | host $(hostname) | start $(date)"
echo "Working directory: $(pwd)"
echo "========================================================================"

# ---- Environment ------------------------------------------------------------
# ALWAYS module purge + load the Rocky9 R/4.4.1 module. The package library at
# /sc/arion/projects/naiklab/ikjot/R is ABI-compiled against THIS module — any
# other R (conda, CentOS7 builds) fails to load packages. Never use conda R.
module purge
module load R/4.4.1
module list

# ---- CRITICAL: cap BLAS / OpenMP threads ------------------------------------
# Any step that requests a single core (e.g. RCTD create.RCTD(max_cores = 1))
# will SEGFAULT if OpenBLAS spawns one thread per physical core — it exhausts
# RLIMIT_NPROC. Pin every BLAS/OMP backend to a single thread. This is the
# single most common cause of "the job died with no R error" on Minerva.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1

# ---- Run --------------------------------------------------------------------
# Heavy compute runs as a plain .R on the cluster. Render an .Rmd only for the
# lighter reporting steps. Edit these two lines for your step:
SCRIPT="<script.R>"        # e.g. Master_RCTD_Deconvolution_8um.R
CONFIG=""                  # optional, e.g. "config_visium_hd.R"

echo "Starting: $SCRIPT ..."
Rscript "$SCRIPT" $CONFIG
EXIT_CODE=$?

echo "========================================================================"
echo "Exit code: $EXIT_CODE | end $(date)"
echo "========================================================================"
exit $EXIT_CODE
