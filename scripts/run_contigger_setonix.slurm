#!/bin/bash --login
#
# Example Pawsey Setonix CPU job for Contigger.
# Replace `your-project` with your Pawsey project code before submitting.
#
#   sbatch scripts/run_contigger_setonix.slurm
#   sbatch scripts/run_contigger_setonix.slurm samples.tsv results/contigger
#
#SBATCH --account=your-project
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=58880M
#SBATCH --time=01:00:00
#SBATCH --job-name=contigger
#SBATCH --output=logs/contigger-%j.out
#SBATCH --error=logs/contigger-%j.err

set -euo pipefail

MANIFEST="${1:-samples.tsv}"
OUTPUT_PREFIX="${2:-results/contigger}"
THREADS="${SLURM_CPUS_PER_TASK:-32}"
EVIDENCE="${CONTIGGER_EVIDENCE:-none}"
INDEX_DIR="${CONTIGGER_INDEX_DIR:-${SCRATCH:-/tmp}/contigger-indexes/${SLURM_JOB_ID:-local}}"

eval "$(conda shell.bash hook)"
conda activate contigger

command -v contigger >/dev/null || {
    echo "contigger is not on PATH; install it in the selected environment" >&2
    exit 127
}
command -v minimap2 >/dev/null || {
    echo "minimap2 is not on PATH; install it in the selected environment" >&2
    exit 127
}
if [[ "${EVIDENCE}" == "alignments" ]]; then
    command -v samtools >/dev/null || {
        echo "samtools is required for --evidence alignments" >&2
        exit 127
    }
fi

mkdir -p "$(dirname "${OUTPUT_PREFIX}")" "${INDEX_DIR}" logs
echo "Job ${SLURM_JOB_ID:-local} on ${HOSTNAME}"
echo "Manifest: ${MANIFEST}; output: ${OUTPUT_PREFIX}"
echo "Threads: ${THREADS}; evidence: ${EVIDENCE}"

srun --nodes=1 --ntasks=1 --cpus-per-task="${THREADS}" -m block:block:block \
    contigger merge \
    --manifest "${MANIFEST}" \
    --output-prefix "${OUTPUT_PREFIX}" \
    --threads "${THREADS}" \
    --index-dir "${INDEX_DIR}" \
    --evidence "${EVIDENCE}"
