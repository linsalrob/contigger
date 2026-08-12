#!/bin/bash --login
# Bounded candidate-only Contigger preflight for a multi-assembly manifest.
#
# Example:
#   sbatch scripts/profile_manifest_candidates_setonix.slurm MANIFEST OUTPUT_JSON
#
# This does not invoke minimap2 or produce biological merge outputs. Both limits are
# required so an unfamiliar collection stops before unbounded seed/pair expansion.
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=contigger-preflight

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 MANIFEST OUTPUT_JSON" >&2
    exit 2
fi

MANIFEST="$1"
OUTPUT="$2"
CONTIGGER_ROOT="${CONTIGGER_ROOT:-/scratch/pawsey1018/edwa0468/GitHubs/contigger}"
PYTHON_BIN="${CONTIGGER_PYTHON:-/scratch/pawsey1018/edwa0468/software/miniconda3/atavide_lite/bin/python}"
SEED_LIMIT="${CONTIGGER_MAX_SEED_PAIR_OBSERVATIONS:-100000000}"
CANDIDATE_LIMIT="${CONTIGGER_MAX_CANDIDATE_PAIRS:-1000000}"

mkdir -p "$(dirname -- "${OUTPUT}")"
echo "manifest=${MANIFEST}"
echo "output=${OUTPUT}"
echo "candidate_shards=64 max_seed_pair_observations=${SEED_LIMIT} max_candidate_pairs=${CANDIDATE_LIMIT}"

srun --nodes=1 --ntasks=1 --cpus-per-task=1 -m block:block:block env \
    PYTHONPATH="${CONTIGGER_ROOT}/src" \
    "${PYTHON_BIN}" "${CONTIGGER_ROOT}/benchmarking/profile_manifest_candidates.py" \
    --manifest "${MANIFEST}" --output "${OUTPUT}" --candidate-shards 64 \
    --max-seed-pair-observations "${SEED_LIMIT}" --max-candidate-pairs "${CANDIDATE_LIMIT}"
