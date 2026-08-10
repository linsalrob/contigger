#!/usr/bin/env bash
# Collect final Slurm accounting fields for a completed Contigger job.
set -euo pipefail

usage() {
    echo "Usage: $0 JOB_ID OUTPUT_PATH" >&2
    echo "Collect final Slurm accounting fields for a completed Contigger job." >&2
}

if [[ $# -ne 2 ]]; then
    usage
    exit 2
fi

JOB_ID="$1"
OUTPUT_PATH="$2"

if ! [[ "${JOB_ID}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "error: JOB_ID must be a numeric Slurm job ID" >&2
    exit 2
fi
if ! command -v sacct >/dev/null 2>&1; then
    echo "error: sacct is required to collect completed Slurm accounting data" >&2
    exit 127
fi

mkdir -p "$(dirname -- "${OUTPUT_PATH}")"
temporary_path="$(mktemp "${OUTPUT_PATH}.tmp.XXXXXX")"
cleanup() {
    rm -f -- "${temporary_path}"
}
trap cleanup EXIT

printf 'job_id\tstate\texit_code\telapsed\talloc_cpus\tmax_rss\treq_mem\treq_tres\talloc_tres\n' \
    > "${temporary_path}"

# -X reports the job allocation rather than individual batch/step records. The
# parsable output is converted to TSV while preserving Slurm's resource strings.
accounting="$(sacct -X -j "${JOB_ID}" \
    --format=JobID,State,ExitCode,Elapsed,AllocCPUS,MaxRSS,ReqMem,ReqTRES,AllocTRES \
    --noheader --parsable2)"
if [[ -z "${accounting//[[:space:]]/}" ]]; then
    echo "error: no completed Slurm accounting record found for ${JOB_ID}" >&2
    exit 1
fi

awk -F'|' -v OFS='\t' -v job="${JOB_ID}" '
    NF >= 9 {
        sub(/[[:space:]]+$/, "", $1)
        print job, $2, $3, $4, $5, $6, $7, $8, $9
    }
' <<< "${accounting}" >> "${temporary_path}"

if [[ "$(wc -l < "${temporary_path}")" -lt 2 ]]; then
    echo "error: Slurm returned no usable accounting record for ${JOB_ID}" >&2
    exit 1
fi

mv -- "${temporary_path}" "${OUTPUT_PATH}"
trap - EXIT
echo "wrote ${OUTPUT_PATH}"
