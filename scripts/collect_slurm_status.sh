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

printf 'job_id\tstate\texit_code\telapsed\talloc_cpus\tmax_rss_kib\tmax_rss\treq_mem\treq_tres\talloc_tres\tsacct_version\tsacct_command\n' \
    > "${temporary_path}"

SACCT_FORMAT='JobID,State,ExitCode,Elapsed,AllocCPUS,MaxRSS,ReqMem,ReqTRES,AllocTRES'
SACCT_COMMAND="sacct -X -j ${JOB_ID} --format=${SACCT_FORMAT} --noheader --parsable2"
SACCT_VERSION="$(sacct --version 2>&1 | tr '\n' ' ' | sed 's/[[:space:]]*$//')"

# Keep allocation state/exit fields from -X, but query all rows separately:
# MaxRSS is normally recorded on the srun .0 step rather than the allocation.
allocation="$(sacct -X -j "${JOB_ID}" --format="${SACCT_FORMAT}" --noheader --parsable2)"
accounting="$(sacct -j "${JOB_ID}" --format="${SACCT_FORMAT}" --noheader --parsable2)"
if [[ -z "${allocation//[[:space:]]/}" || -z "${accounting//[[:space:]]/}" ]]; then
    echo "error: no completed Slurm accounting record found for ${JOB_ID}" >&2
    exit 1
fi

allocation_row="$(awk -F'|' 'NF >= 9 {print; exit}' <<< "${allocation}")"
IFS='|' read -r _ state exit_code elapsed alloc_cpus _ req_mem req_tres alloc_tres <<< "${allocation_row}"

read -r max_rss_kib max_rss <<< "$(awk -F'|' '
    function to_kib(value, number, unit) {
        if (value == "" || value == "N/A" || value == "Unknown") return 0
        match(value, /^[0-9.]+/)
        number = substr(value, RSTART, RLENGTH) + 0
        unit = substr(value, RLENGTH + 1, 1)
        if (unit == "M") return number * 1024
        if (unit == "G") return number * 1024 * 1024
        if (unit == "T") return number * 1024 * 1024 * 1024
        if (unit == "K") return number
        return number / 1024
    }
    NF >= 9 {
        value = to_kib($6)
        if (value > maximum) {
            maximum = value
            raw = $6
        }
    }
    END { print int(maximum + 0.5), raw }
' <<< "${accounting}")"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${JOB_ID}" "${state}" "${exit_code}" "${elapsed}" "${alloc_cpus}" \
    "${max_rss_kib}" "${max_rss}" "${req_mem}" "${req_tres}" "${alloc_tres}" \
    "${SACCT_VERSION}" "${SACCT_COMMAND}" >> "${temporary_path}"

if [[ "$(wc -l < "${temporary_path}")" -lt 2 ]]; then
    echo "error: Slurm returned no usable accounting record for ${JOB_ID}" >&2
    exit 1
fi

mv -- "${temporary_path}" "${OUTPUT_PATH}"
trap - EXIT
echo "wrote ${OUTPUT_PATH}"
