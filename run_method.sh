#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${SCRIPT_DIR}/bin/caie_part1_final"

if [[ $# -lt 1 ]]; then
  echo "Usage: ./run_method.sh milp-s|milp-sh|grasp-pr" >&2
  exit 2
fi

METHOD="$1"

if [[ ! -x "${BIN}" ]]; then
  "${SCRIPT_DIR}/build_cplex.sh" >/dev/null
fi

mkdir -p "${SCRIPT_DIR}/results" "${SCRIPT_DIR}/logs"

case "${METHOD}" in
  milp-s|milps|strengthened)
    SOLVER="milp-s"
    OUT_FILE="${SCRIPT_DIR}/results/final_milp_s.csv"
    LABEL="MILP-S"
    ;;
  milp-sh|milpsh|milp-h|milph)
    SOLVER="milp-heur"
    OUT_FILE="${SCRIPT_DIR}/results/final_milp_sh_light_start.csv"
    LABEL="MILP-SH-light"
    ;;
  grasp-pr|grasp_pr|heuristic)
    SOLVER="grasp-pr"
    OUT_FILE="${SCRIPT_DIR}/results/final_grasp_pr.csv"
    LABEL="GRASP-PR"
    ;;
  *)
    echo "Unknown method: ${METHOD}" >&2
    echo "Usage: ./run_method.sh milp-s|milp-sh|grasp-pr" >&2
    exit 2
    ;;
esac

N_VALUES="${N_VALUES:-50,100,200}"
TAU_VALUES="${TAU_VALUES:-0.5,0.6}"
R_VALUES="${R_VALUES:-0.8}"
GAMMA_VALUES="${GAMMA_VALUES:-3,5}"
P_MIN="${P_MIN:-2}"
P_MAX="${P_MAX:-30}"
E_MIN="${E_MIN:-2}"
E_MAX="${E_MAX:-50}"
ALPHA="${ALPHA:-10}"
BETA="${BETA:-1}"
REPS="${REPS:-10}"

case "${SOLVER}" in
  milp-s|milp-heur)
    TIME_LIMIT="${TIME_LIMIT:-600}"
    ;;
  grasp-pr)
    TIME_LIMIT="${TIME_LIMIT:-30}"
    ;;
esac

LOG_FILE="${SCRIPT_DIR}/logs/${SOLVER}.log"

{
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting ${LABEL}"
  echo "Output: ${OUT_FILE}"
  echo "Grid: n=${N_VALUES}; tau=${TAU_VALUES}; R=${R_VALUES}; gamma=${GAMMA_VALUES}; reps=${REPS}; p=[${P_MIN},${P_MAX}]; e=[${E_MIN},${E_MAX}]; alpha=${ALPHA}; beta=${BETA}; time=${TIME_LIMIT}s"
} | tee "${LOG_FILE}"

"${BIN}" \
  --mode part1 \
  --solver "${SOLVER}" \
  --due-mode taur \
  --n-values "${N_VALUES}" \
  --tau-values "${TAU_VALUES}" \
  --R-values "${R_VALUES}" \
  --gamma-values "${GAMMA_VALUES}" \
  --p-min "${P_MIN}" \
  --p-max "${P_MAX}" \
  --e-min "${E_MIN}" \
  --e-max "${E_MAX}" \
  --alpha "${ALPHA}" \
  --beta "${BETA}" \
  --reps "${REPS}" \
  --time "${TIME_LIMIT}" \
  --out "${OUT_FILE}" 2>&1 | tee -a "${LOG_FILE}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished ${LABEL}" | tee -a "${LOG_FILE}"
