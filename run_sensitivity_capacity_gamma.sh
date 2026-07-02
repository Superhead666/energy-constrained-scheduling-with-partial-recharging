#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${SCRIPT_DIR}/bin/caie_part1_final"

if [[ ! -x "${BIN}" ]]; then
  "${SCRIPT_DIR}/build_cplex.sh" >/dev/null
fi

mkdir -p "${SCRIPT_DIR}/results" "${SCRIPT_DIR}/logs"

N_EXACT="${N_EXACT:-50}"
N_ROBUST="${N_ROBUST:-100}"
TAU="${TAU:-0.5}"
R="${R:-0.4}"
GAMMA_VALUES="${GAMMA_VALUES:-2,2.5,3,3.5,4,4.5,5,5.5,6}"
DUE_REFERENCE_GAMMA="${DUE_REFERENCE_GAMMA:-4}"
P_MIN="${P_MIN:-2}"
P_MAX="${P_MAX:-30}"
E_MIN="${E_MIN:-2}"
E_MAX="${E_MAX:-50}"
ALPHA="${ALPHA:-10}"
BETA="${BETA:-1}"
REPS="${REPS:-30}"
TIME_EXACT="${TIME_EXACT:-600}"
TIME_HEURISTIC="${TIME_HEURISTIC:-30}"

EXACT_OUT="${SCRIPT_DIR}/results/sensitivity_capacity_gamma_milp_sh_n50.csv"
ROBUST_OUT="${SCRIPT_DIR}/results/sensitivity_capacity_gamma_grasp_pr_n100.csv"
LOG_FILE="${SCRIPT_DIR}/logs/sensitivity_capacity_gamma.log"

{
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting capacity sensitivity"
  echo "Exact output: ${EXACT_OUT}"
  echo "Robustness output: ${ROBUST_OUT}"
  echo "Exact grid: n=${N_EXACT}; tau=${TAU}; R=${R}; gamma=${GAMMA_VALUES}; due-reference-gamma=${DUE_REFERENCE_GAMMA}; reps=${REPS}; time=${TIME_EXACT}s"
  echo "Robustness grid: n=${N_ROBUST}; tau=${TAU}; R=${R}; gamma=${GAMMA_VALUES}; due-reference-gamma=${DUE_REFERENCE_GAMMA}; reps=${REPS}; time=${TIME_HEURISTIC}s"
} | tee "${LOG_FILE}"

"${BIN}" \
  --mode part1 \
  --solver milp-heur \
  --due-mode taur \
  --n-values "${N_EXACT}" \
  --tau-values "${TAU}" \
  --R-values "${R}" \
  --gamma-values "${GAMMA_VALUES}" \
  --due-reference-gamma "${DUE_REFERENCE_GAMMA}" \
  --common-base-seeds \
  --p-min "${P_MIN}" \
  --p-max "${P_MAX}" \
  --e-min "${E_MIN}" \
  --e-max "${E_MAX}" \
  --alpha "${ALPHA}" \
  --beta "${BETA}" \
  --reps "${REPS}" \
  --time "${TIME_EXACT}" \
  --out "${EXACT_OUT}" 2>&1 | tee -a "${LOG_FILE}"

"${BIN}" \
  --mode part1 \
  --solver grasp-pr \
  --due-mode taur \
  --n-values "${N_ROBUST}" \
  --tau-values "${TAU}" \
  --R-values "${R}" \
  --gamma-values "${GAMMA_VALUES}" \
  --due-reference-gamma "${DUE_REFERENCE_GAMMA}" \
  --common-base-seeds \
  --p-min "${P_MIN}" \
  --p-max "${P_MAX}" \
  --e-min "${E_MIN}" \
  --e-max "${E_MAX}" \
  --alpha "${ALPHA}" \
  --beta "${BETA}" \
  --reps "${REPS}" \
  --time "${TIME_HEURISTIC}" \
  --out "${ROBUST_OUT}" 2>&1 | tee -a "${LOG_FILE}"

"${SCRIPT_DIR}/summarize_sensitivity_capacity.py" | tee -a "${LOG_FILE}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished capacity sensitivity" | tee -a "${LOG_FILE}"
