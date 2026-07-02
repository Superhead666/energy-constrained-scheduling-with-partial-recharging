#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${SCRIPT_DIR}/bin/caie_part1_final"

if [[ ! -x "${BIN}" ]]; then
  "${SCRIPT_DIR}/build_cplex.sh" >/dev/null
fi

N="${N:-40}"
TAU="${TAU:-0.5}"
R="${R:-0.4}"
GAMMA_VALUES="${GAMMA_VALUES:-3,5}"
RHO_VALUES="${RHO_VALUES:-0.2,0.6}"
REPS="${REPS:-20}"
TIME_LIMIT="${TIME_LIMIT:-180}"
WORKERS="${WORKERS:-4}"
DUE_REFERENCE_GAMMA="${DUE_REFERENCE_GAMMA:-4}"

echo "Running partial-vs-full recharge sensitivity"
echo "Grid: n=${N}; tau=${TAU}; R=${R}; gamma=${GAMMA_VALUES}; rho=${RHO_VALUES}; reps=${REPS}; time=${TIME_LIMIT}s"

N="${N}" \
TAU="${TAU}" \
R="${R}" \
GAMMA_VALUES="${GAMMA_VALUES}" \
RHO_VALUES="${RHO_VALUES}" \
REPS="${REPS}" \
TIME_LIMIT="${TIME_LIMIT}" \
WORKERS="${WORKERS}" \
DUE_REFERENCE_GAMMA="${DUE_REFERENCE_GAMMA}" \
"${SCRIPT_DIR}/run_sensitivity_partial_full.py"
