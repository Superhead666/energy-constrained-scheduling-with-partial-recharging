#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${SCRIPT_DIR}/bin/caie_part1_final"

if [[ ! -x "${BIN}" ]]; then
  "${SCRIPT_DIR}/build_cplex.sh" >/dev/null
fi

mkdir -p "${SCRIPT_DIR}/tuning_results" "${SCRIPT_DIR}/logs"

# Calibration subset: it includes the large instances where the heuristic matters,
# but keeps tuning cheap enough to repeat.
N_VALUES="${N_VALUES:-100,200}"
TAU_VALUES="${TAU_VALUES:-0.5,0.6}"
R_VALUES="${R_VALUES:-0.8}"
GAMMA_VALUES="${GAMMA_VALUES:-3,5}"
REPS="${REPS:-3}"

P_MIN="${P_MIN:-2}"
P_MAX="${P_MAX:-30}"
E_MIN="${E_MIN:-2}"
E_MAX="${E_MAX:-50}"
ALPHA="${ALPHA:-10}"
BETA="${BETA:-1}"
TIME_LIMIT="${TIME_LIMIT:-30}"

configs=(
  "eps015_r075_e10 0.15 0.75 10"
  "eps010_r075_e10 0.10 0.75 10"
  "eps020_r075_e10 0.20 0.75 10"
  "eps015_r050_e10 0.15 0.50 10"
  "eps015_r100_e10 0.15 1.00 10"
  "eps015_r075_e15 0.15 0.75 15"
)

for cfg in "${configs[@]}"; do
  read -r tag eps radius elite <<<"${cfg}"
  out="${SCRIPT_DIR}/tuning_results/${tag}.csv"
  log="${SCRIPT_DIR}/logs/tuning_${tag}.log"
  {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Tuning ${tag}"
    echo "eps=${eps}; radius=${radius}; elite=${elite}; subset n=${N_VALUES}; reps=${REPS}"
  } | tee "${log}"
  GRASP_EPS="${eps}" \
  GRASP_LS_RADIUS_FACTOR="${radius}" \
  GRASP_ELITE_SIZE="${elite}" \
  "${BIN}" \
    --mode part1 \
    --solver grasp-pr \
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
    --out "${out}" 2>&1 | tee -a "${log}"
done

"${SCRIPT_DIR}/summarize_results.py" "${SCRIPT_DIR}/tuning_results"/*.csv | tee "${SCRIPT_DIR}/tuning_results/summary.txt"
