#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

N_VALUES="50,100,200" \
TAU_VALUES="0.5,0.6" \
R_VALUES="0.8" \
GAMMA_VALUES="3,5" \
P_MIN="2" \
P_MAX="30" \
E_MIN="2" \
E_MAX="50" \
ALPHA="10" \
BETA="1" \
REPS="10" \
TIME_LIMIT="600" \
"${SCRIPT_DIR}/run_method.sh" milp-s
