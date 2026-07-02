#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Launching MILP-S, MILP-SH-light, and GRASP-PR in parallel."

"${SCRIPT_DIR}/run_milp_s.sh" &
pid_milp_s=$!

"${SCRIPT_DIR}/run_milp_sh_light_start.sh" &
pid_milp_sh=$!

"${SCRIPT_DIR}/run_grasp_pr.sh" &
pid_grasp=$!

status=0
for item in \
  "MILP-S:${pid_milp_s}" \
  "MILP-SH-light:${pid_milp_sh}" \
  "GRASP-PR:${pid_grasp}"; do
  label="${item%%:*}"
  pid="${item##*:}"
  if wait "${pid}"; then
    echo "${label} finished successfully."
  else
    echo "${label} failed." >&2
    status=1
  fi
done

exit "${status}"
