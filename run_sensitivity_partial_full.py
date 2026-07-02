#!/usr/bin/env python3
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
BIN = SCRIPT_DIR / "bin" / "caie_part1_final"
RESULTS_DIR = SCRIPT_DIR / "results"
LOG_DIR = SCRIPT_DIR / "logs"
INSTANCE_DIR = SCRIPT_DIR / "instances" / "sensitivity_partial_full"


def parse_list(text, cast=float):
    return [cast(x) for x in text.split(",") if x.strip()]


N = int(os.environ.get("N", "40"))
TAU = os.environ.get("TAU", "0.5")
R = os.environ.get("R", "0.4")
GAMMA_VALUES = parse_list(os.environ.get("GAMMA_VALUES", "3,5"), float)
RHO_VALUES = parse_list(os.environ.get("RHO_VALUES", "0.2,0.6"), float)
REPS = int(os.environ.get("REPS", "20"))
TIME_LIMIT = os.environ.get("TIME_LIMIT", "180")
WORKERS = int(os.environ.get("WORKERS", "4"))
P_MIN = os.environ.get("P_MIN", "2")
P_MAX = os.environ.get("P_MAX", "30")
E_MIN = os.environ.get("E_MIN", "2")
E_MAX = os.environ.get("E_MAX", "50")
BETA = os.environ.get("BETA", "1")
DUE_REFERENCE_GAMMA = os.environ.get("DUE_REFERENCE_GAMMA", "4")

RAW_OUT = RESULTS_DIR / "sensitivity_partial_full_exact_n40_raw.csv"
SUMMARY_OUT = RESULTS_DIR / "sensitivity_partial_full_exact_n40_summary.csv"
LOG_FILE = LOG_DIR / "sensitivity_partial_full.log"


COMMON_ARGS = [
    "--mode", "one",
    "--due-mode", "taur",
    "--n", str(N),
    "--tau", TAU,
    "--R", R,
    "--due-reference-gamma", DUE_REFERENCE_GAMMA,
    "--common-base-seeds",
    "--p-min", P_MIN,
    "--p-max", P_MAX,
    "--e-min", E_MIN,
    "--e-max", E_MAX,
    "--beta", BETA,
    "--time", TIME_LIMIT,
]


def gamma_label(value):
    return str(value).replace(".", "p")


def rho_label(value):
    return str(value).replace(".", "p")


def run_case(policy, gamma, rho, rep):
    solver = "milp-heur" if policy == "partial" else "milp-full-heur"
    cmd = [
        str(BIN),
        *COMMON_ARGS,
        "--solver", solver,
        "--gamma", str(gamma),
        "--rho", str(rho),
        "--rep", str(rep),
    ]
    completed = subprocess.run(cmd, check=True, text=True, capture_output=True)
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"No CSV row returned for {policy}, gamma={gamma}, rho={rho}, rep={rep}")
    row = next(csv.DictReader(lines))
    row["policy"] = policy
    return row


def export_instance(gamma, rho, rep):
    out_dir = INSTANCE_DIR / f"gamma{gamma_label(gamma)}_rho{rho_label(rho)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"n{N}_tau{TAU.replace('.', 'p')}_R{R.replace('.', 'p')}_gamma{gamma_label(gamma)}_rho{rho_label(rho)}_rep{rep}.json"
    cmd = [
        str(BIN),
        "--mode", "export-instance",
        "--due-mode", "taur",
        "--n", str(N),
        "--tau", TAU,
        "--R", R,
        "--gamma", str(gamma),
        "--rho", str(rho),
        "--due-reference-gamma", DUE_REFERENCE_GAMMA,
        "--common-base-seeds",
        "--p-min", P_MIN,
        "--p-max", P_MAX,
        "--e-min", E_MIN,
        "--e-max", E_MAX,
        "--beta", BETA,
        "--rep", str(rep),
        "--out", str(path),
    ]
    subprocess.run(cmd, check=True)
    return path


def to_float(row, key):
    return float(row[key])


def summarize(rows):
    by_key = {}
    for row in rows:
        key = (int(float(row["n"])), float(row["theta"]), float(row["R"]),
               float(row["sigma"]), float(row["rho"]), int(float(row["rep"])))
        by_key.setdefault(key, {})[row["policy"]] = row

    fields = [
        "gamma",
        "rho",
        "instances",
        "opt_partial",
        "opt_full",
        "U_partial",
        "U_full",
        "delta_U",
        "late_reduction_pct",
        "recharges_partial",
        "recharges_full",
        "util_partial",
        "util_full",
        "delta_recharge_time",
        "predicted_delta_recharge_time",
        "delta_makespan",
    ]
    groups = {}
    for key, pair in by_key.items():
        if "partial" not in pair or "full" not in pair:
            continue
        _, _, _, gamma, rho, _ = key
        groups.setdefault((gamma, rho), []).append(pair)

    with SUMMARY_OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for (gamma, rho), pairs in sorted(groups.items()):
            def avg(policy, name):
                return sum(to_float(pair[policy], name) for pair in pairs) / len(pairs)

            full_u = avg("full", "objective")
            partial_u = avg("partial", "objective")
            delta_u = full_u - partial_u
            late_reduction = 0.0 if full_u <= 1e-9 else 100.0 * delta_u / full_u
            pred = sum(
                to_float(pair["full"], "beta") *
                (to_float(pair["full"], "recharges") * to_float(pair["full"], "capacity") -
                 to_float(pair["full"], "sum_e"))
                for pair in pairs
            ) / len(pairs)
            writer.writerow({
                "gamma": f"{gamma:.1f}",
                "rho": f"{rho:.2f}",
                "instances": len(pairs),
                "opt_partial": sum(int(float(pair["partial"]["certified"])) for pair in pairs),
                "opt_full": sum(int(float(pair["full"]["certified"])) for pair in pairs),
                "U_partial": f"{partial_u:.2f}",
                "U_full": f"{full_u:.2f}",
                "delta_U": f"{delta_u:.2f}",
                "late_reduction_pct": f"{late_reduction:.2f}",
                "recharges_partial": f"{avg('partial', 'recharges'):.2f}",
                "recharges_full": f"{avg('full', 'recharges'):.2f}",
                "util_partial": f"{avg('partial', 'battery_utilization'):.3f}",
                "util_full": f"{avg('full', 'battery_utilization'):.3f}",
                "delta_recharge_time": f"{avg('full', 'total_recharge_time') - avg('partial', 'total_recharge_time'):.2f}",
                "predicted_delta_recharge_time": f"{pred:.2f}",
                "delta_makespan": f"{avg('full', 'makespan') - avg('partial', 'makespan'):.2f}",
            })


def main():
    if not BIN.exists():
        subprocess.run([str(SCRIPT_DIR / "build_cplex.sh")], check=True)
    RESULTS_DIR.mkdir(exist_ok=True)
    LOG_DIR.mkdir(exist_ok=True)
    INSTANCE_DIR.mkdir(parents=True, exist_ok=True)

    tasks = []
    for gamma in GAMMA_VALUES:
        for rho in RHO_VALUES:
            for rep in range(REPS):
                export_instance(gamma, rho, rep)
                tasks.append(("partial", gamma, rho, rep))
                tasks.append(("full", gamma, rho, rep))

    started = time.strftime("%Y-%m-%d %H:%M:%S")
    LOG_FILE.write_text(
        f"[{started}] Starting partial-vs-full sensitivity\n"
        f"Grid: n={N}; tau={TAU}; R={R}; gamma={GAMMA_VALUES}; rho={RHO_VALUES}; reps={REPS}; time={TIME_LIMIT}s\n"
        f"Output: {RAW_OUT}\n"
    )

    rows = []
    done = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = [pool.submit(run_case, *task) for task in tasks]
        for future in as_completed(futures):
            rows.append(future.result())
            done += 1
            if done % 10 == 0 or done == len(tasks):
                msg = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Completed {done}/{len(tasks)} solves\n"
                with LOG_FILE.open("a") as f:
                    f.write(msg)
                print(msg, end="")

    fieldnames = list(rows[0].keys())
    with RAW_OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted(
            rows,
            key=lambda r: (float(r["sigma"]), float(r["rho"]), int(float(r["rep"])), r["policy"])
        ))
    summarize(rows)
    with LOG_FILE.open("a") as f:
        f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Finished partial-vs-full sensitivity\n")
        f.write(f"Summary: {SUMMARY_OUT}\n")
    print(f"Wrote {RAW_OUT}")
    print(f"Wrote {SUMMARY_OUT}")
    print(f"Wrote instances under {INSTANCE_DIR}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
