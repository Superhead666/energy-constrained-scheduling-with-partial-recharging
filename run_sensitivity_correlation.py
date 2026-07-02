#!/usr/bin/env python3
import csv
import json
import math
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
INSTANCE_DIR = SCRIPT_DIR / "instances" / "sensitivity_correlation"

CORRELATIONS = [x for x in os.environ.get("CORRELATIONS", "positive,independent,negative").split(",") if x]
REPS = int(os.environ.get("REPS", "20"))
WORKERS = int(os.environ.get("WORKERS", "4"))

TAU = os.environ.get("TAU", "0.5")
R = os.environ.get("R", "0.4")
GAMMA = os.environ.get("GAMMA", "4")
DUE_REFERENCE_GAMMA = os.environ.get("DUE_REFERENCE_GAMMA", "4")
DUE_REFERENCE_CORRELATION = os.environ.get("DUE_REFERENCE_CORRELATION", "independent")
P_MIN = os.environ.get("P_MIN", "2")
P_MAX = os.environ.get("P_MAX", "30")
E_MIN = os.environ.get("E_MIN", "2")
E_MAX = os.environ.get("E_MAX", "50")
ALPHA = os.environ.get("ALPHA", "10")
BETA = os.environ.get("BETA", "1")

N_EXACT = int(os.environ.get("N_EXACT", "50"))
N_ROBUST = int(os.environ.get("N_ROBUST", "100"))
TIME_EXACT = os.environ.get("TIME_EXACT", "600")
TIME_HEURISTIC = os.environ.get("TIME_HEURISTIC", "30")

RAW_OUT = RESULTS_DIR / "sensitivity_correlation_raw.csv"
SUMMARY_OUT = RESULTS_DIR / "sensitivity_correlation_summary.csv"
LOG_FILE = LOG_DIR / "sensitivity_correlation.log"


COMMON_ARGS = [
    "--mode", "one",
    "--due-mode", "taur",
    "--tau", TAU,
    "--R", R,
    "--gamma", GAMMA,
    "--due-reference-gamma", DUE_REFERENCE_GAMMA,
    "--due-reference-correlation", DUE_REFERENCE_CORRELATION,
    "--common-base-seeds",
    "--p-min", P_MIN,
    "--p-max", P_MAX,
    "--e-min", E_MIN,
    "--e-max", E_MAX,
    "--alpha", ALPHA,
    "--beta", BETA,
]


def pearson_from_jobs(jobs):
    ps = [job["p"] for job in jobs]
    es = [job["e"] for job in jobs]
    mp = sum(ps) / len(ps)
    me = sum(es) / len(es)
    cov = sum((p - mp) * (e - me) for p, e in zip(ps, es))
    sp = math.sqrt(sum((p - mp) ** 2 for p in ps))
    se = math.sqrt(sum((e - me) ** 2 for e in es))
    if sp == 0.0 or se == 0.0:
        return 0.0
    return cov / (sp * se)


def export_instance(n, corr, rep):
    out_dir = INSTANCE_DIR / f"n{n}" / corr
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"n{n}_tau{TAU.replace('.', 'p')}_R{R.replace('.', 'p')}_gamma{GAMMA.replace('.', 'p')}_{corr}_rep{rep}.json"
    cmd = [
        str(BIN),
        "--mode", "export-instance",
        *COMMON_ARGS[2:],
        "--n", str(n),
        "--rep", str(rep),
        "--pe-correlation", corr,
        "--out", str(path),
    ]
    subprocess.run(cmd, check=True)
    with path.open() as f:
        data = json.load(f)
    return pearson_from_jobs(data["jobs"])


def run_case(size_label, n, corr, rep):
    if size_label == "exact":
        solver = "milp-heur"
        time_limit = TIME_EXACT
    else:
        solver = "grasp-pr"
        time_limit = TIME_HEURISTIC
    cmd = [
        str(BIN),
        *COMMON_ARGS,
        "--solver", solver,
        "--n", str(n),
        "--rep", str(rep),
        "--time", str(time_limit),
        "--pe-correlation", corr,
    ]
    completed = subprocess.run(cmd, check=True, text=True, capture_output=True)
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    row = next(csv.DictReader(lines))
    row["size_label"] = size_label
    row["correlation_label"] = corr
    return row


def summarize(rows, pearson_map):
    groups = {}
    for row in rows:
        key = (row["size_label"], row["correlation_label"])
        groups.setdefault(key, []).append(row)

    fields = [
        "size_label",
        "n",
        "correlation_label",
        "instances",
        "avg_pearson",
        "optimal_instances",
        "avg_obj",
        "avg_dev",
        "max_dev",
        "avg_recharges",
        "avg_battery_utilization",
        "avg_makespan",
        "avg_time",
    ]
    with SUMMARY_OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for key in sorted(groups, key=lambda x: (0 if x[0] == "exact" else 1, CORRELATIONS.index(x[1]))):
            size_label, corr = key
            group = groups[key]
            avg = lambda name: sum(float(row[name]) for row in group) / len(group)
            dev = [max(0.0, float(row["objective"]) - float(row["best_bound"])) for row in group]
            n = int(float(group[0]["n"]))
            pearsons = [pearson_map[(n, corr, int(float(row["rep"])))] for row in group]
            writer.writerow({
                "size_label": size_label,
                "n": n,
                "correlation_label": corr,
                "instances": len(group),
                "avg_pearson": f"{sum(pearsons) / len(pearsons):.3f}",
                "optimal_instances": sum(int(float(row["certified"])) for row in group),
                "avg_obj": f"{avg('objective'):.2f}",
                "avg_dev": f"{sum(dev) / len(dev):.2f}",
                "max_dev": f"{max(dev):.2f}",
                "avg_recharges": f"{avg('recharges'):.2f}",
                "avg_battery_utilization": f"{avg('battery_utilization'):.3f}",
                "avg_makespan": f"{avg('makespan'):.2f}",
                "avg_time": f"{avg('runtime'):.2f}",
            })


def main():
    if not BIN.exists():
        subprocess.run([str(SCRIPT_DIR / "build_cplex.sh")], check=True)
    RESULTS_DIR.mkdir(exist_ok=True)
    LOG_DIR.mkdir(exist_ok=True)
    INSTANCE_DIR.mkdir(parents=True, exist_ok=True)

    tasks = []
    pearson_map = {}
    for n, size_label in [(N_EXACT, "exact"), (N_ROBUST, "robust")]:
        for corr in CORRELATIONS:
            for rep in range(REPS):
                pearson_map[(n, corr, rep)] = export_instance(n, corr, rep)
                tasks.append((size_label, n, corr, rep))

    LOG_FILE.write_text(
        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting processing-energy correlation sensitivity\n"
        f"Exact: n={N_EXACT}; solver=MILP-SH; time={TIME_EXACT}s\n"
        f"Robust: n={N_ROBUST}; solver=GRASP-PR; time={TIME_HEURISTIC}s\n"
        f"Common grid: tau={TAU}; R={R}; gamma={GAMMA}; due-reference-gamma={DUE_REFERENCE_GAMMA}; due-reference-correlation={DUE_REFERENCE_CORRELATION}; reps={REPS}; correlations={CORRELATIONS}\n"
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
            key=lambda r: (
                0 if r["size_label"] == "exact" else 1,
                CORRELATIONS.index(r["correlation_label"]),
                int(float(r["rep"]))
            )
        ))
    summarize(rows, pearson_map)
    with LOG_FILE.open("a") as f:
        f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Finished processing-energy correlation sensitivity\n")
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
