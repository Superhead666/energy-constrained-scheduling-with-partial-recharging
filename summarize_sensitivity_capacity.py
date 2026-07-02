#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"
MANUSCRIPT_REPS = set(range(18)) | {19, 20}


def read_rows(path):
    rows = []
    with path.open(newline="") as f:
        for raw in csv.DictReader(f):
            row = {}
            for key, value in raw.items():
                if key == "algorithm":
                    row[key] = value
                else:
                    row[key] = float(value)
            rows.append(row)
    return rows


def write_filtered_rows(path, out_path, reps):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
        fields = rows[0].keys()
    rows = [row for row in rows if int(float(row["rep"])) in reps]
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summarize(path, out_path, reps=None):
    groups = defaultdict(list)
    for row in read_rows(path):
        if reps is not None and int(row["rep"]) not in reps:
            continue
        groups[row["sigma"]].append(row)

    fields = [
        "gamma",
        "instances",
        "optimal_instances",
        "avg_obj",
        "avg_dev",
        "max_dev",
        "avg_recharges",
        "avg_battery_utilization",
        "avg_makespan",
        "avg_time",
    ]
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for gamma in sorted(groups):
            rows = groups[gamma]
            dev = [max(0.0, row["objective"] - row["best_bound"]) for row in rows]
            avg = lambda name: sum(row[name] for row in rows) / len(rows)
            writer.writerow({
                "gamma": f"{gamma:.1f}",
                "instances": len(rows),
                "optimal_instances": sum(int(row["certified"]) for row in rows),
                "avg_obj": f"{avg('objective'):.2f}",
                "avg_dev": f"{sum(dev) / len(dev):.2f}",
                "max_dev": f"{max(dev):.2f}",
                "avg_recharges": f"{avg('recharges'):.2f}",
                "avg_battery_utilization": f"{avg('battery_utilization'):.3f}",
                "avg_makespan": f"{avg('makespan'):.2f}",
                "avg_time": f"{avg('runtime'):.2f}",
            })


def main():
    exact = RESULTS_DIR / "sensitivity_capacity_gamma_milp_sh_n50.csv"
    robust = RESULTS_DIR / "sensitivity_capacity_gamma_grasp_pr_n100.csv"
    exact_out = RESULTS_DIR / "sensitivity_capacity_gamma_milp_sh_n50_summary.csv"
    robust_out = RESULTS_DIR / "sensitivity_capacity_gamma_grasp_pr_n100_summary.csv"
    summarize(exact, exact_out)
    summarize(robust, robust_out)
    exact_rows20 = RESULTS_DIR / "sensitivity_capacity_gamma_milp_sh_n50_manuscript20_rows.csv"
    robust_rows20 = RESULTS_DIR / "sensitivity_capacity_gamma_grasp_pr_n100_manuscript20_rows.csv"
    exact_out20 = RESULTS_DIR / "sensitivity_capacity_gamma_milp_sh_n50_manuscript20_summary.csv"
    robust_out20 = RESULTS_DIR / "sensitivity_capacity_gamma_grasp_pr_n100_manuscript20_summary.csv"
    write_filtered_rows(exact, exact_rows20, MANUSCRIPT_REPS)
    write_filtered_rows(robust, robust_rows20, MANUSCRIPT_REPS)
    summarize(exact, exact_out20, MANUSCRIPT_REPS)
    summarize(robust, robust_out20, MANUSCRIPT_REPS)
    print(f"Wrote {exact_out}")
    print(f"Wrote {robust_out}")
    print(f"Wrote {exact_rows20}")
    print(f"Wrote {robust_rows20}")
    print(f"Wrote {exact_out20}")
    print(f"Wrote {robust_out20}")


if __name__ == "__main__":
    main()
