#!/usr/bin/env python3
import csv
import os
import sys
from collections import defaultdict


def read_rows(paths):
    rows = []
    for path in paths:
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                row["_file"] = os.path.basename(path)
                row["n"] = int(float(row["n"]))
                row["rep"] = int(float(row["rep"]))
                row["theta"] = float(row["theta"])
                row["R"] = float(row["R"])
                row["sigma"] = float(row["sigma"])
                row["objective"] = float(row["objective"])
                row["best_bound"] = float(row["best_bound"])
                row["certified"] = int(float(row["certified"]))
                row["runtime"] = float(row["runtime"])
                rows.append(row)
    return rows


def main():
    paths = sys.argv[1:]
    if not paths:
        paths = [
            "results/final_milp_s.csv",
            "results/final_milp_sh_light_start.csv",
            "results/final_grasp_pr.csv",
        ]
    rows = read_rows(paths)
    if not rows:
        raise SystemExit("No rows found.")

    key = lambda r: (r["n"], r["theta"], r["R"], r["sigma"], r["rep"])
    best = {}
    for r in rows:
        k = key(r)
        best[k] = min(best.get(k, r["objective"]), r["objective"])

    groups = defaultdict(list)
    for r in rows:
        groups[(r["_file"], r["algorithm"], r["n"])].append(r)
        groups[(r["_file"], r["algorithm"], "All")].append(r)

    print("file,algorithm,n,count,avg_obj,best_count,bk_gap_pct,lb_gap_pct,certified,avg_time")
    for (file_name, alg, n), rs in sorted(groups.items(), key=lambda x: (x[0][0], str(x[0][2]), x[0][1])):
        avg_obj = sum(r["objective"] for r in rs) / len(rs)
        best_count = sum(1 for r in rs if r["objective"] <= best[key(r)] + 1e-9)
        bk_gap = sum((r["objective"] - best[key(r)]) / max(1.0, r["objective"]) * 100.0 for r in rs) / len(rs)
        lb_gap = sum((r["objective"] - r["best_bound"]) / max(1.0, r["objective"]) * 100.0 for r in rs) / len(rs)
        cert = sum(r["certified"] for r in rs)
        avg_time = sum(r["runtime"] for r in rs) / len(rs)
        print(f"{file_name},{alg},{n},{len(rs)},{avg_obj:.3f},{best_count},{bk_gap:.3f},{lb_gap:.3f},{cert},{avg_time:.2f}")


if __name__ == "__main__":
    main()
