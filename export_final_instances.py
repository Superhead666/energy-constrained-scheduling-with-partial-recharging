#!/usr/bin/env python3
import csv
import json
import os
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
BIN = SCRIPT_DIR / "bin" / "caie_part1_final"
RESULT_FILE = SCRIPT_DIR / "results" / "final_milp_s.csv"
OUT_DIR = SCRIPT_DIR / "instances" / "part1"


def key_from_row(row):
    return (
        int(float(row["n"])),
        float(row["theta"]),
        float(row["R"]),
        float(row["sigma"]),
        int(float(row["rep"])),
    )


def fmt_float(value):
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return text.replace(".", "p")


def instance_name(n, tau, r_value, gamma, rep):
    return (
        f"n{n}_tau{fmt_float(tau)}_R{fmt_float(r_value)}_"
        f"gamma{fmt_float(gamma)}_rep{rep}.json"
    )


def main():
    if not BIN.exists():
        subprocess.run([str(SCRIPT_DIR / "build_cplex.sh")], check=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with RESULT_FILE.open(newline="") as f:
        rows = list(csv.DictReader(f))
    keys = sorted({key_from_row(row) for row in rows})

    manifest = []
    for n, tau, r_value, gamma, rep in keys:
        filename = instance_name(n, tau, r_value, gamma, rep)
        out_path = OUT_DIR / filename
        cmd = [
            str(BIN),
            "--mode", "export-instance",
            "--due-mode", "taur",
            "--n", str(n),
            "--tau", str(tau),
            "--R", str(r_value),
            "--gamma", str(gamma),
            "--p-min", "2",
            "--p-max", "30",
            "--e-min", "2",
            "--e-max", "50",
            "--alpha", "10",
            "--beta", "1",
            "--rep", str(rep),
            "--out", str(out_path),
        ]
        subprocess.run(cmd, check=True)
        manifest.append({
            "file": str(out_path.relative_to(SCRIPT_DIR)),
            "n": n,
            "tau": tau,
            "R": r_value,
            "gamma": gamma,
            "rep": rep,
        })

    manifest_path = OUT_DIR / "manifest.json"
    with manifest_path.open("w") as f:
        json.dump({
            "description": "Part 1 computational instances used in the CAIE revision experiments.",
            "source_results": str(RESULT_FILE.relative_to(SCRIPT_DIR)),
            "count": len(manifest),
            "instances": manifest,
        }, f, indent=2)
        f.write("\n")

    print(f"Exported {len(manifest)} instances to {OUT_DIR}")


if __name__ == "__main__":
    main()
