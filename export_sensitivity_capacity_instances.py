#!/usr/bin/env python3
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
BIN = SCRIPT_DIR / "bin" / "caie_part1_final"
OUT_ROOT = SCRIPT_DIR / "instances" / "sensitivity_capacity_gamma"

GAMMA_VALUES = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6]
REPS = range(30)

COMMON_ARGS = [
    "--mode", "export-instance",
    "--due-mode", "taur",
    "--tau", "0.5",
    "--R", "0.4",
    "--due-reference-gamma", "4",
    "--common-base-seeds",
    "--p-min", "2",
    "--p-max", "30",
    "--e-min", "2",
    "--e-max", "50",
    "--alpha", "10",
    "--beta", "1",
]


def gamma_label(value):
    return str(value).replace(".", "p")


def export_group(n, label):
    out_dir = OUT_ROOT / label
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for gamma in GAMMA_VALUES:
        for rep in REPS:
            name = f"n{n}_tau0p5_R0p4_gamma{gamma_label(gamma)}_rep{rep}.json"
            path = out_dir / name
            cmd = [
                str(BIN),
                *COMMON_ARGS,
                "--n", str(n),
                "--gamma", str(gamma),
                "--rep", str(rep),
                "--out", str(path),
            ]
            subprocess.run(cmd, check=True)
            manifest.append(name)
    (out_dir / "manifest.txt").write_text("\n".join(manifest) + "\n")


def main():
    if not BIN.exists():
        subprocess.run([str(SCRIPT_DIR / "build_cplex.sh")], check=True)
    export_group(50, "exact_milp_sh_n50")
    export_group(100, "robustness_grasp_pr_n100")
    print(f"Wrote sensitivity instances under {OUT_ROOT}")


if __name__ == "__main__":
    main()
