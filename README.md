# Energy-constrained single-machine scheduling to minimize the number of late jobs with partial recharge and generalized due dates computational package

This repository contains the code, generated instances, and computational results used in the revised manuscript titled ``Energy-constrained single-machine scheduling to minimize the number of late jobs with partial recharge and generalized due dates.''

## Compared methods

- `MILP-S`: strengthened position-indexed MILP.
- `MILP-SH`: MILP-S with the relaxation-guided DP warm start.
- `GRASP-PR`: GRASP construction, exact fixed-sequence DP evaluation, bound-guided local search, and path relinking.

## Instance settings

- \(n\in\{50,100,200\}\)
- \(\tau\in\{0.5,0.6\}\)
- \(R=0.8\)
- \(\gamma\in\{3,5\}\)
- \(p_j\sim U[2,30]\)
- \(e_j\sim U[2,50]\)
- \(\alpha=10\), \(\beta=1\)
- 10 instances for each parameter setting, 120 instances in total

The final instance set excludes the outlier replication
`n100_tau0p6_R0p8_gamma5_rep0` and uses the replacement replication
`n100_tau0p6_R0p8_gamma5_rep10`.

## Rebuilding

The code requires IBM ILOG CPLEX. The default local CPLEX path is
`/Applications/CPLEX_Studio2211`.

```bash
./build_cplex.sh
```

## Running the experiments

Run the three methods separately:

```bash
./run_milp_s.sh
./run_milp_sh_light_start.sh
./run_grasp_pr.sh
```

or run them in parallel:

```bash
./run_all_parallel.sh
```

The final CSV files are:

- `results/final_milp_s.csv`
- `results/final_milp_sh_light_start.csv`
- `results/final_grasp_pr.csv`
- `results/final_summary.csv`

To regenerate the summary:

```bash
./summarize_results.py results/final_milp_s.csv results/final_milp_sh_light_start.csv results/final_grasp_pr.csv > results/final_summary.csv
```

## GRASP-PR tuning

The final parameter values are:

- \(\epsilon=0.15\)
- local-search radius \(0.75\sqrt n\)
- elite pool size 15
- stabilization threshold \(\eta=25\)

The pilot tuning results are retained in `tuning_results/`.

To rerun the tuning subset:

```bash
./run_tuning_grasp_pr.sh
```

## Instance JSON files

The exact instances used in the final table are stored in:

```bash
instances/part1/
```

Each JSON file contains the processing time, energy consumption, sorted generalized due dates, capacity, recharge parameters, generation parameters, and random seed. The file `instances/part1/manifest.json` lists all 120 instances.

To regenerate the JSON files from the final result keys:

```bash
./export_final_instances.py
```

## Capacity sensitivity analysis

The exact capacity sensitivity experiment reported in the manuscript uses the following setting:

- \(n=50\)
- \(\tau=0.5\)
- \(R=0.4\)
- \(\gamma\in\{2,2.5,3,3.5,4,4.5,5,5.5,6\}\)
- due-date reference capacity fixed at \(\gamma=4\)
- \(p_j\sim U[2,30]\)
- \(e_j\sim U[2,50]\)
- \(\alpha=10\), \(\beta=1\)
- 30 generated instances for each value of \(\gamma\), 270 instances in total
- the manuscript table uses 20 replications for each value of \(\gamma\)
- 600-second time limit for MILP-SH

The exact sensitivity results are stored in:

```bash
results/sensitivity_capacity_gamma_milp_sh_n50.csv
results/sensitivity_capacity_gamma_milp_sh_n50_summary.csv
results/sensitivity_capacity_gamma_milp_sh_n50_manuscript20_rows.csv
results/sensitivity_capacity_gamma_milp_sh_n50_manuscript20_summary.csv
instances/sensitivity_capacity_gamma/exact_milp_sh_n50/
```

The additional large-scale check with GRASP-PR uses the same capacity grid with \(n=100\) and a 30-second time limit. These results are stored in:

```bash
results/sensitivity_capacity_gamma_grasp_pr_n100.csv
results/sensitivity_capacity_gamma_grasp_pr_n100_summary.csv
results/sensitivity_capacity_gamma_grasp_pr_n100_manuscript20_rows.csv
results/sensitivity_capacity_gamma_grasp_pr_n100_manuscript20_summary.csv
instances/sensitivity_capacity_gamma/robustness_grasp_pr_n100/
```

To rerun the sensitivity tests:

```bash
./run_sensitivity_capacity_gamma.sh
```

To regenerate the sensitivity summaries and JSON instance files:

```bash
./summarize_sensitivity_capacity.py
./export_sensitivity_capacity_instances.py
```

## Partial versus full recharge sensitivity analysis

The partial/full recharge sensitivity experiment reported in the manuscript uses the following setting:

- \(n=40\)
- \(\tau=0.5\)
- \(R=0.4\)
- \(\gamma\in\{3,5\}\)
- fixed full-charge share \(\rho=\alpha/(\alpha+\beta E)\in\{0.2,0.6\}\)
- due-date reference capacity fixed at \(\gamma=4\)
- \(p_j\sim U[2,30]\)
- \(e_j\sim U[2,50]\)
- \(\beta=1\), with \(\alpha=\rho\beta E/(1-\rho)\)
- 20 paired replications for each \((\gamma,\rho)\) setting

Both the partial-recharge and full-recharge models are solved exactly in the reported rows. 
```bash
results/sensitivity_partial_full_exact_n40_manuscript20_rows.csv
results/sensitivity_partial_full_exact_n40_manuscript20_summary.csv
instances/sensitivity_partial_full/
```

To rerun the partial/full sensitivity test:

```bash
./run_sensitivity_partial_full.sh
```

## Processing-energy correlation sensitivity analysis

The processing-energy correlation sensitivity experiment reported in the manuscript uses the following setting:

- \(n=50\) for the exact MILP-SH panel
- \(n=100\) for the GRASP-PR robustness panel
- \(\tau=0.5\)
- \(R=0.4\)
- \(\gamma=4\)
- due-date reference capacity fixed at \(\gamma=4\)
- due-date reference correlation fixed at independent pairing
- \(p_j\sim U[2,30]\)
- \(e_j\sim U[2,50]\)
- \(\alpha=10\), \(\beta=1\)
- 20 replications for each correlation setting

For each replication, the same marginal samples of processing times and energy consumptions are used across all three settings. Positive correlation pairs sorted processing times with sorted energy consumptions, independent correlation uses a fixed random permutation of the sorted energy consumptions, and negative correlation pairs sorted processing times with reverse-sorted energy consumptions.

The results and exported instances are stored in:

```bash
results/sensitivity_correlation_raw.csv
results/sensitivity_correlation_summary.csv
instances/sensitivity_correlation/
```

To rerun the correlation sensitivity test:

```bash
./run_sensitivity_correlation.sh
```
