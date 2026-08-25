--Downward Prior Drift in Sequential Bayes Factor Designs--

This repository contains the data, code, and supplementary materials for the manuscript:
Downward Prior Drift in Sequential Bayes Factor Designs: A Monte Carlo Investigation of Error Rates, Power, and Threshold Calibration
The study examines how a systematic narrowing of the Cauchy prior scale over the course of a sequential Bayes factor (SBF) study affects the design's operating characteristics, and provides a parametric correction model for recalibrating Bayes factor thresholds under downward prior drift.

--Repository Structure--
downward_drift_bf/
README.md                          -- This file
codes/
01_sim_engine.R                    -- Main Monte Carlo simulation engine
02_analysis.R                      -- Analysis, tables, and figures
03_adaptive_threshold.R            -- BF threshold calibration via bisection
tables/
table1_error_rates_JZS.csv         -- Error rates for JZS default (r0 = 0.707)
table2_all_error_rates.csv         -- Full condition-level error rates
table3_calibrated_thresholds.csv   -- Calibrated BF* thresholds and CFs
figures/
fig1_typeI.png                     -- Type I error rate by drift magnitude
fig2_power.png                     -- Statistical power by drift magnitude
fig3_calibration.png               -- Calibrated BF* thresholds
results/
error_rates.rds                    -- Condition-level summaries (Type I/II error, power, inconclusive rate, median stopping n)
all_simulations.rds                -- Raw replication-level data (decision, stopping n, terminal BF, prior scale at stop)
correction_model_downward.rds      -- Fitted reduced OLS model object (log(CF) ~ delta r + r0c + interaction)

--Reproducing the Results--
Run the scripts in order:
1. 01_sim_engine.R -- Simulates 81 conditions (K = 2,000 reps each).
2. 02_analysis.R -- Summarizes results and generates tables/figures.
3. 03_adaptive_threshold.R -- Bisection search for BF threshold calibration. 

All scripts are deterministic (seeded by condition index) and write checkpoint files per condition, allowing safe resume after interruption.

--Prerequisites--
R >= 4.6

The following R packages:
    install.packages(c(
      "BayesFactor",
      "doParallel",
      "foreach",
      "dplyr",
      "tidyr",
      "ggplot2",
      "scales",
      "car"
    ))

--Notes--
- Deterministic seeds: base 20,240,602 (simulation), 99,999 (calibration). Exact reproduction requires the same R and BayesFactor versions.
- Design imbalance: Due to the feasibility constraint (r0 - delta r >= 0.1), not all r0 by delta r combinations are present. The correction model should not be extrapolated beyond the observed grid.

--Contact--
For questions about the code or the analysis, please open an issue in this repository or contact the corresponding author.

--License--
This project is licensed under the MIT License. The manuscript text and figures are subject to the journal's copyright upon publication.
