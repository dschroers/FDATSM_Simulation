# FDATSM Simulation

## Description

This repository contains the code for the simulation study and empirical study presented in the article ["Dynamically Consistent Analysis of Realized Covariations in Term Structure Models"](https://arxiv.org/abs/2406.19412) by Dennis Schroers. The core methodology is implemented and documented in the [`FDATSM`](https://github.com/dschroers/FDATSM) R package.

The simulation study investigates the finite-sample performance of the truncated covariation estimator for difference returns, as outlined in the article. It aims to provide insights into the robustness of the estimator against jumps and its ability to recover the dimension of the original quadratic covariation.

---

## Files

### Simulation Study
- `Simulation Study.R`: Main R script for replicating the simulation study.
- `simulation_results.csv`: Output of the simulation study (500 Monte Carlo iterations).

### Simulation Study with Noise-Volatility Dependence
- `Simulation Study with noise-volatility dep.R`: Simulation study with noise-volatility dependence (as in the appendix).
- `simulation_results_with_dep.csv`: Results of this simulation study (500 Monte Carlo iterations).

### Empirical Study
- `Empirical Study (Sec 5.1,5.2).R`: Reproduces empirical results for Sections 5.1 and 5.2 of the article.
- `Empirical Study (Sec 5.3).R`: Reproduces empirical results for Section 5.3.
- `empirical_results_jumps.csv`: Output for Section 5.1 (jump analysis).
- `empirical_results_dimensions.csv`: Output for Section 5.2 (dimension estimation).
- `empirical_results_outofsample.csv`: Output for Section 5.3 (out-of-sample analysis).
---

## Installation

To run the code, you must have the following R packages installed:

- `MASS`
- `matrixcalc`
- `pracma`
- `npreg`
- `splines`

Install them using:

```r
install.packages(c("MASS", "matrixcalc", "pracma", "npreg", "splines"))
```
