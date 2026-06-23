# Multivariate Probit Regression Using Antedependence Structures

Bayesian multivariate probit models for correlated binary longitudinal data, with an emphasis on first-order antedependence structures. The repository contains the Stan models, simulation studies, real-data analyses, and figures associated with the manuscript **“Bayesian Multivariate Probit Regression for Longitudinal Binary Data via Structured Latent Dependence Models.”**

## Overview

For subject \(i\) at occasion \(t\), the model introduces a latent Gaussian response

\[
Z_{it} = x_{it}^{\mathsf T}\beta + \varepsilon_{it},
\qquad
Y_{it} = \mathbb{1}(Z_{it} > 0),
\]

where the within-subject errors are jointly normal and their correlation matrix describes the dependence among repeated binary outcomes.

The implementations cover the following residual correlation structures:

- independence;
- AR(1);
- Toeplitz (lag-specific correlations);
- first-order antedependence, AD(1);
- ARMA(1,1); and
- an unstructured correlation matrix.

Under AD(1), correlations are built from local dependence parameters. For \(s<t\),

\[
\operatorname{Corr}(Z_{is},Z_{it}) = \prod_{r=s}^{t-1}\rho_r.
\]

This parameterization is useful when dependence evolves across unequally related consecutive occasions rather than being governed by one stationary correlation coefficient.

## Repository structure

```text
.
├── Programs/
│   ├── Aux_Functions.R
│   ├── FitMultProbit_*.stan
│   ├── FitProbit_Ind*.stan
│   └── simMultProbit*.stan
├── SimulationStudies/
│   ├── simLatinSquare.R
│   ├── RunRepStudy_MultProbit.R
│   ├── CompRepStudy.R
│   ├── RunSensStudy_MultProbit.R
│   └── CompSensStudy_MultProbit.R
├── RealDataAnalysis/
│   ├── Data/
│   ├── Madras/
│   ├── RealDataAnal_HPC.R
│   ├── RealDataAnal_Comp.R
│   └── RealDataAnal_PostSum.R
└── Figures/
```

### Stan programs

The files in `Programs/` are organized by model and data pattern:

- `FitMultProbit_Structured.stan` implements Toeplitz, AR(1), AD(1), and ARMA(1,1), selected through `cor_type`;
- `FitMultProbit_Lcorr.stan` fits an unstructured correlation matrix;
- `FitMultProbit_AR1.stan` and `FitMultProbit_Toep.stan` provide structure-specific implementations;
- `FitProbit_Ind.stan` fits the independence model;
- files containing `_Missing` support outcomes encoded as `-1` for missing values;
- files containing `_Dropout` target monotone dropout patterns; and
- files containing `_Diag` include additional diagnostic quantities.

The `simMultProbit*.stan` files generate multivariate probit data. Shared R utilities for correlation matrices, initialization, parameter selection, data splitting, and posterior calculations are defined in `Programs/Aux_Functions.R`.

## Requirements

- R with a working C++ toolchain;
- Stan through the `rstan` package; and
- the R packages used by the analysis scripts.

Install the main dependencies with:

```r
install.packages(c(
  "rstan", "here", "dplyr", "tidyr", "ggplot2", "patchwork",
  "loo", "knitr", "xtable", "shinystan", "boa", "remotes"
))

remotes::install_github("ewan-keith/rstansim")
```

Follow the platform-specific `rstan` installation instructions if Stan models cannot be compiled. Model fitting can be computationally intensive, particularly for the unstructured model and scenarios with many repeated measurements.

## Getting started

Clone the repository and start R from its root directory:

```bash
git clone https://github.com/silvadosjr/Multivariate-Probit-regression-antedependence.git
cd Multivariate-Probit-regression-antedependence
```

Then load the shared functions:

```r
library(rstan)
library(here)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source(here("Programs", "Aux_Functions.R"))
```

Several scripts currently contain the original development path:

```r
setwd("~/GitHub/Multivariate-Probit-regression-antedependence/")
```

Before running them, either launch R from the repository root and remove/comment out that line, or replace it with the location of your clone. The legacy exploratory script `RealDataAnalysis/Main.R` contains additional machine-specific absolute paths; the `RealDataAnal_*.R` scripts are the recommended entry points.

## Simulation studies

The simulation workflow has three stages.

### 1. Generate data

Open `SimulationStudies/simLatinSquare.R`, choose the generating structure,

```r
corr_type <- "AR1"       # or "Toeplitz"
```

review the scenario matrix and number of Monte Carlo replications, and run:

```r
source("SimulationStudies/simLatinSquare.R")
```

The scenarios vary the number of repeated measurements, sample size, and correlation strength.

> **Path note:** `simLatinSquare.R` currently writes to `LatinSquareScenarios/<corr_type>` at the repository root, whereas the fitting scripts read from `SimulationStudies/LatinSquareScenarios/<corr_type>`. Make these paths consistent before starting a full run.

### 2. Fit competing models

For the robustness study, configure `corr_type` and `fit_models` in `RunRepStudy_MultProbit.R`, then run:

```r
source("SimulationStudies/RunRepStudy_MultProbit.R")
```

For the prior-sensitivity study, use:

```r
source("SimulationStudies/RunSensStudy_MultProbit.R")
```

Both scripts save posterior summaries as RDS files. Output directories must exist or be created before a run if the selected script does not create them automatically.

### 3. Summarize results

```r
source("SimulationStudies/CompRepStudy.R")
source("SimulationStudies/CompSensStudy_MultProbit.R")
```

The summaries include bias, absolute bias, RMSE, posterior standard deviation, effective sample size, \(\widehat R\), and Stein’s loss for estimated correlation matrices. Tables and EPS figures are produced for reporting.

## Real-data analyses

Two longitudinal binary datasets are included under `RealDataAnalysis/Data/`:

- **Six Cities:** 537 subjects measured at ages 7–10, with maternal smoking as a covariate;
- **Madras:** 12 monthly measurements, with age, gender, and time interactions; the corresponding workflow supports incomplete longitudinal responses.

### Six Cities

```r
source("RealDataAnalysis/RealDataAnal_HPC.R")
source("RealDataAnalysis/RealDataAnal_Comp.R")
source("RealDataAnalysis/RealDataAnal_PostSum.R")
```

The scripts fit independence, Toeplitz, AR(1), AD(1), ARMA(1,1), and unstructured models. They perform a subject-level train/test split, compute held-out predictive criteria, compare models with LOO and WAIC, and generate posterior summaries and figures.

### Madras

```r
source("RealDataAnalysis/Madras/RealDataAnal_HPC.R")
source("RealDataAnalysis/Madras/RealDataAnal_Comp.R")
source("RealDataAnalysis/Madras/RealDataAnal_PostSum.R")
```

This workflow compares independence, AD(1), and unstructured models using Stan programs that accommodate missing observations and monotone dropout.

The `*_HPC.R` filenames refer to **held-out predictive checks**, not to a required high-performance computing environment.

## Reproducibility notes

- The main scripts use seed `258` for train/test splitting where applicable.
- MCMC settings are declared near the top of each fitting script and should be reviewed before execution.
- Some study scripts default to one chain. Use at least two chains—preferably four—when \(\widehat R\) is required as a convergence diagnostic.
- Generated fits and simulation results can be large; RDS outputs are intended to be produced locally.
- Run scripts from the repository root so that `here()` and relative data paths resolve correctly.

## Citation

If you use this code in academic work, please cite the associated manuscript:

> Santos, J. R. S., Farias, R. B. A., & Santos-Neto, M. (2026). *Bayesian multivariate probit regression for longitudinal binary data via structured latent dependence models* [Manuscript in preparation for submission to the *Journal of Statistical Computation and Simulation*].

BibTeX:

```bibtex
@unpublished{dosSantos2026BayesianMultivariateProbit,
  author = {dos Santos, Jos\'{e} Roberto Silva and Farias, Rafael Br\'{a}z Azevedo and Santos-Neto, Manoel},
  title  = {Bayesian Multivariate Probit Regression for Longitudinal Binary Data via Structured Latent Dependence Models},
  year   = {2026},
  note   = {Manuscript in preparation for submission to the Journal of Statistical Computation and Simulation}
}
```

Update this entry with the DOI and publication details after acceptance.

## License

No license file is currently included. Unless a license is added, the source code remains under the copyright of its authors and reuse is not automatically granted.
