# Competition and Investment Model of Wealth Distribution

Replication code and data for:

> Ono, Y., & Ishida, A. ([YEAR]). Competition and investment model of wealth distribution. *The Journal of Mathematical Sociology*. [DOI]

## Files

### Julia code (`julia_simulation/`)

| File | Description | Figures |
| --- | --- | --- |
| `code/01_run_simulation.ipynb` | Simulation of the unified model: mean and Gini coefficient, time evolution of the Gini coefficient, and generation of the wealth distributions used in the analysis | Figures 1, 2 |

### R code (`r_analysis/script/`)

| File | Description | Figures |
| --- | --- | --- |
| `02_Fitting.R` | CCDF plots with gamma and log-normal fits | Figures 3, 4 |
| `03_Goodness_of_Fit.R` | Vuong tests: power-law vs. log-normal for the upper tail, gamma vs. log-normal for the bulk (about 3 minutes) | Figures 5, 6 |
| `04_Powerlaw_Bootstrap.R` | Power-law goodness-of-fit test with bootstrap p-values (about 1 hour) | Section 4 (text) |
| `05_Plot_Heatmaps.R` | Heatmaps of the test results | Figures 5, 6 |

### Data (`data/`)

| File | Description |
| --- | --- |
| `sorted_mean_incomes_99.feather` | Simulated wealth distributions for each (ω, α) |
| `sorted_mean_incomes_99.rda` | The same data in R format (`mean_income`), created by `02_Fitting.R` and used by `03` and `04` |
| `vuong_pvalue.rda` | Vuong test results (Figures 5 and 6), created by `03_Goodness_of_Fit.R` |
| `powerlaw_bootstrap.rda` | Power-law bootstrap p-values, created by `04_Powerlaw_Bootstrap.R` |

`sorted_mean_incomes_99.feather` contains 99 rows, one for each parameter set, with the following variables:

- `ω`: Proportion of wealth lost by the loser (0.1--0.9)
- `α`: Mixing parameter between competition and investment (0.0--1.0)
- `mean_incomes`: Wealth of the N = 1000 agents at T = 100, sorted in ascending order and averaged by rank over 100 trials (JSON string)

`vuong_pvalue.rda` and `powerlaw_bootstrap.rda` contain the same data (`mean_income`) and the test results:

- `pl_vs_lnorm_stat`: One-sided p-value of the Vuong test, power-law vs. log-normal (Figure 5)
- `lnorm_vs_gamma_stat`: One-sided p-value of the Vuong test, gamma vs. log-normal (Figure 6)
- `p_value`: Bootstrap p-value of the power-law goodness-of-fit test

### Figures (`results/`)

| Paper | File | Description | Created by |
| --- | --- | --- | --- |
| Figure 1 | `Fig1_gini_convergence.pdf` | Gini coefficient over time for each (ω, α), with α = 1.0 highlighted | `01_run_simulation.ipynb` |
| Figure 2(a) | `Fig2a_mean.pdf` | Mean wealth (T = 10) | `01_run_simulation.ipynb` |
| Figure 2(b) | `Fig2b_gini.pdf` | Gini coefficient (T = 100) | `01_run_simulation.ipynb` |
| Figure 3(a) | `Fig3a_ccdf_gamma.pdf` | CCDF and gamma fits (α = 0.0) | `02_Fitting.R` |
| Figure 3(b) | `Fig3b_ccdf_lnorm.pdf` | CCDF and log-normal fits (α = 1.0) | `02_Fitting.R` |
| Figure 4(a) | `Fig4a_ccdf_alpha0.2.pdf` | CCDF with gamma and log-normal fits (α = 0.2) | `02_Fitting.R` |
| Figure 4(b) | `Fig4b_ccdf_alpha0.4.pdf` | CCDF with gamma and log-normal fits (α = 0.4) | `02_Fitting.R` |
| Figure 4(c) | `Fig4c_ccdf_alpha0.6.pdf` | CCDF with gamma and log-normal fits (α = 0.6) | `02_Fitting.R` |
| Figure 4(d) | `Fig4d_ccdf_alpha0.8.pdf` | CCDF with gamma and log-normal fits (α = 0.8) | `02_Fitting.R` |
| Figure 5 | `Fig5_vuong_pl_vs_lnorm.pdf` | Vuong test p-values, power-law vs. log-normal for the upper tail | `05_Plot_Heatmaps.R` |
| Figure 6 | `Fig6_vuong_lnorm_vs_gamma.pdf` | Vuong test p-values, gamma vs. log-normal for the bulk | `05_Plot_Heatmaps.R` |
| Section 4 (text) | `powerlaw_pvalue.pdf` | Bootstrap p-values of the power-law goodness-of-fit test for the upper tail | `05_Plot_Heatmaps.R` |

## Requirements

- Julia 1.10.6 (packages pinned in `julia_simulation/Project.toml` and `julia_simulation/Manifest.toml`)
- R 4.4.2 (packages pinned with [renv](https://rstudio.github.io/renv/) in `renv.lock`), including
  - poweRlaw 1.0.0
  - fitdistrplus 1.2-2
  - ggplot2 4.0.0
  - arrow 22.0.0

Tested on macOS (Apple silicon).

## Usage

### Simulation (Julia)

```sh
julia +1.10.6 --project=julia_simulation -e 'using Pkg; Pkg.instantiate()'
```

Run `julia_simulation/code/01_run_simulation.ipynb` with this environment
(a Julia 1.10.6 kernel from [IJulia](https://github.com/JuliaLang/IJulia.jl) is required).

### Analysis (R)

Start R 4.4.2 in the repository root (renv is activated by `.Rprofile`), restore the packages, and run the scripts in numerical order:

```r
renv::restore()

source("r_analysis/script/02_Fitting.R")
source("r_analysis/script/03_Goodness_of_Fit.R")
source("r_analysis/script/04_Powerlaw_Bootstrap.R")
source("r_analysis/script/05_Plot_Heatmaps.R")
```

Figures are saved to `results/` and the test results to `data/`.

## License

MIT
