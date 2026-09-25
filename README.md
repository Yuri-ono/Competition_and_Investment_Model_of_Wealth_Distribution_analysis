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

### Data (`data/`, `results/`)

| File | Description |
| --- | --- |
| `data/sorted_mean_incomes_99.feather` | Simulated wealth distributions for each (ω, α) |
| `data/sorted_mean_incomes_99.rda` | The same data in R format (`mean_income`), created by `02_Fitting.R` and used by `03` and `04` |
| `results/vuong_pvalue.rda` | Vuong test results (Figures 5 and 6), created by `03_Goodness_of_Fit.R` |
| `results/powerlaw_bootstrap.rda` | Power-law bootstrap p-values, created by `04_Powerlaw_Bootstrap.R` |

`data/sorted_mean_incomes_99.feather` contains 99 rows, one for each parameter set, with the following variables:

- `ω`: Proportion of wealth lost by the loser (0.1--0.9)
- `α`: Mixing parameter between competition and investment (0.0--1.0)
- `mean_incomes`: Wealth of the N = 1000 agents at T = 100, sorted in ascending order and averaged by rank over 100 trials (JSON string)

The `.rda` files in `results/` contain the same data (`mean_income`) and the test results:

- `pl_vs_lnorm_stat`: One-sided p-value of the Vuong test, power-law vs. log-normal (Figure 5)
- `lnorm_vs_gamma_stat`: One-sided p-value of the Vuong test, gamma vs. log-normal (Figure 6)
- `p_value`: Bootstrap p-value of the power-law goodness-of-fit test

### Figures (`results/`)

| Paper | File | Created by |
| --- | --- | --- |
| Figure 1 | `Fig1_gini_convergence.pdf` | `01_run_simulation.ipynb` |
| Figure 2(a) | `Fig2a_mean.pdf` | `01_run_simulation.ipynb` |
| Figure 2(b) | `Fig2b_gini.pdf` | `01_run_simulation.ipynb` |
| Figure 3(a) | `CCDF_gamma.pdf` (α = 0.0) | `02_Fitting.R` |
| Figure 3(b) | `CCDF_lnorm.pdf` (α = 1.0) | `02_Fitting.R` |
| Figure 4(a) | `CCDF_gl_2.pdf` (α = 0.2) | `02_Fitting.R` |
| Figure 4(b) | `CCDF_gl_4.pdf` (α = 0.4) | `02_Fitting.R` |
| Figure 4(c) | `CCDF_gl_6.pdf` (α = 0.6) | `02_Fitting.R` |
| Figure 4(d) | `CCDF_gl_8.pdf` (α = 0.8) | `02_Fitting.R` |
| Figure 5 | `Fig5_vuong_pl_vs_lnorm.pdf` | `05_Plot_Heatmaps.R` |
| Figure 6 | `Fig6_vuong_lnorm_vs_gamma.pdf` | `05_Plot_Heatmaps.R` |
| Section 4 (text) | `powerlaw_pvalue.pdf` | `05_Plot_Heatmaps.R` |

- `powerlaw_pvalue.pdf` shows the bootstrap p-values of the power-law goodness-of-fit test,
  which are reported in the text of Section 4 but not shown as a figure in the paper.
- Figure 3 shows maximum likelihood fits, and Figure 4 shows moment-matching fits (`method = "mme"`).
- In Figure 4(d), wealth is divided by a power of ten before fitting and plotting when its maximum
  has more than eight digits (ω >= 0.4), so the x-axis of these panels is not on the original scale.
  This shifts the curves but does not change their shape.

## Requirements

- Julia 1.10.6 (packages pinned in `julia_simulation/Project.toml` and `julia_simulation/Manifest.toml`)
- R 4.4.2 (packages pinned with [renv](https://rstudio.github.io/renv/) in `renv.lock`), including
  - poweRlaw 1.0.0
  - fitdistrplus 1.2-2
  - ggplot2 4.0.0
  - arrow 22.0.0

## Usage

### Simulation (Julia)

```sh
julia +1.10.6 --project=julia_simulation -e 'using Pkg; Pkg.instantiate()'
```

Run `julia_simulation/code/01_run_simulation.ipynb` with this environment
(a Julia 1.10.6 kernel from [IJulia](https://github.com/JuliaLang/IJulia.jl) is required).
Running the "Data Generate" section overwrites `data/sorted_mean_incomes_99.feather`.

### Analysis (R)

Start R 4.4.2 in the repository root (renv is activated by `.Rprofile`), restore the packages, and run the scripts in numerical order:

```r
renv::restore()

source("r_analysis/script/02_Fitting.R")
source("r_analysis/script/03_Goodness_of_Fit.R")
source("r_analysis/script/04_Powerlaw_Bootstrap.R")
source("r_analysis/script/05_Plot_Heatmaps.R")
```

Figures are saved to `results/`.

## Notes on reproducibility

- In the bulk comparison (Figure 6), observations with x >= xmin are right-censored at xmin,
  so that they contribute only the information that the value is at least xmin.
- The fits in `03` (SANN) and the bootstrap in `04` are stochastic, but both scripts fix the seed,
  so they reproduce the values in the figures with the R and package versions pinned here.
- `data/sorted_mean_incomes_99.feather` is the simulation output used in the paper.
  Running the notebook with Julia 1.10.6 regenerates it exactly (verified on macOS, Apple silicon).

## License

MIT
