# Competition and Investment Model of Wealth Distribution

Replication code and data for:

> [AUTHORS] ([YEAR]). Competition and investment model of wealth distribution. *The Journal of Mathematical Sociology*. [DOI]

## Files

### Julia code (`julia_simulation/`)

| File | Description | Figures |
| --- | --- | --- |
| `code/01_run_simulation.ipynb` | Simulation of the unified model: mean and Gini coefficient, time evolution of the Gini coefficient, and generation of the wealth distributions used in the analysis | Figures 1, 2 |

### R code (`r_analysis/script/`)

| File | Description | Figures |
| --- | --- | --- |
| `02_Fitting.R` | CCDF plots with gamma and log-normal fits | Figures 3, 4 |
| `03_Goodness_of_Fit.R` | Vuong tests: power-law vs. log-normal for the upper tail, gamma vs. log-normal for the bulk (recomputation, about 3 minutes) | Figures 5, 6 |
| `04_Powerlaw_Bootstrap.R` | Power-law goodness-of-fit test with bootstrap p-values (recomputation, about 1 hour) | Section 4 (text) |
| `05_Plot_Heatmaps.R` | Heatmaps of the test results | Figures 5, 6 |

### Data (`data/`, `results/`)

| File | Description |
| --- | --- |
| `data/sorted_mean_incomes_99.feather` | Simulated wealth distributions for each (ω, α) |
| `results/merge_01_vuong_pvalue_main.rda` | Vuong test results reported in the paper (Figures 5 and 6) |
| `results/mergedata_01_main_pvalue_toward_powerlaw.rda` | Power-law bootstrap p-values reported in the paper |

`data/sorted_mean_incomes_99.feather` contains 99 rows, one for each parameter set, with the following variables:

- `ω`: Proportion of wealth lost by the loser (0.1--0.9)
- `α`: Mixing parameter between competition and investment (0.0--1.0)
- `mean_incomes`: Wealth of the N = 1000 agents at T = 100, sorted in ascending order and averaged by rank over 100 trials (JSON string)

The `.rda` files in `results/` contain the same data (`mean_income`) and the test results:

- `pl_vs_lnorm_stat`: One-sided p-value of the Vuong test, power-law vs. log-normal (Figure 5)
- `lnorm_vs_gamma_stat`: One-sided p-value of the Vuong test, gamma vs. log-normal (Figure 6)
- `p_value`: Bootstrap p-value of the power-law goodness-of-fit test

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

- `05_Plot_Heatmaps.R` draws Figures 5 and 6 from the stored results in `results/`.
  If these files are removed, it uses the results recomputed by `03` and `04` instead.
- The gamma vs. log-normal comparison (Figure 6) and the bootstrap p-values are stochastic,
  and the stored results were computed without a fixed seed,
  so the values recomputed by `03` and `04` differ slightly from the paper.
- `data/sorted_mean_incomes_99.feather` is the simulation output used in the paper.
  Re-running the simulation reproduces its statistical properties, but not the exact values.

## License

MIT
