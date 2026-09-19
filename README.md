# Competition_and_Investment_Model_of_Wealth_Distribution_analysis

This simulations were performed using Julia 1.10.6, and subsequent data analysis were conducted using R 4.4.2.

## Environment

### Julia (1.10.6)

Package versions are pinned in `julia_simulation/Project.toml` and `julia_simulation/Manifest.toml`.

```sh
julia +1.10.6 --project=julia_simulation -e 'using Pkg; Pkg.instantiate()'
```

Run `julia_simulation/code/01_run_simulation.ipynb` with this environment activated
(a Julia 1.10.6 kernel from [IJulia](https://github.com/JuliaLang/IJulia.jl) is required).

### R (4.4.2)

Package versions are pinned with [renv](https://rstudio.github.io/renv/) in `r_analysis/renv.lock`.
Start R 4.4.2 in `r_analysis/` (renv is activated by `r_analysis/.Rprofile`) and restore the packages:

```r
renv::restore()
```

Then run the scripts in `r_analysis/script/` in numerical order.

## Scripts and outputs

| Script | Output | Paper |
| --- | --- | --- |
| `julia_simulation/code/01_run_simulation.ipynb` | `data/sorted_mean_incomes_99.feather`, `results/Fig1_gini_convergence.pdf`, `results/Fig2a_mean.pdf`, `results/Fig2b_gini.pdf` | Figures 1 and 2 |
| `r_analysis/script/02_Fitting.R` | `mergedf_01.rda`, `results/CCDF_gamma.pdf`, `results/CCDF_lnorm.pdf` | Figure 3 |
| | `results/CCDF_gl_2.pdf`, `results/CCDF_gl_4.pdf`, `results/CCDF_gl_6.pdf`, `results/CCDF_gl_8.pdf` | Figure 4 |
| `r_analysis/script/03_Goodness_of_Fit.R` | `results/merge_01_vuong_pvalue_rerun.rda` | Recomputes the Vuong tests (about 3 minutes) |
| `r_analysis/script/04_Powerlaw_Bootstrap.R` | `results/powerlaw_pvalue_rerun.rda` | Recomputes the power-law bootstrap p-values (about 1 hour) |
| `r_analysis/script/05_Plot_Heatmaps.R` | `results/Fig5_vuong_pl_vs_lnorm.pdf`, `results/Fig6_vuong_lnorm_vs_gamma.pdf`, `results/powerlaw_pvalue.pdf` | Figures 5 and 6, power-law p-values in the text |

The values reported in the paper for the goodness-of-fit tests are stored in
`results/merge_01_vuong_pvalue_main.rda` (Figures 5 and 6) and
`results/mergedata_01_main_pvalue_toward_powerlaw.rda` (power-law p-values),
and `05_Plot_Heatmaps.R` draws the figures from these files.
The gamma vs. log-normal comparison (Figure 6) and the bootstrap p-values are stochastic
and the stored results were computed without a fixed seed,
so the values recomputed by `03` and `04` differ slightly from the paper.

`data/sorted_mean_incomes_99.feather` is the simulation output used in the paper.
Re-running the simulation reproduces its statistical properties, but not the exact values.
