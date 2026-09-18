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
