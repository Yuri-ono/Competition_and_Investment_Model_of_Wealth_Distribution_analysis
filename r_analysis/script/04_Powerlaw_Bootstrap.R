library(here)
i_am("r_analysis/script/04_Powerlaw_Bootstrap.R")
library(poweRlaw)

# Power-law goodness-of-fit test for the upper tail (Clauset et al., 2009).
# These are the p-values reported in Section 4 and shown by 05_Plot_Heatmaps.R.
# bootstrap_p() is stochastic, so the seed below fixes the results.
# Runtime: about 1 hour in total.

load(file = here("mergedf_01.rda"))  # created by 02_Fitting.R

set.seed(1)

p_values <- numeric(nrow(df))
for (i in 1:nrow(df)) {
  data <- df$mean_income[[i]]
  # Rescale to at most 4 integer digits to stabilize the estimation
  max_digits <- nchar(as.character(floor(max(data, na.rm = TRUE))))
  if (max_digits > 4) {
    scale_factor <- 10^(max_digits - 4)
    data <- data / scale_factor
  }
  # Estimate xmin and the exponent by minimizing the KS statistic
  pl_model <- conpl$new(data)
  est <- estimate_xmin(pl_model, xmax = max(data))
  pl_model$setXmin(est)

  # Bootstrap p-value (p < 0.1 rejects the power-law hypothesis)
  booty <- bootstrap_p(pl_model)
  p_values[i] <- booty$p
}

df$p_value <- p_values

save(df, file = here("results", "powerlaw_bootstrap.rda"))
