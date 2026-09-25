library(here)
i_am("r_analysis/script/03_Goodness_of_Fit.R")
library(fitdistrplus)
library(poweRlaw)

vuong <- function(censdata,ln_fit,gamma_fit){
  
  shape <- gamma_fit$estimate["shape"]
  rate  <- gamma_fit$estimate["rate"]
  
  # Likelihood of each observation under the gamma fit
  likelihoods <- mapply(function(l, r) {
    if (is.na(r)) {
      # Right-censored: S(x) = 1 - F(x)
      return(1 - pgamma(l, shape = shape, rate = rate))
    } else if (l == r) {
      # Uncensored: f(x)
      return(dgamma(l, shape = shape, rate = rate))
    } else {
      # Interval-censored (not used)
      return(pgamma(r, shape = shape, rate = rate) - pgamma(l, shape = shape, rate = rate))
    }
  }, censdata$left, censdata$right)
  
  # Log-likelihoods
  galoglikelihoods <- log(pmax(likelihoods, .Machine$double.eps))
  
  # Log-normal distribution
  #ln_fit <- fitdistcens(censdata, "lnorm",optim.method = "SANN")
  meanlog <- ln_fit$estimate["meanlog"]
  sdlog   <- ln_fit$estimate["sdlog"]
  
  # Likelihood of each observation under the log-normal fit
  likelihoods <- mapply(function(l, r) {
    if (is.na(r)) {
      # Right-censored: S(x) = 1 - F(x)
      return(1 - plnorm(l, meanlog = meanlog, sdlog = sdlog))
    } else if (l == r) {
      # Uncensored: f(x)
      return(dlnorm(l, meanlog = meanlog, sdlog = sdlog))
    } else {
      # Interval-censored (not used)
      return(plnorm(r, meanlog = meanlog, sdlog = sdlog) -
               plnorm(l, meanlog = meanlog, sdlog = sdlog))
    }
  }, censdata$left, censdata$right)
  
  # Log-likelihoods
  lnloglikelihoods <- log(pmax(likelihoods, .Machine$double.eps))
  
  # Vuong test statistic
  LR <- lnloglikelihoods - galoglikelihoods # pointwise log-likelihood ratios
  
  omega_est <- sd(LR) # as in poweRlaw::compare_distributions()
  
  n <- sqrt(length(likelihoods)) # square root of the sample size
  return (n*mean(LR)/omega_est) # normalized log-likelihood ratio (Vuong, 1989)
  # positive values favor the log-normal distribution
}

load(file = here("data", "sorted_mean_incomes_99.rda"))

set.seed(1)

# Vectors to store the results
test_stats <- numeric(nrow(df))
vuongst <- numeric(nrow(df))

for (i in 1:nrow(df)) {
  data <- df$mean_income[[i]]
  
  # Rescale large values to avoid overflow
  max_digits <- nchar(as.character(floor(max(data, na.rm = TRUE))))
  if (max_digits > 13) {
    scale_factor <- 10^(max_digits - 13)
    data <- data / scale_factor
  }
  
  # Fit power-law and log-normal models to the upper tail
  pl_model <- conpl$new(data)
  est <- estimate_xmin(pl_model, xmax = max(data))
  pl_model$setXmin(est)
  
  pl_model_ln <- conlnorm$new(data)
  pl_model_ln$setXmin(pl_model$getXmin())
  pl_model_ln$setPars(estimate_pars(pl_model_ln))
  
  # Vuong test: power-law vs. log-normal
  compd <- compare_distributions(pl_model, pl_model_ln)
  #test_stats[i] <- compd$test_statistic
  test_stats[i] <- compd$p_one_sided
  
  # ---- Bulk of the distribution ----
  xmin <- pl_model$xmin  # lower bound of the power-law tail
  censdata <- data.frame(
    left  = ifelse(data < xmin, data, xmin),  # observations x >= xmin are right-censored at xmin
    right = ifelse(data < xmin, data, NA)     # so they contribute only P(X >= xmin)
  )

  gamma_fit <- fitdistcens(censdata, "gamma", optim.method = "SANN")
  ln_fit <- fitdistcens(censdata, "lnorm", optim.method = "SANN")

  vuongst[i] <- 1 - pnorm(vuong(censdata, ln_fit, gamma_fit))
}
# Add the results as columns of df
df$pl_vs_lnorm_stat <- test_stats
df$lnorm_vs_gamma_stat <- vuongst

# These are the values shown in Figures 5 and 6, drawn by 05_Plot_Heatmaps.R.
save(df, file = here("data", "vuong_pvalue.rda"))

df$pl_vs_lnorm_stat |> summary()
df$lnorm_vs_gamma_stat |> summary()
