library(here)
i_am("r_analysis/script/03_Goodness_of_Fit.R")
library(fitdistrplus)
library(jsonlite)
library(dplyr)
library(poweRlaw)
library(gridExtra)
library(ggplot2)
library(tidyr)

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

load(file = here("mergedf_01.rda"))

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
  #censdata <- data.frame(
  #  left  = ifelse(data < xmin, data, xmin),  
  #  right = ifelse(data < xmin, data, NA)
  #)
  censdata <- data.frame(
    left  = data,   # all observations
    right = ifelse(data >= xmin, NA, data)  # observations x >= xmin are right-censored (NA)
  )
  
  gamma_fit <- fitdistcens(censdata, "gamma", optim.method = "SANN")
  ln_fit <- fitdistcens(censdata, "lnorm", optim.method = "SANN")

  vuongst[i] <- 1 - pnorm(vuong(censdata, ln_fit, gamma_fit))
}
# Add the results as columns of df
df$pl_vs_lnorm_stat <- test_stats
df$lnorm_vs_gamma_stat <- vuongst

# The values reported in the paper are stored in
# results/merge_01_vuong_pvalue_main.rda. SANN is stochastic and the stored
# results were computed without a fixed seed, so lnorm_vs_gamma_stat differs
# slightly from the paper; the output is saved under a different name.
save(df, file = here("results", "merge_01_vuong_pvalue_rerun.rda"))



# Heatmap
ggplot(df, aes(x = factor(ω), y = factor(α), fill = pl_vs_lnorm_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5, name = "Test Statistic") +
  labs(
    title = "Power-law vs Lognormal p_one_sided by (ω, α)",
    x = expression(omega),
    y = expression(alpha)
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
#ggsave("merge_vuong_pvalue_high_main.pdf")
df$pl_vs_lnorm_stat |> summary()


# Heatmap
ggplot(df, aes(x = factor(ω), y = factor(α), fill = lnorm_vs_gamma_stat)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5, name = "Test Statistic") +
  labs(
    title = "Lognormal vs Gamma p_one_sided by (ω, α)",
    x = expression(omega),
    y = expression(alpha)
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
#ggsave("merge_vuong_pvalue_low_main.pdf")
df$lnorm_vs_gamma_stat |> summary()
# Red favors gamma, blue favors log-normal