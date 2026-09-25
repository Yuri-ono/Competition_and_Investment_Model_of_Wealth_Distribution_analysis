library(here)
i_am("r_analysis/script/02_Fitting.R")
library(jsonlite)
library(poweRlaw)
library(fitdistrplus)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(arrow)
library(tibble)
library(patchwork)
library(purrr)

# Plotting helpers
plot_ccdf <- function(data) {
  ggplot(data.frame(x = data), aes(x = x, y = 1 - ecdf(x)(x))) +
    geom_point() +
    scale_y_continuous(trans = "log10", breaks = 10^(0:-10)) +
    labs(x = "Income", y = "CCDF (1 - CDF)") +
    theme_minimal()
}
plot_ccdf_log_log <- function(data, title) {
  ggplot(data.frame(x = data), aes(x = x, y = 1 - ecdf(x)(x))) +
    geom_point() +
    scale_x_continuous(trans = "log10", breaks = 10^(0:10))+
    scale_y_continuous(trans = "log10", breaks = 10^(0:-10))+
    labs(x = "Income", y = "CCDF (1 - CDF)") +
    theme_minimal()
}

######################
####----Data----####
######################
data_path <- here("data", "sorted_mean_incomes_99.feather")

df <- read_feather(data_path)
df$mean_incomes <- lapply(df$mean_incomes, fromJSON)
df$mean_incomes
df <- df %>%
  rename(mean_income = mean_incomes)
save(df,file = here("data", "sorted_mean_incomes_99.rda"))

load(file = here("data", "sorted_mean_incomes_99.rda"))
mergedf <- df
df <- mergedf %>% filter(α == 0.0)
n_series <- length(df$mean_income)

all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  
  gamma_fit <- fitdist(data, "gamma", method = "mle")
  
  sorted_data <- sort(data)
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  gamma_ccdf <- 1 - pgamma(sorted_data,
                           shape = gamma_fit$estimate["shape"],
                           rate  = gamma_fit$estimate["rate"])
  
  data.frame(
    Income = rep(sorted_data[1:length(sorted_data)], 2),
    CCDF = c(emp_ccdf[1:length(sorted_data)], gamma_ccdf[1:length(sorted_data)]),
    Type = rep(c("Simulated Data", "Gamma Fit"), each = length(sorted_data[1:length(sorted_data)])),
    Index = df$ω[i]
  )
})

label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # plotmath string for parse = TRUE

ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Index, ncol = 3, scales = "free_x") +
  labs(x = "Income",
       y = "CCDF") +
  scale_color_manual(values = c(
    "Simulated Data" = "#1f77b4",
    "Gamma Fit" = "#2ca02c"
  )) +
  scale_linetype_manual(values = c(
    "Simulated Data" = "solid",
    "Gamma Fit" = "dashed"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,
    label.size = 0.7,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  )

ggsave(
  here("results", "Fig3a_ccdf_gamma.pdf"),
  width = 12,
  height = 10,
  units = "in"
)

###1.0##################
load(file = here("data", "sorted_mean_incomes_99.rda"))
mergedf <- df
df <- mergedf %>% filter(α == 1.0)

n_series <- length(df$mean_income)

all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  
  lnorm_fit <- fitdist(data, "lnorm", method = "mle")
  
  sorted_data <- sort(data)
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  lnorm_ccdf <- 1 - plnorm(sorted_data,
                           meanlog = lnorm_fit$estimate["meanlog"],
                           sdlog = lnorm_fit$estimate["sdlog"])
  
  data.frame(
    Income = rep(sorted_data[1:length(sorted_data)], 2),
    CCDF = c(emp_ccdf[1:length(sorted_data)], lnorm_ccdf[1:length(sorted_data)]),
    Type = rep(c("Simulated Data", "Log-Normal Fit"), each = length(sorted_data[1:length(sorted_data)])),
    Index = df$ω[i]
  )
})

label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # plotmath string for parse = TRUE
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Index, ncol = 3, scales = "free_x") +
  labs(x = "Income",
       y = "CCDF") +
  scale_color_manual(values = c(
    "Simulated Data" = "#1f77b4",
    "Log-Normal Fit" = "#d62728"
  )) +
  scale_linetype_manual(values = c(
    "Simulated Data" = "solid",
    "Log-Normal Fit" = "dotdash"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,
    label.size = 0.7,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  )

ggsave(
  here("results", "Fig3b_ccdf_lnorm.pdf"),
  width = 12,
  height = 10,
  units = "in"
)

###0.2----------------------------

load(file = here("data", "sorted_mean_incomes_99.rda"))
mergedf <- df
df <- mergedf %>% filter(α == 0.2)

n_series <- length(df$mean_income)

all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  
  sorted_data <- sort(data)
  
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  
  gamma_fit <- fitdist(data, "gamma", method = "mme")
  gamma_ccdf <- 1 - pgamma(sorted_data,
                           shape = gamma_fit$estimate["shape"],
                           rate  = gamma_fit$estimate["rate"])
  
  lnorm_fit <- fitdist(data, "lnorm", method = "mme")
  lnorm_ccdf <- 1 - plnorm(sorted_data,
                           meanlog = lnorm_fit$estimate["meanlog"],
                           sdlog = lnorm_fit$estimate["sdlog"])
  
  data.frame(
    Income = rep(sorted_data, 3),
    CCDF = c(emp_ccdf, gamma_ccdf, lnorm_ccdf),
    Type = rep(c("Simulated Data", "Gamma Fit", "Log-Normal Fit"), each = length(sorted_data)),
    Index = df$ω[i]
  )
})

label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # plotmath string for parse = TRUE
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Index, ncol = 3, scales = "free_x") +
  labs(x = "Income",
       y = "CCDF") +
  scale_color_manual(values = c(
    "Simulated Data" = "#1f77b4",
    "Gamma Fit" = "#2ca02c",
    "Log-Normal Fit" = "#d62728"
  )) +
  scale_linetype_manual(values = c(
    "Simulated Data" = "solid",
    "Gamma Fit" = "dashed",
    "Log-Normal Fit" = "dotdash"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,
    label.size = 0.7,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  ) -> a2

###0.4----------------------------

df <- mergedf %>% filter(α == 0.4)

n_series <- length(df$mean_income)

all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  
  sorted_data <- sort(data)
  
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  
  gamma_fit <- fitdist(data, "gamma", method = "mme")
  gamma_ccdf <- 1 - pgamma(sorted_data,
                           shape = gamma_fit$estimate["shape"],
                           rate  = gamma_fit$estimate["rate"])
  
  lnorm_fit <- fitdist(data, "lnorm", method = "mme")
  lnorm_ccdf <- 1 - plnorm(sorted_data,
                           meanlog = lnorm_fit$estimate["meanlog"],
                           sdlog = lnorm_fit$estimate["sdlog"])
  
  data.frame(
    Income = rep(sorted_data, 3),
    CCDF = c(emp_ccdf, gamma_ccdf, lnorm_ccdf),
    Type = rep(c("Simulated Data", "Gamma Fit", "Log-Normal Fit"), each = length(sorted_data)),
    Index = df$ω[i]
  )
})

label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # plotmath string for parse = TRUE
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Index, ncol = 3, scales = "free_x") +
  labs(x = "Income",
       y = "CCDF") +
  scale_color_manual(values = c(
    "Simulated Data" = "#1f77b4",
    "Gamma Fit" = "#2ca02c",
    "Log-Normal Fit" = "#d62728"
  )) +
  scale_linetype_manual(values = c(
    "Simulated Data" = "solid",
    "Gamma Fit" = "dashed",
    "Log-Normal Fit" = "dotdash"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,
    label.size = 0.7,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  ) -> a4

###0.6----------------------------

df <- mergedf %>% filter(α == 0.6)

n_series <- length(df$mean_income)

all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  
  sorted_data <- sort(data)
  
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  
  gamma_fit <- fitdist(data, "gamma", method = "mme")
  gamma_ccdf <- 1 - pgamma(sorted_data,
                           shape = gamma_fit$estimate["shape"],
                           rate  = gamma_fit$estimate["rate"])
  
  lnorm_fit <- fitdist(data, "lnorm", method = "mme")
  lnorm_ccdf <- 1 - plnorm(sorted_data,
                           meanlog = lnorm_fit$estimate["meanlog"],
                           sdlog = lnorm_fit$estimate["sdlog"])
  
  data.frame(
    Income = rep(sorted_data, 3),
    CCDF = c(emp_ccdf, gamma_ccdf, lnorm_ccdf),
    Type = rep(c("Simulated Data", "Gamma Fit", "Log-Normal Fit"), each = length(sorted_data)),
    Index = df$ω[i]
  )
})

label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # plotmath string for parse = TRUE
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Index, ncol = 3, scales = "free_x") +
  labs(x = "Income",
       y = "CCDF") +
  scale_color_manual(values = c(
    "Simulated Data" = "#1f77b4",
    "Gamma Fit" = "#2ca02c",
    "Log-Normal Fit" = "#d62728"
  )) +
  scale_linetype_manual(values = c(
    "Simulated Data" = "solid",
    "Gamma Fit" = "dashed",
    "Log-Normal Fit" = "dotdash"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,
    label.size = 0.7,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  ) -> a6

###0.8----------------------------

df <- mergedf %>% filter(α == 0.8)

n_series <- length(df$mean_income)

all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  max_digits <- nchar(as.character(floor(max(data, na.rm = TRUE))))
  if (max_digits > 8) {
    scale_factor <- 10^(max_digits - 8)
    data <- data / scale_factor
  }
  
  sorted_data <- sort(data)
  
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  
  gamma_fit <- fitdist(data, "gamma", method = "mme")
  gamma_ccdf <- 1 - pgamma(sorted_data,
                           shape = gamma_fit$estimate["shape"],
                           rate  = gamma_fit$estimate["rate"])
  
  lnorm_fit <- fitdist(data, "lnorm", method = "mme")
  lnorm_ccdf <- 1 - plnorm(sorted_data,
                           meanlog = lnorm_fit$estimate["meanlog"],
                           sdlog = lnorm_fit$estimate["sdlog"])
  
  data.frame(
    Income = rep(sorted_data, 3),
    CCDF = c(emp_ccdf, gamma_ccdf, lnorm_ccdf),
    Type = rep(c("Simulated Data", "Gamma Fit", "Log-Normal Fit"), each = length(sorted_data)),
    Index = df$ω[i]
  )
})

label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # plotmath string for parse = TRUE
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Index, ncol = 3, scales = "free_x") +
  labs(x = "Income",
       y = "CCDF") +
  scale_color_manual(values = c(
    "Simulated Data" = "#1f77b4",
    "Gamma Fit" = "#2ca02c",
    "Log-Normal Fit" = "#d62728"
  )) +
  scale_linetype_manual(values = c(
    "Simulated Data" = "solid",
    "Gamma Fit" = "dashed",
    "Log-Normal Fit" = "dotdash"
  )) +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,
    label.size = 0.7,
    label.r = unit(0.15, "lines"),
    fill = "white",
    color = "black"
  ) -> a8

a2
ggsave(
  here("results", "Fig4a_ccdf_alpha0.2.pdf"),
  plot = a2,
  width = 12,
  height = 10,
  units = "in"
)
a4
ggsave(
  here("results", "Fig4b_ccdf_alpha0.4.pdf"),
  plot = a4,
  width = 12,
  height = 10,
  units = "in"
)
a6
ggsave(
  here("results", "Fig4c_ccdf_alpha0.6.pdf"),
  plot = a6,
  width = 12,
  height = 10,
  units = "in"
)
a8
ggsave(
  here("results", "Fig4d_ccdf_alpha0.8.pdf"),
  plot = a8,
  width = 12,
  height = 10,
  units = "in"
)
