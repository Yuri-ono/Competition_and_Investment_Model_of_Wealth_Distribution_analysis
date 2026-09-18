#load data
library(here)
i_am("r_analysis/script/02_Fitting.R")
library(jsonlite)
# fit data
library(poweRlaw)
library(fitdistrplus)
# graphics
library(dplyr)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(arrow)
library(tibble)
library(patchwork)
library(purrr)  # for map and related functions

# プロット作成用の関数
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
####----データ----####
######################
data_path <- here("data", "sorted_mean_incomes_99.feather")

df <- read_feather(data_path)
df$mean_incomes <- lapply(df$mean_incomes, fromJSON)  # JSON文字列をリストに変換
df$mean_incomes
df <- df %>%
  rename(mean_income = mean_incomes)
save(df,file = here("mergedf_01.rda"))

load(file = here("mergedf_01.rda"))
mergedf <- df
df <- mergedf %>% filter(α == 0.0)
# データ系列の数（たとえば9）
n_series <- length(df$mean_income)

# データ作成
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
    Index = df$ω[i]  # ここでは数値だけ保持
  )
})

# 図に書き込む用のラベルデータ
label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,  # 左下あたりに配置する（少し余裕を持たせる）
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # expression形式にするために文字列作成

ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +  # 線も少し太め
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
  theme_minimal(base_size = 12) +  # 基本フォントサイズ
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),  # 全体の基本文字サイズ
    axis.title = element_text(size = 16, face = "bold"),  # 軸タイトルは大きく太く
    axis.text = element_text(size = 12),  # 軸目盛り
    legend.text = element_text(size = 12),  # 凡例文字
    legend.title = element_blank()  # 凡例タイトルはいらない or 小さく
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,         # ★ラベルの文字を大きめに
    label.size = 0.7, # ★ボックス枠線もややしっかり
    label.r = unit(0.15, "lines"),
    fill = "white",   # ★背景白
    color = "black"   # ★枠線黒
  )

ggsave(
  here("results", "CCDF_gamma.pdf"),
  width = 12,
  height = 10,
  units = "in"
)

###1.0##################
# データ読み込みとフィルター
load(file = here("mergedf_01.rda"))
mergedf <- df
df <- mergedf %>% filter(α == 1.0)

# データ系列の数（たとえば9）
n_series <- length(df$mean_income)

# 各系列について処理
all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  
  # フィット
  lnorm_fit <- fitdist(data, "lnorm", method = "mle")
  
  sorted_data <- sort(data)
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  lnorm_ccdf <- 1 - plnorm(sorted_data,
                           meanlog = lnorm_fit$estimate["meanlog"],
                           sdlog = lnorm_fit$estimate["sdlog"])
  
  # データフレームにまとめて返す
  data.frame(
    Income = rep(sorted_data[1:length(sorted_data)], 2),
    CCDF = c(emp_ccdf[1:length(sorted_data)], lnorm_ccdf[1:length(sorted_data)]),
    Type = rep(c("Simulated Data", "Log-Normal Fit"), each = length(sorted_data[1:length(sorted_data)])),
    Index = df$ω[i]  # ファセットタイトル用
  )
})

# 図に書き込む用のラベルデータ
label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,  # 左下あたりに配置する（少し余裕を持たせる）
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # expression形式にするために文字列作成
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +  # 線も少し太め
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
  theme_minimal(base_size = 12) +  # 基本フォントサイズ
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),  # 全体の基本文字サイズ
    axis.title = element_text(size = 16, face = "bold"),  # 軸タイトルは大きく太く
    axis.text = element_text(size = 12),  # 軸目盛り
    legend.text = element_text(size = 12),  # 凡例文字
    legend.title = element_blank()  # 凡例タイトルはいらない or 小さく
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,         # ★ラベルの文字を大きめに
    label.size = 0.7, # ★ボックス枠線もややしっかり
    label.r = unit(0.15, "lines"),
    fill = "white",   # ★背景白
    color = "black"   # ★枠線黒
  )

ggsave(
  here("results", "CCDF_lnorm.pdf"),
  width = 12,
  height = 10,
  units = "in"
)

###0.2----------------------------

load(file = here("mergedf_01.rda"))
mergedf <- df
df <- mergedf %>% filter(α == 0.2)

n_series <- length(df$mean_income)

all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  
  sorted_data <- sort(data)
  
  # 経験的CCDF
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  
  # フィッティング
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

# 図に書き込む用のラベルデータ
label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,  # 左下あたりに配置する（少し余裕を持たせる）
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # expression形式にするために文字列作成
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +  # 線も少し太め
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
  theme_minimal(base_size = 12) +  # 基本フォントサイズ
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),  # 全体の基本文字サイズ
    axis.title = element_text(size = 16, face = "bold"),  # 軸タイトルは大きく太く
    axis.text = element_text(size = 12),  # 軸目盛り
    legend.text = element_text(size = 12),  # 凡例文字
    legend.title = element_blank()  # 凡例タイトルはいらない or 小さく
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,         # ★ラベルの文字を大きめに
    label.size = 0.7, # ★ボックス枠線もややしっかり
    label.r = unit(0.15, "lines"),
    fill = "white",   # ★背景白
    color = "black"   # ★枠線黒
  ) -> a2

###0.4----------------------------

df <- mergedf %>% filter(α == 0.4)

n_series <- length(df$mean_income)

all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  
  sorted_data <- sort(data)
  
  # 経験的CCDF
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  
  # フィッティング
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

# 図に書き込む用のラベルデータ
label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,  # 左下あたりに配置する（少し余裕を持たせる）
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # expression形式にするために文字列作成
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +  # 線も少し太め
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
  theme_minimal(base_size = 12) +  # 基本フォントサイズ
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),  # 全体の基本文字サイズ
    axis.title = element_text(size = 16, face = "bold"),  # 軸タイトルは大きく太く
    axis.text = element_text(size = 12),  # 軸目盛り
    legend.text = element_text(size = 12),  # 凡例文字
    legend.title = element_blank()  # 凡例タイトルはいらない or 小さく
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,         # ★ラベルの文字を大きめに
    label.size = 0.7, # ★ボックス枠線もややしっかり
    label.r = unit(0.15, "lines"),
    fill = "white",   # ★背景白
    color = "black"   # ★枠線黒
  ) -> a4

###0.6----------------------------

df <- mergedf %>% filter(α == 0.6)

n_series <- length(df$mean_income)

all_data <- map_dfr(1:n_series, function(i) {
  data <- df$mean_income[[i]]
  
  sorted_data <- sort(data)
  
  # 経験的CCDF
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  
  # フィッティング
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

# 図に書き込む用のラベルデータ
label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,  # 左下あたりに配置する（少し余裕を持たせる）
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # expression形式にするために文字列作成
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +  # 線も少し太め
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
  theme_minimal(base_size = 12) +  # 基本フォントサイズ
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),  # 全体の基本文字サイズ
    axis.title = element_text(size = 16, face = "bold"),  # 軸タイトルは大きく太く
    axis.text = element_text(size = 12),  # 軸目盛り
    legend.text = element_text(size = 12),  # 凡例文字
    legend.title = element_blank()  # 凡例タイトルはいらない or 小さく
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,         # ★ラベルの文字を大きめに
    label.size = 0.7, # ★ボックス枠線もややしっかり
    label.r = unit(0.15, "lines"),
    fill = "white",   # ★背景白
    color = "black"   # ★枠線黒
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
  
  # 経験的CCDF
  emp_ccdf <- 1 - ecdf(data)(sorted_data)
  
  # フィッティング
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

# 図に書き込む用のラベルデータ
label_data <- all_data %>%
  group_by(Index) %>%
  summarise(
    Income = min(Income) * 1.2,  # 左下あたりに配置する（少し余裕を持たせる）
    CCDF = min(CCDF) * 1.2
  ) %>%
  mutate(Label = paste0("omega == ", Index))  # expression形式にするために文字列作成
ggplot(all_data, aes(x = Income, y = CCDF, color = Type, linetype = Type)) +
  geom_line(size = 1.2) +  # 線も少し太め
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
  theme_minimal(base_size = 12) +  # 基本フォントサイズ
  theme(
    aspect.ratio = 4/5,
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),  # 全体の基本文字サイズ
    axis.title = element_text(size = 16, face = "bold"),  # 軸タイトルは大きく太く
    axis.text = element_text(size = 12),  # 軸目盛り
    legend.text = element_text(size = 12),  # 凡例文字
    legend.title = element_blank()  # 凡例タイトルはいらない or 小さく
  ) +
  geom_label(
    data = label_data,
    aes(x = Income, y = CCDF, label = Label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -3,
    size = 6,         # ★ラベルの文字を大きめに
    label.size = 0.7, # ★ボックス枠線もややしっかり
    label.r = unit(0.15, "lines"),
    fill = "white",   # ★背景白
    color = "black"   # ★枠線黒
  ) -> a8

a2
ggsave(
  here("results", "CCDF_gl_2.pdf"),
  width = 12,
  height = 10,
  units = "in"
)
a4
ggsave(
  here("results", "CCDF_gl_4.pdf"),
  width = 12,
  height = 10,
  units = "in"
)
a6
ggsave(
  here("results", "CCDF_gl_6.pdf"),
  width = 12,
  height = 10,
  units = "in"
)
a8
ggsave(
  here("results", "CCDF_gl_8.pdf"),
  width = 12,
  height = 10,
  units = "in"
)
