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
  rate  <- 1/gamma_fit$estimate["scale"]
  #rate  <- gamma_fit$estimate["rate"]
  
  # 3. 各点の尤度・対数尤度を計算
  likelihoods <- mapply(function(l, r) {
    if (is.na(r)) {
      # 右検閲：S(x) = 1 - F(x)
      return(1 - pgamma(l, shape = shape, rate = rate))
    } else if (l == r) {
      # 完全観測：f(x)
      return(dgamma(l, shape = shape, rate = rate))
    } else {
      # 区間検閲（必要に応じて）
      return(pgamma(r, shape = shape, rate = rate) - pgamma(l, shape = shape, rate = rate))
    }
  }, censdata$left, censdata$right)
  
  # 4. 対数尤度
  galoglikelihoods <- log(pmax(likelihoods, .Machine$double.eps))
  
  # 対数正規分布
  #ln_fit <- fitdistcens(censdata, "lnorm",optim.method = "SANN")
  meanlog <- ln_fit$estimate["meanlog"]
  sdlog   <- ln_fit$estimate["sdlog"]
  
  # --- Step 2: 尤度を計算（検閲のタイプ別） ---
  likelihoods <- mapply(function(l, r) {
    if (is.na(r)) {
      # 右検閲: S(x) = 1 - F(x)
      return(1 - plnorm(l, meanlog = meanlog, sdlog = sdlog))
    } else if (l == r) {
      # 完全観測: f(x)
      return(dlnorm(l, meanlog = meanlog, sdlog = sdlog))
    } else {
      # 区間検閲（今回は出てこないかも）
      return(plnorm(r, meanlog = meanlog, sdlog = sdlog) -
               plnorm(l, meanlog = meanlog, sdlog = sdlog))
    }
  }, censdata$left, censdata$right)
  
  # --- Step 3: 対数尤度に変換 ---
  lnloglikelihoods <- log(pmax(likelihoods, .Machine$double.eps))
  
  ##-----あわせて計算
  LR <- lnloglikelihoods - galoglikelihoods #対数尤度比
  
  omega_est <- sd(LR)#パッケージはこっちで計算している
  
  n <- sqrt(length(likelihoods))#サンプルサイズ由来
  return (n*mean(LR)/omega_est)#これが統計量のはず
  #正ならば，対数正規分布
}

load(file = here("mergedf_01.rda"))

# 結果を保存するベクトル
test_stats <- numeric(nrow(df))
vuongst <- numeric(nrow(df))

for (i in 1:nrow(df)) {
  data <- df$mean_income[[i]]
  
  # 桁数を確認してスケーリング（オーバーフロー防止）
  max_digits <- nchar(as.character(floor(max(data, na.rm = TRUE))))
  if (max_digits > 13) {
    scale_factor <- 10^(max_digits - 13)
    data <- data / scale_factor
  }
  
  # パワーローモデルとログ正規モデルを適合
  pl_model <- conpl$new(data)
  est <- estimate_xmin(pl_model, xmax = max(data))
  pl_model$setXmin(est)
  
  pl_model_ln <- conlnorm$new(data)
  pl_model_ln$setXmin(pl_model$getXmin())
  pl_model_ln$setPars(estimate_pars(pl_model_ln))
  
  # 比較
  compd <- compare_distributions(pl_model, pl_model_ln)
  #test_stats[i] <- compd$test_statistic
  test_stats[i] <- compd$p_one_sided
  
  # ---- 中低所得者フィット ----
  xmin <- pl_model$xmin  # 高所得者の閾値
  #censdata <- data.frame(
  #  left  = ifelse(data < xmin, data, xmin),  
  #  right = ifelse(data < xmin, data, NA)
  #)
  censdata <- data.frame(#よくわからんでこれでやってみる
    left  = data,   # すべてのデータを left に設定
    right = ifelse(data >= xmin, NA, data)  # xmin 以上のデータは右検閲 (NA)
  )
  
  gamma_fit <- fitdistcens(censdata, "gamma",
    optim.method = "L-BFGS-B",   # 制約付き推定
  lower = c(1e-6, 1e-6),       # shape, scale > 0
  start = list(shape = 2, scale = mean(data[data < xmin])/2) # 初期値
  )
  ln_fit <- fitdistcens(censdata, "lnorm",
    optim.method = "Nelder-Mead",
  start = list(
    meanlog = log(mean(data[data < xmin])),
    sdlog   = sd(log(data[data < xmin]))
  )
)
  #SANN
  vuongst[i] <- 1 - pnorm(vuong(censdata, ln_fit, gamma_fit))
}
# dfに列として追加
df$pl_vs_lnorm_stat <- test_stats
df$lnorm_vs_gamma_stat <- vuongst



# ヒートマップの作成
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


# ヒートマップの作成
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
#赤だとガンマ分布，青だと対数正規