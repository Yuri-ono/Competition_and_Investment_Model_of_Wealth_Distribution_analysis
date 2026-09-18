library(here)
i_am("r_analysis/script/05_Plot_Heatmaps.R")
library(dplyr)
library(ggplot2)

# Heatmaps of the test results reported in the paper, drawn from the stored
# results in results/:
#   merge_01_vuong_pvalue_main.rda             Vuong tests (Figures 5 and 6)
#   mergedata_01_main_pvalue_toward_powerlaw.rda  power-law bootstrap p-values

plot_heatmap <- function(df, value) {
  df$label_stat <- sprintf("%.3f", df[[value]])
  ggplot(df, aes(x = factor(ω), y = factor(α), fill = .data[[value]])) +
    geom_tile(color = "white") +
    geom_text(aes(label = label_stat), size = 3) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0.5, name = "Test Statistic"
    ) +
    labs(
      x = expression(omega),
      y = expression(alpha)
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14)
    )
}

###Vuong test----------------------------

load(here("results", "merge_01_vuong_pvalue_main.rda"))

# Figure 5: power-law vs log-normal for the upper tail (x >= xmin)
# p < 0.1 favors power-law, p > 0.9 favors log-normal
plot_heatmap(df, "pl_vs_lnorm_stat")
ggsave(here("results", "Fig5_vuong_pl_vs_lnorm.pdf"), width = 7, height = 7, units = "in")
df$pl_vs_lnorm_stat |> summary()

# Figure 6: gamma vs log-normal for the bulk (x < xmin, right-censored at xmin)
# p < 0.1 favors log-normal, p > 0.9 favors gamma
plot_heatmap(df, "lnorm_vs_gamma_stat")
ggsave(here("results", "Fig6_vuong_lnorm_vs_gamma.pdf"), width = 7, height = 7, units = "in")
df$lnorm_vs_gamma_stat |> summary()

###Power-law bootstrap p-value----------------------------

load(here("results", "mergedata_01_main_pvalue_toward_powerlaw.rda"))

df <- df %>%
  mutate(p_over_0_1 = ifelse(p_value >= 0.1, "p >= 0.1", "p < 0.1"),
         p_value_label = sprintf("%.3f", p_value))

# Parameter sets where the power-law hypothesis is rejected (p < 0.1)
df %>% filter(p_value < 0.1) %>% select(ω, α, p_value)

ggplot(df, aes(x = factor(ω), y = factor(α), fill = p_over_0_1)) +
  geom_tile(color = "white") +
  geom_text(aes(label = p_value_label), size = 3) +
  scale_fill_manual(
    values = c("p < 0.1" = "blue", "p >= 0.1" = "red"),
    name = expression(italic(p)*"-value range"),
    labels = c(expression(italic(p) < 0.1), expression(italic(p) >= 0.1))
  ) +
  labs(x = expression(omega),
       y = expression(alpha)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12)
  )
ggsave(here("results", "powerlaw_pvalue.pdf"), width = 7, height = 7, units = "in")
