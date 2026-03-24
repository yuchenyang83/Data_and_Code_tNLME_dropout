################################################################################
#
#   Filename    :    Fig2.R
#   Project     :    BiomJ article "Robust HIV Viral Dynamics: A Nonlinear Mixed-Effects Framework for
#                    Heavy-Tailed Data with Informative Dropout"
#   Authors     :    Yu-Chen Yang and Tsung-I Lin and Luis M. Castro and Wan-Lun Wang
#   Date        :    14.03.2026
#   Purpose     :    produce Figure 2 for AIDS data
#
#   Input data files  :  Data/fixed_alpha.txt
#   Output data files :  Result/Figure2.eps
#
#   R Version   :    R-4.3.1
#   Required R packages : ggplot2; dplyr; grid
#
################################################################################
library(ggplot2)
library(dplyr)
library(grid)

kk <- read.table(paste0(PATH, "/Data/fixed_alpha.txt"), header = TRUE, na.strings = "NA")

kk <- as.data.frame(kk)
kk$beta <- gsub("^paste\\((beta\\[[0-9]+\\])\\)$", "\\1", as.character(kk$beta))

kk$alpha <- as.numeric(kk$alpha)
key_alpha <- c(-6, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 6)

kk_key <- kk %>%
  group_by(beta) %>%
  group_modify(~{
    do.call(
      rbind,
      lapply(key_alpha, function(a0) {
        .x %>%
          slice_min(order_by = abs(alpha - a0), n = 1, with_ties = FALSE)
      })
    )
  }) %>%
  ungroup() %>%
  distinct(beta, alpha, .keep_all = TRUE)


ya1 <- ggplot(kk, aes(x = alpha, y = est.beta)) +
  geom_line(aes(color = "Estimated value"), linewidth = 0.8) +
  geom_errorbar(data = kk_key, aes(ymin = est.lower, ymax = est.upper, color = "95% confidence interval"),
                width = 0.9, linetype = 1, linewidth = 0.65) +
  ## all points plotted on top
  geom_point(aes(color = "Estimated value"), size = 3) +
  facet_wrap(. ~ beta, labeller = label_parsed, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = c(-6, -4, -2, -1, 0, 1, 2, 4, 6)) +
  scale_color_manual(
    values = c("Estimated value" = "black", "95% confidence interval" = "#084594"),
    breaks = c("Estimated value", "95% confidence interval"),
    labels = c("Estimated value", "95% Confidence interval")
  ) +
  labs(x = expression(alpha[2]), y = "Estimated fixed effects", color = NULL) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top",
    legend.key.width = unit(2.0, "cm"),
    legend.text = element_text(size = 22),
    strip.text.x = element_text(size = 24),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(
      reverse = FALSE,
      override.aes = list(
        linewidth = c(0.9, 0.8),
        shape = c(16, NA)
      )
    )
  )

print(ya1)

# postscript(paste0(PATH, "/Result/Figure2.eps"), width = 20, height = 10, paper = "special")
# print(ya1)
# dev.off()
