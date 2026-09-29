################################################################################
#
#   Filename    :    FigD1.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.09.2026
#   Purpose     :    produce Figure D.1 for the stopping-criterion sensitivity
#                    analysis by displaying the approximated observed-data
#                    log-likelihood trajectory of the selected Scenario (III)
#                    AR(1) tNLME model under MNAR over 5000 SAEM iterations
#
#   Input data files  :  Data_and_Code/Data/fit.t.III.ARp.MNAR5000.RData
#
#   Output data files :  Data_and_Code/Result/FigureD1.pdf
#
#   R Version   :    R-4.6.0
#   Required R packages : ggplot2
#
################################################################################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
load(paste0(PATH, "/Data/fit.t.III.ARp.MNAR5000.RData"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

lnL <- fit.t.III.ARp.MNAR.5000$model.inf$iter.lnL
ref <- -3742.522

plot.data <- data.frame(
  Iter = 0:(length(lnL) - 1),
  LogLik = lnL
)

base.breaks <- pretty(range(plot.data$LogLik, na.rm = TRUE), n = 6)

p <- ggplot(plot.data, aes(x = Iter, y = LogLik)) +
  geom_line(linewidth = 0.35) +
  geom_hline(yintercept = ref, color = "red", linewidth = 0.6) +
  annotate(
    "text",
    x = -Inf,
    y = ref,
    label = formatC(ref, format = "f", digits = 3),
    color = "red",
    fontface = "bold",
    hjust = -0.15,
    vjust = -0.8,
    size = 5
  ) +
  scale_y_continuous(
    breaks = base.breaks,
    labels = formatC(base.breaks, format = "f", digits = 0),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  scale_x_continuous(
    breaks = seq(0, 5000, by = 1000),
    limits = c(0, 5000),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Iter",
    y = "Observed log-likelihood value"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 19, margin = margin(t = 10)),
    axis.title.y = element_text(size = 19, margin = margin(r = 10)),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 20, 10, 55)
  )


# postscript(paste0(PATH, "/Result/FigureD1.eps"), width = 15, height = 10, paper = "special")
# print(p)
# dev.off()

pdf(paste0(PATH, "/Result/FigureD1.pdf"), width = 15, height = 10)
print(p)
dev.off()


