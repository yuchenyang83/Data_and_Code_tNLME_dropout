################################################################################
#
#   Filename    :    FigB1.R
#   Project     :    "Robust HIV Viral Dynamics: A Nonlinear Mixed-Effects Framework for
#                    Heavy-Tailed Data with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    14.03.2026
#   Purpose     :    produce Figure B.1 for Simulation
#
#   Input data files  :  None
#   Output data files :  Data_and_Code/results/FigureB1.pdf
#
#   R Version   :    R-4.3.1
#   Required R packages : ggplot2; ggtext; mvtnorm; nlme
#
################################################################################
set.seed(8970)
source(paste0(PATH, "/function/simulation/simulate_dropout_study.R"))


# Parameter setting:
n <- 100
p <- 3
si <- 10
q <- 1
g <- 1
Beta <- matrix(c(27, 2, -0.5), ncol = g)
DD <- matrix(0, q, q)
DD[1, 1] <- 2
sigma <- 5
nu <- 5
Phi <- 1e-6

### 75 dropout
alpha <- c(-2.96, -0.69, -0.30, 0.30) ## intercept; treat; y_si-1; y_si
para <- list(alpha = alpha, Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, nu = nu)

cor.type <- "UNC"
gen.Data <- gen.tnlmm(n, para, cor.type = "UNC", si, q)
Data <- gen.Data$Data
Ymat <- gen.Data$Ymat
b <- gen.Data$b

gen.Data.miss <- add.MNAR.tNLMM(gen.Data$Data, gen.Data$Ymat, b, si, para$alpha)
Data.miss <- gen.Data.miss$Data.na
Data.miss$R <- rep(0, nrow(Data.miss))
Data.miss$R[is.na(Data.miss$Var1)] <- 1
Data <- Data.miss
table(Data$treat[Data$Time == 9]) / (n / 2)


Data1 <- data.frame(x = 0:9, y = mu.fn(Beta, 0:9))
matrix(c(27, 1, -0.5))

library(ggplot2)
k1 <- ggplot(data = Data[-which(Data$R == 1), ], aes(x = Time, y = Var1, group = as.factor(Subject))) +
  geom_line() +
  geom_point(size = 3) +
  xlab("Time") +
  ylab(NULL) +
  # geom_line(data = Data1, aes(x = 0:9, y = mu.fn(matrix(c(20, 0.5, -0.3)), (0:9*3)-10)+6, group = 1), colour = "red", size = 4) +
  theme(legend.position = "top") + ## 圖標位置
  scale_color_discrete(labels = c("Placebo", "Drug")) +
  scale_shape_discrete(labels = c("Placebo", "Drug")) +
  scale_x_continuous(breaks = 0:9) +
  theme(text = element_text(size = 20)) +
  theme(
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm")
  ) +
  theme(legend.title = element_blank())
k1


k2 <- ggplot(data = Data, aes(x = as.factor(Time), y = Var1, grop = as.factor(Time))) +
  geom_boxplot(na.rm = T) +
  xlab(NULL) +
  ylab(NULL) +
  theme(legend.position = "top") + ## 圖標位置
  scale_color_discrete(labels = c("Placebo", "Drug")) +
  scale_shape_discrete(labels = c("Placebo", "Drug")) +
  theme(text = element_text(size = 20)) +
  theme(
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm")
  ) +
  theme(legend.title = element_blank())
k2

k3 <- ggplot(data = Data[-which(is.na(Data$Var1)), ], aes(x = as.factor(Time), fill = as.factor(treat))) +
  geom_bar(position = position_dodge()) +
  xlab("Time") +
  ylab("Number of observed responses") +
  geom_text(
    stat = "count", aes(label = paste0(round(after_stat(count) / n, 3) * 100, "%")), vjust = 1.6, color = "black", size = 5,
    position = position_dodge(0.9)
  ) +
  theme(legend.position = "top") +
  scale_fill_discrete(name = NULL, labels = c("Placebo", "Drug")) +
  theme(text = element_text(size = 20))
k3

source(paste0(PATH, "/function/multiplot.R"))
layout <- matrix(c(1, 2, 1, 3), nrow = 2, byrow = TRUE)
multiplot(plotlist = list(k1, k2, k3), layout = layout)

pdf(paste0(PATH, "/Result/FigureB1.pdf"), width = 14, height = 10)
multiplot(plotlist = list(k1, k2, k3), layout = layout)
dev.off()
