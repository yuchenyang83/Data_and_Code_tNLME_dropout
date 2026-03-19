################################################################################
#
#   Filename    :    Fig1.R
#   Project     :    BiomJ article "Robust HIV Viral Dynamics: A Nonlinear Mixed-Effects Framework for
#                    Heavy-Tailed Data with Informative Dropout"
#   Authors     :    Yu-Chen Yang and Tsung-I Lin and Luis M. Castro and Wan-Lun Wang
#   Date        :    14.03.2026
#   Purpose     :    produce Figure 1 for AIDS data
#
#   Input data files  :  Data_and_Code/Data/source/actg398.txt
#   Output data files :  Data_and_Code/results/Figure1.eps
#
#   R Version   :    R-4.3.1
#   Required R packages : ggplot2; ggtext; nlme; mvtnorm
#
################################################################################
actg398 <- read.table(paste0(PATH, "/Data/source/actg398.txt"), header = T)
table(actg398$calwk)

setdiff(1:481, actg398$patid[actg398$calwk == 0])
actg398$patid[actg398$calwk == 0]

actg398[actg398$patid == 56, ]
actg398[actg398$patid == 229, ]

actg398[actg398$patid == 56, ]$calwk <- c(0, 2, 8)
actg398[actg398$patid == 56, ]$txday <- actg398[actg398$patid == 56, ]$txday - 13
actg398[actg398$patid == 229, ]$calwk <- c(0, 2, 8, 16, 24)
actg398[actg398$patid == 229, ]$txday <- actg398[actg398$patid == 229, ]$txday - 15

table(actg398$trtarm[actg398$calwk == 0])
table(actg398$trtarm[actg398$calwk == 0], actg398$nnrti[actg398$calwk == 0])
actg398$trtarm[actg398$trtarm == 1] <- 0
actg398$trtarm[actg398$trtarm == 2] <- 0
actg398$trtarm[actg398$trtarm == 3] <- 0
actg398$trtarm[actg398$trtarm == 4] <- 1
table(actg398$trtarm[actg398$calwk == 0])

iid <- which(actg398$txday[actg398$calwk == 0] != 0)
for (i in 1:length(iid))
{
  actg398$txday[actg398$patid == iid[i]] <- actg398$txday[actg398$patid == iid[i]] - min(actg398$txday[actg398$patid == iid[i]])
}
table(actg398$calwk)
table(actg398$txday)

ll <- c(0, 2, 4, 8, 16, 24)
actg398$Subject <- actg398$patid
actg398$Time <- actg398$calwk

Subject <- unique(actg398$Subject)
N <- length(Subject)
ni <- numeric(N)
for (i in 1:N) ni[i] <- length(actg398$Time[actg398$Subject == Subject[i]])

actg398
actg398$D <- NA
for (i in 1:N)
{
  actg398$D[actg398$Subject == Subject[i]] <- max(actg398$Time[actg398$Subject == Subject[i]])
}

D.max <- 24
actg398.miss <- NULL
ll <- c(0, 2, 4, 8, 16, 24)
for (i in 1:N)
{
  if (unique(actg398[actg398$Subject == i, ]$D) < D.max) {
    actg398.i <- actg398[actg398$Subject == i, ]
    actg398.i <- rbind(actg398.i, actg398.i[dim(actg398.i)[1], ])
    kk <- which(ll == actg398.i$D[1])
    actg398.i$Time[length(actg398.i$Time)] <- ll[kk + 1]
    actg398.i$logrna[which(actg398.i$Time == ll[kk + 1])] <- NA
    actg398.i$D <- max(actg398.i$Time)
    actg398.i <- cbind(actg398.i, miss = 1)
    actg398.miss <- rbind(actg398.miss, actg398.i)
  } else {
    actg398.i <- actg398[actg398$Subject == i, ]
    actg398.i <- cbind(actg398.i, miss = 0)
    actg398.miss <- rbind(actg398.miss, actg398.i)
  }
}


actg398.miss[which(actg398.miss$txday == -1 | actg398.miss$txday == 0), ]
library(ggplot2)
kk <- ggplot(data = actg398.miss, aes(x = txday, y = logrna, grop = as.factor(patid))) +
  geom_line() +
  xlab("Day") +
  ylab("log(RNA)") +
  ggtitle("actg398") +
  theme(legend.position = "top") +
  theme(text = element_text(size = 20))
kk

## change to matrix
n <- max(ni)
N <- length(unique(actg398$Subject))

Subject <- unique(actg398.miss$Subject)
n <- length(Subject)
nj <- numeric(n)
for (i in 1:n) nj[i] <- length(actg398.miss$Dayt[actg398$Subject == Subject[i]])


ya1 <- ggplot(data = actg398.miss, aes(x = Time, y = logrna, group = as.factor(Subject))) +
  geom_line(aes(group = as.factor(Subject))) +
  geom_point() +
  xlab("Time") +
  ylab(NULL) +
  ggtitle("actg398") +
  theme(legend.position = "top") +
  theme(text = element_text(size = 20))
ya1


library(ggplot2)
ggplot() +
  geom_boxplot(
    data = actg398.miss[-which(is.na(actg398.miss$logrna)), ], aes(y = logrna, x = as.factor(calwk), group = as.factor(calwk)), alpha = 0,
    outlier.colour = "red", outlier.shape = 2, outlier.size = 3, outlier.alpha = 1, linewidth = 0.9
  ) +
  xlab("Month") +
  ylab(expression(sqrt("CD4"))) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 20))


library(ggplot2)
kk <- ggplot() +
  geom_line(data = actg398.miss[-which(is.na(actg398.miss$logrna)), ], aes(x = as.factor(calwk), y = logrna, group = as.factor(Subject))) +
  geom_point(data = actg398.miss[-which(is.na(actg398.miss$logrna)), ], aes(x = as.factor(calwk), y = logrna, group = as.factor(Subject))) +
  geom_boxplot(
    data = actg398.miss[-which(is.na(actg398.miss$logrna)), ], aes(y = logrna, x = as.factor(calwk), group = as.factor(calwk)), alpha = 0,
    outlier.colour = "red", colour = "red", lwd = 1, notchwidth = 0.5, outlier.shape = 2, outlier.size = 3, outlier.alpha = 1
  ) +
  xlab("Day") +
  ylab(expression(log[10]("RNA"))) +
  scale_x_discrete(breaks = ll) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 20))
kk


k5 <- ggplot(data = actg398.miss[-which(is.na(actg398.miss$logrna)), ], aes(x = as.factor(Time))) +
  geom_bar(position = position_dodge()) +
  xlab("Time") +
  ylab("Number of observed responses") +
  geom_text(
    stat = "count", aes(label = paste0(round(after_stat(count) / n, 3) * 100, "%")), vjust = -0.1, color = "red", size = 7,
    position = position_dodge(0.9)
  ) +
  scale_y_continuous(limits = c(0, 500)) +
  theme(legend.position = "top") +
  scale_fill_discrete(name = NULL, labels = c("Placebo", "Therapy")) +
  theme(text = element_text(size = 20))
k5

source(paste0(PATH, "/function/multiplot.R"))
layout <- matrix(c(1, 2), ncol = 1, byrow = TRUE)
multiplot(plotlist = list(kk, k5), layout = layout)

postscript(paste0(PATH, "/Result/Figure1.eps"), width = 10, height = 8, paper = "special")
multiplot(plotlist = list(kk, k5), layout = layout)
dev.off()
