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
#   Required R packages : ggplot2; ggtext; nlme; mvtnorm; dplyr; scales; tidyr
#
################################################################################
library(grid)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
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

## change to matrix
n <- max(ni)
N <- length(unique(actg398$Subject))

Subject <- unique(actg398.miss$Subject)
n <- length(Subject)
nj <- numeric(n)
for (i in 1:n) nj[i] <- length(actg398.miss$Dayt[actg398$Subject == Subject[i]])

## observed data only
obs <- actg398.miss %>%
  filter(!is.na(logrna)) %>%
  mutate(
    calwk  = as.numeric(calwk),
    trtarm = factor(trtarm, levels = c(0, 1), labels = c("Therapy", "Placebo"))
  )

day_levels <- ll

obs <- obs %>%
  mutate(
    offset = ifelse(trtarm == "Therapy", -0.24, 0.24),
    xpos   = calwk + offset
  )

sum_arm <- obs %>%
  group_by(trtarm, calwk) %>%
  summarise(med = mean(logrna, na.rm = TRUE), .groups = "drop") %>%
  mutate(offset = ifelse(trtarm == "Therapy", -0.4, 0.4), xpos   = calwk + offset)

col_arm  <- c("Therapy" = "#9E1F1F", "Placebo" = "#1F4E79")
fill_arm <- c("Therapy" = "#F3DEDE", "Placebo" = "#DCE6F1")

theme_pub <- theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14, color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )

p_top <- ggplot() +
  geom_line(data = obs, aes(x = calwk, y = logrna, group = Subject),
            color = "grey85",
            linewidth = 0.35) +
  geom_boxplot(data = obs, aes(x = xpos, y = logrna, group = interaction(calwk, trtarm), linetype = trtarm, fill = trtarm),
               width = 1.5,
               color = "grey30",
               linewidth = 1.1,
               show.legend = FALSE,
               outlier.shape = 2,
               outlier.fill = "red",
               outlier.colour = "red3",
               outlier.size = 2.8,
               outlier.stroke = 0.9,
               outlier.alpha = 1
  ) +
  geom_line(data = sum_arm, aes(x = xpos, y = med, group = trtarm, color = trtarm, linetype = trtarm), linewidth = 1.3) +
  geom_point(data = sum_arm, aes(x = xpos, y = med, color = trtarm, shape = trtarm), size = 5.2, stroke = 1.2) +
  scale_x_continuous(breaks = day_levels,labels = day_levels) +
  scale_color_manual(values = col_arm) +
  scale_shape_manual(values = c("Therapy" = 1, "Placebo" = 0)) +
  scale_fill_manual(values = fill_arm) +
  labs(x = "Week", y = expression(log[10](RNA))) +
  theme(legend.key.width = unit(2.0, "cm")) +
  guides(
    color = guide_legend(keywidth = unit(2.0, "cm")),
    shape = guide_legend(keywidth = unit(2.0, "cm")),
    linetype = guide_legend(keywidth = unit(2.0, "cm"))) + 
  theme_pub
# p_top

n_total_all <- actg398.miss %>%
  distinct(Subject) %>%
  nrow()

bar_df <- actg398.miss %>%
  filter(!is.na(logrna)) %>%
  distinct(Subject, trtarm, Time) %>%
  mutate(
    trtarm = factor(trtarm, levels = c(0, 1), labels = c("Therapy", "Placebo")),
    Time   = as.numeric(Time)
  ) %>%
  count(trtarm, Time, name = "n_obs") %>%
  complete(
    trtarm,
    Time = day_levels,
    fill = list(n_obs = 0)
  ) %>%
  mutate(
    Time_f    = factor(Time, levels = day_levels),
    prop_bar  = n_obs / n_total_all,
    label_bar = sprintf("%.1f%%", 100 * prop_bar)
  )

line_df <- actg398.miss %>%
  filter(!is.na(logrna)) %>%
  distinct(Subject, Time) %>%
  mutate(Time = as.numeric(Time)) %>%
  count(Time, name = "n_obs_all") %>%
  complete(Time = day_levels, fill = list(n_obs_all = 0)) %>%
  mutate(prop_all   = n_obs_all / n_total_all, label_line = sprintf("%.1f%%", 100 * prop_all),
         Time_f = factor(Time, levels = day_levels))


ymax_bottom <- max(c(bar_df$n_obs, line_df$n_obs_all)) * 1.08

p_bottom <- ggplot(bar_df, aes(x = Time_f, y = n_obs, fill = trtarm, linetype = trtarm)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = "grey35",
    linewidth = 1.2
  ) +
  geom_text(
    aes(label = label_bar),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 2.5,
    color = "red3"
  ) +
  geom_text(
    aes(y = n_obs, label = n_obs),
    position = position_dodge(width = 0.72),
    vjust = 1.8,
    size = 2.5,
    color = "black"
  ) +
  geom_line(
    data = line_df,
    aes(x = Time_f, y = n_obs_all, group = 1),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.9
  ) +
  geom_point(
    data = line_df,
    aes(x = Time_f, y = n_obs_all),
    inherit.aes = FALSE,
    color = "black",
    size = 2.6
  ) +
  geom_text(
    data = line_df,
    aes(x = Time_f, y = n_obs_all, label = label_line),
    inherit.aes = FALSE,
    vjust = -0.75,
    size = 4.6,
    color = "red3"
  ) +
  geom_text(
    data = line_df,
    aes(x = Time_f, y = n_obs_all, label = n_obs_all),
    inherit.aes = FALSE,
    vjust = 1.6,
    size = 4.2,
    color = "black"
  ) +
  scale_fill_manual(values = fill_arm) +
  scale_y_continuous(
    limits = c(0, ymax_bottom),
    breaks = c(0, 144, 288, 432, 481),
    name = "Number of observed responses",
    sec.axis = sec_axis(
      trans = ~ . / n_total_all,
      breaks = c(0.30, 0.60, 0.90, 1.00),
      labels = c("30%", "60%", "90%", "100%"),
      name = "Observed proportion")) +
  labs(x = "Week") +
  theme_pub +
  theme(
    axis.title.y.left  = element_text(color = "black"),
    axis.text.y.left   = element_text(color = "black"),
    axis.ticks.y.left  = element_line(color = "black"),
    axis.title.y.right = element_text(color = "red3"),
    axis.text.y.right  = element_text(color = "red3"),
    axis.ticks.y.right = element_line(color = "red3")
  )

source(paste0(PATH, "/function/multiplot.R"))
layout <- matrix(c(1, 2, 1, 2), ncol = 2, byrow = TRUE)
multiplot(plotlist = list(p_top, p_bottom), layout = layout)

# postscript(paste0(PATH, "/Result/Figure1.eps"), width = 12, height = 8, paper = "special")
# multiplot(plotlist = list(p_top, p_bottom), layout = layout)
# dev.off()

