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
load(paste0(PATH, "/Data/fit_result.RData"))
PATH <- getwd()
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

box_offset  <- 0.42
line_offset <- 0.42

obs <- obs %>%
  mutate(
    offset = ifelse(trtarm == "Therapy", -box_offset, box_offset),
    xpos   = calwk + offset
  )

sum_arm <- obs %>%
  group_by(trtarm, calwk) %>%
  summarise(med = mean(logrna, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    offset = ifelse(trtarm == "Therapy", -line_offset, line_offset),
    xpos   = calwk + offset
  )

out_df <- obs %>%
  group_by(calwk, trtarm, xpos) %>%
  summarise(
    q1    = quantile(logrna, 0.25, na.rm = TRUE),
    q3    = quantile(logrna, 0.75, na.rm = TRUE),
    iqr   = IQR(logrna, na.rm = TRUE),
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr,
    .groups = "drop"
  ) %>%
  left_join(obs, by = c("calwk", "trtarm", "xpos")) %>%
  filter(logrna < lower | logrna > upper) %>%
  mutate(outlier_legend = "Outlier")

col_arm <- c(
  "Therapy" = "#0072B2",
  "Placebo" = "#009E73",
  "Subject outlier" = "yellow3"
)
fill_arm <- c("Therapy" = "#F3DEDE", "Placebo" = "red3")
# fill_arm <- c("Therapy" = "#DCE6F1", "Placebo" = "#1F4E79")

legend_box <- data.frame(
  trtarm = factor(c("Therapy", "Placebo"), levels = c("Therapy", "Placebo")),
  x = day_levels[1],
  y = min(obs$logrna, na.rm = TRUE)
)

theme_pub <- theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14, color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank()
  )

## ------------------------------------------------------------
## fitted parameter estimates
## ------------------------------------------------------------
Beta  <- as.matrix(est.MNAR.ARp$para.est$Beta)
sigma <- est.MNAR.ARp$para.est$sigma
DD    <- est.MNAR.ARp$para.est$D
Phi   <- est.MNAR.ARp$para.est$Phi
b.hat <- as.matrix(est.MNAR.ARp$para.est$b)
nu    <- est.MNAR.ARp$para.est$nu

## pseudo-response after convergence
ytilde <- as.numeric(est.MNAR.ARp$y.c)
p <- length(Beta)
q <- nrow(DD)

ni <- numeric(N)
for (i in 1:N) ni[i] <- length(Data$Subject[Data$Subject == i])
n <- sum(ni)

na.ind <- which(is.na(as.vector(t(Data$Var1))))
cumsum.ni <- cumsum(ni)
cumsum.q <- cumsum(rep(q, N))
ni.o <- numeric(N)
for (i in 1:N) ni.o[i] <- sum(!is.na(Data$Var1[Data$Subject == i]))
cumsum.ni.o <- cumsum(ni.o)

A <- diag(p)
B <- matrix(c(
  0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, 0
), ncol = q)

# evaluate new log-likelihood
MU <- NULL
TXtilde <- 0
TZtilde <- matrix(0, ncol = N * q, nrow = n)
TLam <- TLam.inv <- TCor <- TCor.inv <- matrix(0, n, n)
cor.type = "ARp"
for (i in 1:N) {
  if (i == 1) {
    idx1 <- 1:cumsum.ni[1]
    idx2 <- 1:cumsum.q[1]
  } else {
    idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    idx2 <- (cumsum.q[i - 1] + 1):cumsum.q[i]
  }
  eta.ij <- A %*% Beta + B %*% b.hat[idx2, ]
  MU <- c(MU, mu.fn(eta.ij, Data$Time[Data$Subject == i], Data$nnrti[Data$Subject == i], Data$trtarm[Data$Subject == i]))
  dmu.ij <- dmu(eta.ij, Data$Time[Data$Subject == i], Data$nnrti[Data$Subject == i], Data$trtarm[Data$Subject == i])
  TXtilde <- rbind(TXtilde, dmu.ij %*% A)
  Ztilde <- dmu.ij %*% B
  TZtilde[idx1, ((i - 1) * q + 1):(i * q)] <- Ztilde
  Lam <- Ztilde %*% DD %*% t(Ztilde) + sigma * cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$Time[Data$Subject == i], ga = ga)
  TLam[idx1, idx1] <- Lam
  TLam.inv[idx1, idx1] <- solve(Lam)
  Cor <- cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$Time[Data$Subject == i], ga = ga)
  TCor[idx1, idx1] <- Cor
  TCor.inv[idx1, idx1] <- solve(Cor)
}
TXtilde <- TXtilde[-1, ]
# Ytilde = y.hat - MU + TXtilde %*% Beta + TZtilde %*% as.vector(b.hat)

ytilde.o <- c(ytilde[-na.ind])
Xtilde.o <- TXtilde[-na.ind, ]

TLam.oo <- TLam[-na.ind, -na.ind]
TLam.mo <- TLam[na.ind, -na.ind]
TLam.mm <- TLam[na.ind, na.ind]
# if (num.na == 1) TLam.mo <- t(TLam.mo)
TLam.oo.inv <- matrix(0, ncol = sum(ni.o), nrow = sum(ni.o))
for (i in 1:N) {
  if (i == 1) {
    idx1 <- 1:cumsum.ni.o[1]
  } else {
    idx1 <- (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]
  }
  TLam.oo.inv[idx1, idx1] <- solve(TLam.oo[idx1, idx1])
}

yo.cent <- ytilde.o - Xtilde.o %*% Beta
Delta <- t(yo.cent[1:cumsum.ni.o[1], ]) %*% TLam.oo.inv[1:cumsum.ni.o[1], 1:cumsum.ni.o[1]] %*% yo.cent[1:cumsum.ni.o[1], ]
for (i in 2:N) Delta <- c(Delta, t(yo.cent[(cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i], ]) %*% TLam.oo.inv[(cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i], (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]] %*% yo.cent[(cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i], ])
tau.hat = rep(NA, N)
for(i in 1:N){
  tau.hat[i] <- (nu + ni[i]) / (nu + Delta[i])
}
beta.cut = rep(NA, N)
for(i in 1:N){
  beta.cut[i] <- (1 + ni[i] / nu) *
    qbeta(0.025, shape1 = nu / 2, shape2 = ni[i] / 2)
}
out.beta <- tau.hat <= beta.cut

which(out.beta)
sum(out.beta)

out.beta = which(out.beta)

out.flag <- tau.hat <= beta.cut

which(out.flag)
sum(out.flag)

Delta[c(35,60 ,61,178,200,230,319,384,477)]
Delta[out.flag]
# out_subj <- c(35,60 ,61,178,200,230,319,384,477)
# tau.hat[c(35,60 ,61,178,200,230,319,384,477)]
tau.hat[out_subj]

## actual outlier subject IDs
out_subj <- which(out.flag)

obs_outsubj <- obs %>%
  dplyr::filter(Subject %in% out_subj) %>%
  dplyr::filter(!is.na(calwk), !is.na(logrna)) %>%
  dplyr::arrange(Subject, calwk) %>%
  dplyr::mutate(line_type = "Subject outlier")

## label at the last observed point of each outlier subject
lab_outsubj <- obs_outsubj %>%
  dplyr::group_by(Subject) %>%
  dplyr::slice_max(order_by = calwk, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    subj_lab = as.character(Subject),
    x_lab = calwk + 0.8,
    y_lab = logrna
  )

lab_outsubj[4, 17] = lab_outsubj[4, 17] + 0.1
lab_outsubj[10, 16] = lab_outsubj[10, 16] - 1
lab_outsubj[1, 17] = lab_outsubj[1, 17] - 0.05


p_top <- ggplot() +
  geom_line(data = obs, aes(x = calwk, y = logrna, group = Subject), color = "grey85", linewidth = 0.35) +
  geom_line(
    data = obs_outsubj,
    aes(x = calwk, y = logrna, group = Subject,
        color = line_type, linetype = line_type),
    linewidth = 0.8,
    lineend = "round",
    show.legend = TRUE
  ) +
  geom_boxplot(
    data = obs,
    aes(x = xpos, y = logrna, group = interaction(calwk, trtarm), fill = trtarm),
    width = 0.62,
    color = "grey30",
    linewidth = 1.1,
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  geom_point(
    data = out_df,
    aes(x = xpos, y = logrna, alpha = outlier_legend),
    shape = 2,
    color = "red3",
    size = 2.8,
    stroke = 0.9,
    show.legend = TRUE
  ) +
  geom_point(
    data = legend_box,
    aes(x = x, y = y, fill = trtarm),
    inherit.aes = FALSE,
    shape = 22,
    size = 0,
    stroke = 1.1,
    colour = "grey30",
    show.legend = TRUE
  ) +
  geom_text(
    data = lab_outsubj,
    aes(x = x_lab, y = y_lab, label = subj_lab),
    color = "#2F5BFF",
    size = 4,
    hjust = 0,
    vjust = 0.35,
    show.legend = FALSE,
    check_overlap = FALSE
  ) +
  geom_line(data = sum_arm, aes(x = xpos, y = med, group = trtarm, color = trtarm, linetype = trtarm), linewidth = 1.3) +
  geom_point(
    data = sum_arm,
    aes(x = xpos, y = med, color = trtarm, shape = trtarm),
    size = 5.2,
    stroke = 1.2
  ) +
  scale_x_continuous(
    breaks = day_levels,
    labels = day_levels
  ) +
  scale_color_manual(
    name = NULL,
    values = col_arm,
    breaks = c("Therapy", "Placebo", "Subject outlier")
  ) +
  scale_shape_manual(
    name = "Mean profile",
    values = c("Therapy" = 1, "Placebo" = 0),
    guide = "none"
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Therapy" = "solid",
      "Placebo" = "dashed",
      "Subject outlier" = "solid"
    ),
    guide = "none"
  ) +
  scale_fill_manual(
    name = "Boxplot",
    values = fill_arm
  ) +
  scale_alpha_manual(
    name = "Outlier",
    values = c("Outlier" = 1)
  ) +
  labs(
    x = "Week",
    y = expression(log[10](RNA))
  ) +
  guides(
    color = guide_legend(
      order = 1,
      title = NULL,
      keywidth = unit(0.5, "cm"),
      legend.spacing.x = unit(0.01, "cm"),
      legend.text = element_text(size = 9),
      override.aes = list(
        linetype = c("solid", "dashed", "solid"),
        shape = c(1, 0, NA),
        linewidth = c(1.3, 1.3, 0.8),
        size = c(5.2, 5.2, NA),
        stroke = c(1.2, 1.2, NA),
        fill = NA
      )
    ),
    fill = guide_legend(
      order = 2,
      title = NULL,
      override.aes = list(
        shape = 22,
        size = 6,
        stroke = 1.1,
        colour = "grey30",
        alpha = 1
      )
    ),
    alpha = guide_legend(
      order = 3,
      title = NULL,
      override.aes = list(
        shape = 2,
        color = "red3",
        size = 3.0,
        stroke = 0.9,
        linetype = 0
      )
    )
  ) +
  theme_pub +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.box = "horizontal",
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.spacing.x = unit(0.01, "cm"),
    legend.box.spacing = unit(0.02, "cm"),
    # legend.margin = margin(0, 0, 0, -2),
    legend.box.margin = margin(0, 0, 13, 0),
    legend.text = element_text(size = 9),
    plot.margin = margin(t = 18, r = 8, b = 8, l = 12)
  )
# p_top

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

n_total_all <- actg398.miss %>%
  distinct(Subject) %>%
  nrow()

bar_offset <- 0.22

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
    Time_id   = match(Time, day_levels),
    xpos      = Time_id + ifelse(trtarm == "Therapy", -bar_offset, bar_offset),
    prop_bar  = n_obs / n_total_all,
    label_bar = sprintf("%.1f%%", 100 * prop_bar)
  )

line_df <- actg398.miss %>%
  filter(!is.na(logrna)) %>%
  distinct(Subject, Time) %>%
  mutate(Time = as.numeric(Time)) %>%
  count(Time, name = "n_obs_all") %>%
  complete(Time = day_levels, fill = list(n_obs_all = 0)) %>%
  mutate(
    prop_all   = n_obs_all / n_total_all,
    label_line = sprintf("%.1f%%", 100 * prop_all),
    Time_id    = match(Time, day_levels)
  )


ymax_bottom <- max(c(bar_df$n_obs, line_df$n_obs_all)) * 1.08

pd_bar <- position_dodge2(width = 0.90, preserve = "single", padding = 0.20)

p_bottom <- ggplot() +
  geom_col(
    data = bar_df,
    aes(x = xpos, y = n_obs, fill = trtarm),
    width = 0.35,
    color = "grey35",
    linewidth = 1.2
  ) +
  geom_text(
    data = bar_df,
    aes(x = xpos, y = n_obs, label = label_bar),
    vjust = -0.75,
    hjust = 0.5,
    size = 2.5,
    color = "red3"
  ) +
  geom_text(
    data = bar_df,
    aes(x = xpos, y = n_obs, label = n_obs),
    vjust = 1.8,
    hjust = 0.5,
    size = 2.5,
    color = "black"
  ) +
  geom_line(
    data = line_df,
    aes(x = Time_id, y = n_obs_all, group = 1),
    color = "black",
    linewidth = 0.9
  ) +
  geom_point(
    data = line_df,
    aes(x = Time_id, y = n_obs_all),
    color = "black",
    size = 2.6
  ) +
  geom_text(
    data = line_df,
    aes(x = Time_id, y = n_obs_all, label = label_line),
    vjust = -0.8,
    hjust = 0.5,
    size = 4.6,
    color = "red3"
  ) +
  geom_text(
    data = line_df,
    aes(x = Time_id, y = n_obs_all, label = n_obs_all),
    vjust = 1.6,
    hjust = 0.5,
    size = 4.2,
    color = "black"
  ) +
  scale_fill_manual(values = fill_arm) +
  scale_x_continuous(
    breaks = seq_along(day_levels),
    labels = day_levels,
    name = "Week"
  ) +
  scale_y_continuous(
    limits = c(0, ymax_bottom),
    breaks = c(0, 144, 288, 432, 481),
    name = "Number of observed responses",
    sec.axis = sec_axis(
      trans = ~ . / n_total_all,
      breaks = c(0.30, 0.60, 0.90, 1.00),
      labels = c("30%", "60%", "90%", "100%"),
      name = "Observed proportion"
    )
  ) +
  theme_pub +
  theme(
    axis.title.y.left  = element_text(color = "black"),
    axis.text.y.left   = element_text(color = "black"),
    axis.ticks.y.left  = element_line(color = "black"),
    axis.title.y.right = element_text(color = "red3"),
    axis.text.y.right  = element_text(color = "red3"),
    axis.ticks.y.right = element_line(color = "red3")
  )+
  theme(
    legend.position = "top",
    # legend.justification = "left",
    legend.box = "horizontal",
    legend.key.width = unit(0.75, "cm"),
    legend.spacing.x = unit(0.12, "cm"),
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 18, r = 8, b = 8, l = 12)
  )

source(paste0(PATH, "/function/multiplot.R"))
layout <- matrix(c(2,1,2,1), ncol = 2, byrow = TRUE)
multiplot(plotlist = list(p_bottom, p_top), layout = layout)

postscript(paste0(PATH, "/Result/Figure1.eps"), width = 15, height = 10, paper = "special")
multiplot(plotlist = list(p_bottom, p_top), layout = layout)
dev.off()

