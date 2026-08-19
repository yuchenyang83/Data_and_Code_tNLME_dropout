################################################################################
#
#   Filename    :    Fig1.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Figure 1 for the ACTG 398 data, including the
#                    longitudinal viral-load trajectories, visit-specific
#                    distributions, observed-response frequencies, and
#                    subject-level outliers identified under the selected
#                    tNLME model with Scenario (III), AR(1) errors, and
#                    MNAR dropout
#
#   Input data files  :  Data_and_Code/Data/source/actg398.txt
#                         Data_and_Code/Data/Figure1_outlier_subjects.txt
#
#   Output data files :  Data_and_Code/Result/Figure1.pdf
#
#   R Version   :    R-4.6.0
#   Required R packages : ggplot2; dplyr; ggrepel; scales; tidyr; grid
#
################################################################################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
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
fill_arm <- c("Therapy" = "#F3DEDE", "Placebo" = "#D95F4F")
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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    
    ## make axis labels closer to the black panel border
    axis.text.y.left   = element_text(color = "black", margin = margin(r = 1)),
    axis.title.y.left  = element_text(color = "black", margin = margin(r = 2)),
    axis.text.y.right  = element_text(margin = margin(l = 1)),
    axis.title.y.right = element_text(margin = margin(l = 2)),
    axis.text.x        = element_text(color = "black", margin = margin(t = 1)),
    axis.title.x       = element_text(margin = margin(t = 2)),
    
    ## shorten ticks so labels visually attach closer to the frame
    axis.ticks.length = unit(2, "pt")
  )

panel_box <- annotate(
  "rect",
  xmin = -Inf, xmax = Inf,
  ymin = -Inf, ymax = Inf,
  fill = NA,
  colour = "black",
  linewidth = 0.6
)


## ------------------------------------------------------------
## subject-level outliers extracted once from the selected fitted model
## ------------------------------------------------------------
outlier.file <- file.path(PATH, "Data", "Figure1_outlier_subjects.txt")
outlier.data <- read.table(outlier.file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!"Subject" %in% names(outlier.data)) stop("Figure1_outlier_subjects.txt must contain a Subject column.")
out_subj <- as.integer(outlier.data$Subject)

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
    x_lab = calwk + 0.4,
    y_lab = logrna
  )

# lab_outsubj[4, 17] = lab_outsubj[4, 17] + 0.1
# lab_outsubj[10, 16] = lab_outsubj[10, 16] - 1
# lab_outsubj[1, 17] = lab_outsubj[1, 17] - 0.05

legend_common <- theme(
  legend.position = "top",
  legend.justification = "center",
  legend.box.just = "center",
  legend.box = "horizontal",
  legend.key.width = unit(0.75, "cm"),
  legend.key.height = unit(0.65, "cm"),
  legend.spacing.x = unit(0.08, "cm"),
  legend.box.spacing = unit(0.02, "cm"),
  legend.box.margin = margin(0, 0, 6, 0),
  legend.text = element_text(size = 11)
)


p_top <- ggplot() +
  # geom_line(data = obs, aes(x = calwk, y = logrna, group = Subject), color = "grey85", linewidth = 0.35) +
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
    color = "#D95F4F",
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
  ggrepel::geom_text_repel(
    data = lab_outsubj,
    aes(x = x_lab, y = y_lab, label = subj_lab),
    color = "#2F5BFF",
    size = 4,
    direction = "y",
    nudge_x = 0.45,
    hjust = 0,
    segment.color = NA,
    box.padding = 0.25,
    point.padding = 0.15,
    max.overlaps = Inf,
    show.legend = FALSE
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
    labels = day_levels,
    limits = c(-1, 26.2)
  ) +
  coord_cartesian(xlim = c(-1, 26.2), ylim = c(0.5, 7.4), clip = "off") +
  # geom_segment(aes(x = -1,   xend = 26.2, y = 0.5, yend = 0.5),
  #              inherit.aes = FALSE, color = "black", linewidth = 0.6) +
  # geom_segment(aes(x = -1,   xend = 26.2, y = 7.4, yend = 7.4),
  #              inherit.aes = FALSE, color = "black", linewidth = 0.6) +
  # geom_segment(aes(x = -1,   xend = -1,   y = 0.5, yend = 7.4),
  #              inherit.aes = FALSE, color = "black", linewidth = 0.6) +
  # geom_segment(aes(x = 26.2, xend = 26.2, y = 0.5, yend = 7.4),
  #              inherit.aes = FALSE, color = "black", linewidth = 0.6) +
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
      # legend.text = element_text(size = 9),
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
        color = "#D95F4F",
        size = 3.0,
        stroke = 0.9,
        linetype = 0
      )
    )
  ) +
  theme_pub +
  legend_common +
  theme(
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

## observed counts by treatment group and week
obs_by_group <- actg398.miss %>%
  filter(!is.na(logrna)) %>%
  distinct(Subject, trtarm, Time) %>%
  mutate(
    trtarm = factor(trtarm, levels = c(0, 1),
                    labels = c("Therapy observed", "Placebo observed")),
    Time = as.numeric(Time)
  ) %>%
  count(Time, trtarm, name = "n_obs") %>%
  complete(
    Time = day_levels,
    trtarm = factor(c("Therapy observed", "Placebo observed"),
                    levels = c("Therapy observed", "Placebo observed")),
    fill = list(n_obs = 0)
  )

## total missing at each week
miss_by_week <- obs_by_group %>%
  group_by(Time) %>%
  summarise(
    n_obs_all = sum(n_obs),
    .groups = "drop"
  ) %>%
  mutate(
    n_miss = n_total_all - n_obs_all
  )

## stacked data
bar_stack <- bind_rows(
  obs_by_group %>%
    transmute(
      Time,
      Component = as.character(trtarm),
      Count = n_obs
    ),
  miss_by_week %>%
    transmute(
      Time,
      Component = "Missing",
      Count = n_miss
    )
) %>%
  mutate(
    Component = factor(
      Component,
      levels = c("Therapy observed", "Placebo observed", "Missing")
    ),
    Time_id = match(Time, day_levels),
    Prop = Count / n_total_all
  ) %>%
  arrange(Time_id, Component) %>%
  group_by(Time_id) %>%
  mutate(
    ymin = cumsum(Count) - Count,
    ymax = cumsum(Count),
    ymid = (ymin + ymax) / 2,
    label_pct = sprintf("%.1f%%", 100 * Prop),
    label_n = paste0(Count)
  ) %>%
  ungroup()

bar_label <- bar_stack %>%
  filter(Count > 0) %>%
  mutate(
    y_pct = ymid + 7,
    y_n   = ymid - 7
  )

fill_stack <- c(
  "Therapy observed" = fill_arm["Therapy"],
  "Placebo observed" = fill_arm["Placebo"],
  "Missing" = "white"
)

p_bottom <- ggplot() +
  geom_rect(
    data = bar_stack,
    aes(
      xmin = Time_id - 0.28,
      xmax = Time_id + 0.28,
      ymin = ymin,
      ymax = ymax,
      fill = Component
    ),
    color = "grey35",
    linewidth = 1.1
  ) +
  geom_text(
    data = bar_label,
    aes(x = Time_id, y = y_pct, label = label_pct),
    color = "#0072B2",
    size = 4
  ) +
  geom_text(
    data = bar_label,
    aes(x = Time_id, y = y_n, label = label_n),
    color = "black",
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Therapy observed" = unname(fill_arm["Therapy"]),
      "Placebo observed" = unname(fill_arm["Placebo"]),
      "Missing" = "white"
    ),
    breaks = c("Therapy observed", "Placebo observed", "Missing"),
    labels = c("Therapy observed", "Placebo observed", "Missing"),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq_along(day_levels),
    labels = day_levels,
    name = "Week",
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(0.00, 7.00), clip = "off") +
  scale_y_continuous(
    limits = c(0, n_total_all * 1.05),
    breaks = c(0, 144, 288, 432, 481),
    name = "Number of responses",
    sec.axis = sec_axis(
      transform = ~ . / n_total_all,
      breaks = c(0, 0.30, 0.60, 0.90, 1.00),
      labels = c("0%", "30%", "60%", "90%", "100%"),
      name = "Proportion"
    )
  ) +
  # geom_segment(aes(x = 0.45, xend = 6.55, y = 0, yend = 0),
  #              inherit.aes = FALSE, color = "black", linewidth = 0.6) +
  # geom_segment(aes(x = 0.45, xend = 6.55, y = 481 * 1.05, yend = 481 * 1.05),
  #              inherit.aes = FALSE, color = "black", linewidth = 0.6) +
  # geom_segment(aes(x = 0.45, xend = 0.45, y = 0, yend = 481 * 1.05),
  #              inherit.aes = FALSE, color = "black", linewidth = 0.6) +
  # geom_segment(aes(x = 6.55, xend = 6.55, y = 0, yend = 481 * 1.05),
  #              inherit.aes = FALSE, color = "black", linewidth = 0.6) +
  labs(fill = NULL) +
  theme_pub +
  legend_common +
  theme(
    axis.title.y.left  = element_text(color = "black", margin = margin(r = 2)),
    axis.text.y.left   = element_text(color = "black", margin = margin(r = 1)),
    axis.ticks.y.left  = element_line(color = "black"),
    axis.title.y.right = element_text(color = "#0072B2", margin = margin(l = 2)),
    axis.text.y.right  = element_text(color = "#0072B2", margin = margin(l = 1)),
    axis.ticks.y.right = element_line(color = "#0072B2"),
    plot.margin = margin(t = 18, r = 8, b = 8, l = 12)
  )

# p_bottom

source(paste0(PATH, "/function/multiplot.R"))
layout <- matrix(c(2,1,2,1), ncol = 2, byrow = TRUE)
# multiplot(plotlist = list(p_bottom, p_top), layout = layout)

# postscript(paste0(PATH, "/Result/Figur1.eps"), width = 15, height = 10, paper = "special")
# multiplot(plotlist = list(p_bottom, p_top), layout = layout)
# dev.off()

pdf(paste0(PATH, "/Result/Figure1.pdf"), width = 15, height = 10)
multiplot(plotlist = list(p_bottom, p_top), layout = layout)
dev.off()

