################################################################################
#
#   Filename    :    FigB3.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Figure B.3 for the additional simulation study by
#                    comparing the MSEs of parameter estimates from the NLME
#                    and tNLME models under the MNAR mechanism with nearly
#                    Gaussian data generation across sample sizes and dropout
#                    levels
#
#   Input data files  :  Data_and_Code/Data/simulation/SS-simulation-t25-50/SIM1-SIM5/
#                         Data_and_Code/Data/simulation/SS-simulation-t50-50/SIM1-SIM5/
#                         Data_and_Code/Data/simulation/SS-simulation-t75-50/SIM1-SIM5/
#
#   Output data files :  Data_and_Code/Result/FigureB3.pdf
#
#   R Version   :    R-4.6.0
#   Required R packages : ggplot2; grid
#
################################################################################

################################################################################
# 25% dropout
################################################################################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
PATH1 <- paste(PATH, "/Data/simulation/SS-simulation-t25-50/", sep = "")

realpara <- colMeans(as.matrix(read.table(paste(PATH1, "SIM1/realpara.txt", sep = ""), na.strings = "NA", sep = "")))

t.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

n.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

names(realpara) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")

colnames(t.est1) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")
colnames(t.est2) <- colnames(t.est1)
colnames(t.est3) <- colnames(t.est1)
colnames(t.est4) <- colnames(t.est1)
colnames(t.est5) <- colnames(t.est1)

colnames(n.est1) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "alpha00", "alpha01", "alpha1", "alpha2")
colnames(n.est2) <- colnames(n.est1)
colnames(n.est3) <- colnames(n.est1)
colnames(n.est4) <- colnames(n.est1)
colnames(n.est5) <- colnames(n.est1)

realt <- realpara[c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]

value <- c(
  apply((t(n.est1[, c("Beta1", "Beta2", "Beta3", "alpha00", "Beta2", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est1)[1],
  apply((t(n.est2[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est2)[1],
  apply((t(n.est3[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est3)[1],
  apply((t(n.est4[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est4)[1],
  apply((t(n.est5[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est5)[1],
  apply((t(t.est1[, c("Beta1", "Beta2", "Beta3", "alpha00", "Beta2", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est1)[1],
  apply((t(t.est2[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est2)[1],
  apply((t(t.est3[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est3)[1],
  apply((t(t.est4[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est4)[1],
  apply((t(t.est5[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est5)[1]
)

beta_labels <- c(expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])), expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2])))

ppc1 <- data.frame(value = value, sample = rep(rep(c("25", "50", "100", "200", "400"), each = 7), 2), beta = factor(as.character(rep(gl(7, 1, labels = beta_labels), 10))), model = rep(c("NLME", "tNLME"), each = 7 * 5), missing = "25% dropout")

################################################################################
# 50% dropout
################################################################################

PATH1 <- paste(PATH, "/Data/simulation/SS-simulation-t50-50/", sep = "")

realpara <- colMeans(as.matrix(read.table(paste(PATH1, "SIM1/realpara.txt", sep = ""), na.strings = "NA", sep = "")))

t.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

n.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

names(realpara) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")

colnames(t.est1) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")
colnames(t.est2) <- colnames(t.est1)
colnames(t.est3) <- colnames(t.est1)
colnames(t.est4) <- colnames(t.est1)
colnames(t.est5) <- colnames(t.est1)

colnames(n.est1) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "alpha00", "alpha01", "alpha1", "alpha2")
colnames(n.est2) <- colnames(n.est1)
colnames(n.est3) <- colnames(n.est1)
colnames(n.est4) <- colnames(n.est1)
colnames(n.est5) <- colnames(n.est1)

realt <- realpara[c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]

value <- c(
  apply((t(n.est1[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est1)[1],
  apply((t(n.est2[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est2)[1],
  apply((t(n.est3[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est3)[1],
  apply((t(n.est4[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est4)[1],
  apply((t(n.est5[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est5)[1],
  apply((t(t.est1[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est1)[1],
  apply((t(t.est2[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est2)[1],
  apply((t(t.est3[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est3)[1],
  apply((t(t.est4[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est4)[1],
  apply((t(t.est5[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est5)[1]
)

ppc2 <- data.frame(value = value, sample = rep(rep(c("25", "50", "100", "200", "400"), each = 7), 2), beta = factor(as.character(rep(gl(7, 1, labels = beta_labels), 10))), model = rep(c("NLME", "tNLME"), each = 7 * 5), missing = "50% dropout")

################################################################################
# 75% dropout
################################################################################

PATH1 <- paste(PATH, "/Data/simulation/SS-simulation-t75-50/", sep = "")

realpara <- colMeans(as.matrix(read.table(paste(PATH1, "SIM1/realpara.txt", sep = ""), na.strings = "NA", sep = "")))

t.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
t.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

n.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
n.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.fit.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

names(realpara) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")

colnames(t.est1) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")
colnames(t.est2) <- colnames(t.est1)
colnames(t.est3) <- colnames(t.est1)
colnames(t.est4) <- colnames(t.est1)
colnames(t.est5) <- colnames(t.est1)

colnames(n.est1) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "alpha00", "alpha01", "alpha1", "alpha2")
colnames(n.est2) <- colnames(n.est1)
colnames(n.est3) <- colnames(n.est1)
colnames(n.est4) <- colnames(n.est1)
colnames(n.est5) <- colnames(n.est1)

realt <- realpara[c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]

value <- c(
  apply((t(n.est1[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est1)[1],
  apply((t(n.est2[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est2)[1],
  apply((t(n.est3[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est3)[1],
  apply((t(n.est4[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est4)[1],
  apply((t(n.est5[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(n.est5)[1],
  apply((t(t.est1[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est1)[1],
  apply((t(t.est2[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est2)[1],
  apply((t(t.est3[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est3)[1],
  apply((t(t.est4[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est4)[1],
  apply((t(t.est5[, c("Beta1", "Beta2", "Beta3", "alpha00", "alpha01", "alpha1", "alpha2")]) - realt)^2, 1, sum) / dim(t.est5)[1]
)

ppc3 <- data.frame(value = value, sample = rep(rep(c("25", "50", "100", "200", "400"), each = 7), 2), beta = factor(as.character(rep(gl(7, 1, labels = beta_labels), 10))), model = rep(c("NLME", "tNLME"), each = 7 * 5), missing = "75% dropout")

################################################################################
# Combine three dropout settings
################################################################################

ppcTotal <- rbind(ppc1, ppc2, ppc3)

ppcTotal$sample <- factor(ppcTotal$sample, levels = c("25", "50", "100", "200", "400"))
ppcTotal$beta <- factor(ppcTotal$beta, levels = c("paste(beta[1])", "paste(beta[2])", "paste(beta[3])", "paste(alpha[0][0])", "paste(alpha[0][1])", "paste(alpha[1])", "paste(alpha[2])"))
ppcTotal$model <- factor(ppcTotal$model, levels = c("NLME", "tNLME"))
ppcTotal$missing <- factor(ppcTotal$missing, levels = c("25% dropout", "50% dropout", "75% dropout"))


################################################################################
# Figure B.3
################################################################################

mse.t.normal <- ggplot(data = ppcTotal, aes(x = sample, y = value, colour = missing, shape = model, linetype = model, group = interaction(model, missing))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.4, stroke = 1.2, fill = NA) +
  facet_wrap(~beta, labeller = "label_parsed", scales = "free_y", ncol = 3) +
  scale_y_continuous(labels = function(x) format(x, trim = TRUE, scientific = FALSE, digits = 4)) +
  scale_colour_manual(name = NULL, values = c("25% dropout" = "#F8766D", "50% dropout" = "#00BA38", "75% dropout" = "#619CFF")) +
  scale_shape_manual(name = NULL, values = c("NLME" = 1, "tNLME" = 0)) +
  scale_linetype_manual(name = NULL, values = c("NLME" = "solid", "tNLME" = "dotted")) +
  labs(x = "sample size N", y = "MSE") +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 12)),
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),
    axis.text.x = element_text(size = 13, colour = "grey30"),
    axis.text.y = element_text(size = 13, colour = "grey30"),
    strip.background = element_rect(fill = "grey80", colour = "grey40", linewidth = 0.5),
    strip.text = element_text(size = 15),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.5),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.3),
    panel.spacing = unit(0.12, "cm"),
    legend.position = c(0.84, 0.16),
    legend.justification = c(0.5, 0.5),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.width = unit(1.4, "cm"),
    legend.key.height = unit(0.65, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  ) +
  guides(
    colour = guide_legend(order = 1, title = NULL),
    shape = guide_legend(order = 2, title = NULL),
    linetype = guide_legend(order = 2, title = NULL)
  )

mse.t.normal

pdf(paste0(PATH, "/Result/FigureB3.pdf"), width = 14, height = 10)
print(mse.t.normal)
dev.off()


