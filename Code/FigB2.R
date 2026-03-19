################################################################################
#
#   Filename    :    FigB2.R
#   Project     :    BiomJ article "Robust HIV Viral Dynamics: A Nonlinear Mixed-Effects Framework for
#                    Heavy-Tailed Data with Informative Dropout"
#   Authors     :    Yu-Chen Yang and Tsung-I Lin and Luis M. Castro and Wan-Lun Wang
#   Date        :    14.03.2026
#   Purpose     :    produce Figure B.2 for Simulation
#
#   Input data files  :  None
#   Output data files :  Data_and_Code/results/FigureB2.eps
#
#   R Version   :    R-4.3.1
#   Required R packages : ggplot2; ggtext; cowplot; rlang
#
################################################################################
element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL) {
  structure(
    c(
      element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}

element_grob.element_textbox_highlight <- function(element, label = "", ...) {
  if (label %in% element$hi.labels) {
    element$fill <- element$hi.fill %||% element$fill
    element$colour <- element$hi.col %||% element$colour
    element$box.colour <- element$hi.box.col %||% element$box.colour
  }
  NextMethod()
}
############################################################################
############################################################################  sim 1
############################################################################
PATH1 <- paste(PATH, "/Data/Simulation/SS-simulation-t25/", sep = "")

realpara <- colMeans(as.matrix(read.table(paste(PATH1, "SIM1/realpara.txt", sep = ""), na.strings = "NA", sep = "")))
MAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

names(realpara) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")
colnames(MNAR.est1) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")

realt <- realpara
# realn = realpara[-9]

colMeans(MCAR.est2)
colMeans(MAR.est2)
colMeans(MNAR.est2)
# MNAR.est2[,2] = MNAR.est2[,2]-0.023
value <- c(
  apply((t(MAR.est1[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est1)[1],
  apply((t(MAR.est2[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est2)[1],
  apply((t(MAR.est3[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est3)[1],
  apply((t(MAR.est4[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est4)[1],
  apply((t(MAR.est5[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est5)[1],
  apply((t(MCAR.est1[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est1)[1],
  apply((t(MCAR.est2[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est2)[1],
  apply((t(MCAR.est3[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est3)[1],
  apply((t(MCAR.est4[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est4)[1],
  apply((t(MCAR.est5[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est5)[1],
  apply((t(MNAR.est1) - realt)^2, 1, sum) / dim(MNAR.est1)[1],
  apply((t(MNAR.est2) - realt)^2, 1, sum) / dim(MNAR.est2)[1],
  apply((t(MNAR.est3) - realt)^2, 1, sum) / dim(MNAR.est3)[1],
  apply((t(MNAR.est4) - realt)^2, 1, sum) / dim(MNAR.est4)[1],
  apply((t(MNAR.est5) - realt)^2, 1, sum) / dim(MNAR.est5)[1]
)

beta_labels <- c(
  expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
  expression(paste(sigma^2)), expression(paste(d)), expression(paste(nu))
)

beta_labels2 <- c(
  expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
  expression(paste(sigma^2)), expression(paste(d)), expression(paste(nu)),
  expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2]))
)

rep(gl(6, 1, labels = beta_labels), 10)
rep(gl(10, 1, labels = beta_labels2), 5)
ppc <- data.frame(
  value = value,
  sample = c(rep(rep(c("25", "50", "100", "200", "400"), each = 6), 2), rep(rep(c("25", "50", "100", "200", "400"), each = 10), 1)),
  beta = factor(c(
    as.character(rep(gl(6, 1, labels = beta_labels), 10)),
    as.character(rep(gl(10, 1, labels = beta_labels2), 5))
  )),
  G = c(rep(c("MAR", "MCAR"), each = 6 * 5), rep(c("MNAR"), each = 10 * 5))
)

ppc1.1 <- rbind(cbind(ppc, model = "tLME"))



ppc1.1$sample <- factor(ppc1.1$sample, c("25", "50", "100", "200", "400"))
ppc1.1$beta <- factor(ppc1.1$beta, levels = c(
  "paste(beta[1])", "paste(beta[2])", "paste(beta[3])",
  "paste(sigma^2)",
  "paste(d)", "paste(nu)",
  "paste(alpha[0][0])", "paste(alpha[0][1])", "paste(alpha[1])", "paste(alpha[2])"
))
ppc1.1$G <- factor(ppc1.1$G, c("MCAR", "MAR", "MNAR"))
ppc1.1$model <- factor(ppc1.1$model, levels = c("tLME"))

# ppc1.1$GG = rep()
ppc1.1$GG <- c(c(rep(c("MAR t", "MCAR t"), each = 6 * 5), rep(c("MNAR t"), each = 10 * 5)))


ppc1.1 <- cbind(ppc1.1, missing = "25% dropout")


############################################################################
############################################################################  sim 1
############################################################################
PATH1 <- paste(PATH, "/Data/Simulation/SS-simulation-t50/", sep = "")

realpara <- colMeans(as.matrix(read.table(paste(PATH1, "SIM1/realpara.txt", sep = ""), na.strings = "NA", sep = "")))
MAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

names(realpara) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")
colnames(MNAR.est1) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")

realt <- realpara
# realn = realpara[-9]

colMeans(MCAR.est2)
colMeans(MAR.est2)
colMeans(MNAR.est2)
# MNAR.est2[,2] = MNAR.est2[,2]-0.023
value <- c(
  apply((t(MAR.est1[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est1)[1],
  apply((t(MAR.est2[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est2)[1],
  apply((t(MAR.est3[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est3)[1],
  apply((t(MAR.est4[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est4)[1],
  apply((t(MAR.est5[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est5)[1],
  apply((t(MCAR.est1[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est1)[1],
  apply((t(MCAR.est2[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est2)[1],
  apply((t(MCAR.est3[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est3)[1],
  apply((t(MCAR.est4[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est4)[1],
  apply((t(MCAR.est5[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est5)[1],
  apply((t(MNAR.est1) - realt)^2, 1, sum) / dim(MNAR.est1)[1],
  apply((t(MNAR.est2) - realt)^2, 1, sum) / dim(MNAR.est2)[1],
  apply((t(MNAR.est3) - realt)^2, 1, sum) / dim(MNAR.est3)[1],
  apply((t(MNAR.est4) - realt)^2, 1, sum) / dim(MNAR.est4)[1],
  apply((t(MNAR.est5) - realt)^2, 1, sum) / dim(MNAR.est5)[1]
)

beta_labels <- c(
  expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
  expression(paste(sigma^2)), expression(paste(d)), expression(paste(nu))
)

beta_labels2 <- c(
  expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
  expression(paste(sigma^2)), expression(paste(d)), expression(paste(nu)),
  expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2]))
)

rep(gl(6, 1, labels = beta_labels), 10)
rep(gl(10, 1, labels = beta_labels2), 5)
ppc <- data.frame(
  value = value,
  sample = c(rep(rep(c("25", "50", "100", "200", "400"), each = 6), 2), rep(rep(c("25", "50", "100", "200", "400"), each = 10), 1)),
  beta = factor(c(
    as.character(rep(gl(6, 1, labels = beta_labels), 10)),
    as.character(rep(gl(10, 1, labels = beta_labels2), 5))
  )),
  G = c(rep(c("MAR", "MCAR"), each = 6 * 5), rep(c("MNAR"), each = 10 * 5))
)

ppc2.1 <- rbind(cbind(ppc, model = "tLME"))



ppc2.1$sample <- factor(ppc1.1$sample, c("25", "50", "100", "200", "400"))
ppc2.1$beta <- factor(ppc1.1$beta, levels = c(
  "paste(beta[1])", "paste(beta[2])", "paste(beta[3])",
  "paste(sigma^2)",
  "paste(d)", "paste(nu)",
  "paste(alpha[0][0])", "paste(alpha[0][1])", "paste(alpha[1])", "paste(alpha[2])"
))
ppc2.1$G <- factor(ppc1.1$G, c("MCAR", "MAR", "MNAR"))
ppc2.1$model <- factor(ppc1.1$model, levels = c("tLME"))

ppc2.1$GG <- c(c(rep(c("MAR t", "MCAR t"), each = 6 * 5), rep(c("MNAR t"), each = 10 * 5)))


ppc2.1 <- cbind(ppc2.1, missing = "50% dropout")


############################################################################
############################################################################  sim 3
############################################################################
PATH1 <- paste(PATH, "/Data/Simulation/SS-simulation-t75/", sep = "")

realpara <- colMeans(as.matrix(read.table(paste(PATH1, "SIM1/realpara.txt", sep = ""), na.strings = "NA", sep = "")))
MAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est1 <- as.matrix(read.table(paste(PATH1, "SIM1/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est2 <- as.matrix(read.table(paste(PATH1, "SIM2/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est3 <- as.matrix(read.table(paste(PATH1, "SIM3/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est4 <- as.matrix(read.table(paste(PATH1, "SIM4/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MCAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MCAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]
MNAR.est5 <- as.matrix(read.table(paste(PATH1, "SIM5/MNAR.para.est.txt", sep = ""), na.strings = "NA", sep = ""))[, -1]

names(realpara) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")
colnames(MNAR.est1) <- c("Beta1", "Beta2", "Beta3", "sigma", "d", "nu", "alpha00", "alpha01", "alpha1", "alpha2")

realt <- realpara

colMeans(MCAR.est2)
colMeans(MAR.est2)
colMeans(MNAR.est2)
value <- c(
  apply((t(MAR.est1[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est1)[1],
  apply((t(MAR.est2[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est2)[1],
  apply((t(MAR.est3[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est3)[1],
  apply((t(MAR.est4[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est4)[1],
  apply((t(MAR.est5[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MAR.est5)[1],
  apply((t(MCAR.est1[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est1)[1],
  apply((t(MCAR.est2[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est2)[1],
  apply((t(MCAR.est3[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est3)[1],
  apply((t(MCAR.est4[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est4)[1],
  apply((t(MCAR.est5[, 1:6]) - realt[1:6])^2, 1, sum) / dim(MCAR.est5)[1],
  apply((t(MNAR.est1) - realt)^2, 1, sum) / dim(MNAR.est1)[1],
  apply((t(MNAR.est2) - realt)^2, 1, sum) / dim(MNAR.est2)[1],
  apply((t(MNAR.est3) - realt)^2, 1, sum) / dim(MNAR.est3)[1],
  apply((t(MNAR.est4) - realt)^2, 1, sum) / dim(MNAR.est4)[1],
  apply((t(MNAR.est5) - realt)^2, 1, sum) / dim(MNAR.est5)[1]
)

beta_labels <- c(
  expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
  expression(paste(sigma^2)), expression(paste(d)), expression(paste(nu))
)

beta_labels2 <- c(
  expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
  expression(paste(sigma^2)), expression(paste(d)), expression(paste(nu)),
  expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2]))
)

rep(gl(6, 1, labels = beta_labels), 10)
rep(gl(10, 1, labels = beta_labels2), 5)
ppc <- data.frame(
  value = value,
  sample = c(rep(rep(c("25", "50", "100", "200", "400"), each = 6), 2), rep(rep(c("25", "50", "100", "200", "400"), each = 10), 1)),
  beta = factor(c(
    as.character(rep(gl(6, 1, labels = beta_labels), 10)),
    as.character(rep(gl(10, 1, labels = beta_labels2), 5))
  )),
  G = c(rep(c("MAR", "MCAR"), each = 6 * 5), rep(c("MNAR"), each = 10 * 5))
)


ppc3.1 <- rbind(cbind(ppc, model = "tLME"))

ppc3.1$sample <- factor(ppc1.1$sample, c("25", "50", "100", "200", "400"))
ppc3.1$beta <- factor(ppc1.1$beta, levels = c(
  "paste(beta[1])", "paste(beta[2])", "paste(beta[3])",
  "paste(sigma^2)",
  "paste(d)", "paste(nu)",
  "paste(alpha[0][0])", "paste(alpha[0][1])", "paste(alpha[1])", "paste(alpha[2])"
))
ppc3.1$G <- factor(ppc1.1$G, c("MCAR", "MAR", "MNAR"))
ppc3.1$model <- factor(ppc1.1$model, levels = c("tLME"))
ppc3.1$GG <- c(c(rep(c("MAR t", "MCAR t"), each = 6 * 5), rep(c("MNAR t"), each = 10 * 5)))

ppc3.1 <- cbind(ppc3.1, missing = "75% dropout")

ppcTotal <- rbind(ppc1.1, ppc2.1, ppc3.1)

ppcTotal <- ppcTotal[ppcTotal$model == "tLME", ]
ppcTotal <- ppcTotal[ppcTotal$G == "MNAR", ]
ppcTotal$missing <- factor(ppcTotal$missing, levels = c("25% dropout", "50% dropout", "75% dropout"))
ppcTotal$missing <- factor(ppcTotal$missing)

mse.sim2 <- ggplot(data = ppcTotal, aes(x = sample, y = value, colour = missing, group = missing)) +
  geom_line(aes(linetype = missing), linewidth = 1) +
  geom_point(aes(shape = missing), size = 4, stroke = 1.5) +
  xlab("sample size N") +
  ylab("MSE") +
  facet_wrap(~beta, labeller = "label_parsed", scales = "free_y", ncol = 3, dir = "h") +
  guides(
    shape = guide_legend(reverse = F, title = NULL, ncol = 3),
    linetype = guide_legend(reverse = F, title = NULL, ncol = 3),
    color = guide_legend(reverse = F, title = NULL, ncol = 3)
  ) +
  scale_shape_manual(values = c(2, 0, 1)) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  theme(legend.position = c(0.67, 0.1)) + 
  theme(text = element_text(size = 23)) +
  theme(
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm")
  ) +
  ylab("MSE") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

source(paste(PATH, "/function/multiplot.R", sep = ""))
layout <- matrix(c(1), nrow = 1, byrow = TRUE)
multiplot(plotlist = list(mse.sim2), layout = layout)

postscript(paste0(PATH, "/Result/FigureB2.eps"), width = 12, height = 10, paper = "special")
multiplot(plotlist = list(mse.sim2), layout = layout)
dev.off()
