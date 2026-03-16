################################################################################
#                                                                                     
#   Filename    :    Fig2.R  												  
#   Project     :    "Robust HIV Viral Dynamics: A Nonlinear Mixed-Effects Framework for
#                    Heavy-Tailed Data with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    14.03.2026
#   Purpose     :    produce Figure 2 for AIDS data
#
#   Input data files  :  Data_and_Code/Data/fixed_alpha.txt
#   Output data files :  Data_and_Code/results/Figure2.pdf
#
#   R Version   :    R-4.3.1                                                              
#   Required R packages : ggplot2; ggtext; pROC; rlang; cowplot
#
################################################################################ 
### ### ### ### ### ### ### ### ### 
kk <- as.data.frame(read.table(paste0(PATH, "/Data/fixed_alpha.txt"), na.strings = "NA", sep = ""))

ya1 = ggplot(kk, aes(x=alpha, y=est.beta, color = "black")) + 
  geom_point(size = 2) + geom_line(linewidth = 1) +
  facet_wrap(.~beta, labeller = "label_parsed", scales = "free", ncol = 4) + 
  xlab(expression(paste(alpha[2]))) + 
  geom_line(data = kk, aes(x=alpha, y=est.upper, color = "red"), linetype = 2, linewidth = 1) +
  geom_line(data = kk, aes(x=alpha, y=est.lower, color = "red"), linetype = 2, linewidth = 1) +
  scale_x_continuous(breaks = c(-6,-4, -2, -1, 0, 1, 2, 4,6)) +
  theme(legend.key.width = unit(2, "cm")) +
  guides(shape = guide_legend(reverse=F,title=NULL),
         linetype = guide_legend(reverse=F,title=NULL),
         color = guide_legend(reverse=F,title=NULL),
         fill = guide_legend(reverse=F,title=NULL)) +
  theme(legend.position="top") +
  scale_y_continuous(name = "Estimated fixed effects") +
  scale_color_manual(values = c("black", "red"),labels = c("Estimated value", "95% Confidence bound")) + 
  theme(text = element_text(size=25), 
        strip.text.x = element_text(size = 25), 
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 25)) + 
  theme(axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(), 
        axis.title.y.right = element_text(color = "red", vjust = 2))
print(ya1)

source(paste(PATH, "/function/multiplot.R", sep = ""))
layout <- matrix(c(1), nrow = 1, byrow = TRUE)
multiplot(plotlist = list(ya1), layout = layout)

pdf(paste0(PATH, "/Result/Figure2.pdf"), width = 18, height = 8, paper = "special")
multiplot(plotlist = list(ya1), layout = layout)
dev.off()

