################################################################################
#
#   Filename    :    master.R
#   Purpose     :    reproduce all figures and tables reported in the manuscript
#                    and Supplementary Materials from the stored data and fitted
#                    model objects
#
#   R Version   :    R-4.6.0
#
################################################################################

rm(list = ls())

library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales)
library(tidyr)
library(grid)
library(mvtnorm)
library(nlme)

PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Figure 1
source(file.path(PATH, "Code", "Fig1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Figure 2
source(file.path(PATH, "Code", "Fig2.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Figure D.1
source(file.path(PATH, "Code", "FigD1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Figure F.1
source(file.path(PATH, "Code", "FigF1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Figure F.2
source(file.path(PATH, "Code", "FigF2.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Figure F.3
source(file.path(PATH, "Code", "FigF3.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table 1
source(file.path(PATH, "Code", "Tab1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table 2
source(file.path(PATH, "Code", "Tab2.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table C.1
source(file.path(PATH, "Code", "TabC1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table C.2
source(file.path(PATH, "Code", "TabC2.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table E.1
source(file.path(PATH, "Code", "TabE1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table F.1
source(file.path(PATH, "Code", "TabF1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table F.2
source(file.path(PATH, "Code", "TabF2.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table F.3
source(file.path(PATH, "Code", "TabF3.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

