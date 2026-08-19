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

# Re-produce Figure B.1
source(file.path(PATH, "Code", "FigB1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Figure B.2
source(file.path(PATH, "Code", "FigB2.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Figure B.3
source(file.path(PATH, "Code", "FigB3.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Figure E.1
source(file.path(PATH, "Code", "FigE1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table 1
source(file.path(PATH, "Code", "Tab1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table 2
source(file.path(PATH, "Code", "Tab2.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table B.1
source(file.path(PATH, "Code", "TabB1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table B.2
source(file.path(PATH, "Code", "TabB2.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table B.3
source(file.path(PATH, "Code", "TabB3.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table D.1
source(file.path(PATH, "Code", "TabD1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table D.2
source(file.path(PATH, "Code", "TabD2.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Re-produce Table F.1
source(file.path(PATH, "Code", "TabF1.R"))
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

