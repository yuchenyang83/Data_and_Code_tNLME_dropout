################################################################################
#
#   Filename    :    Fig2.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    produce Figure 2 for the ACTG 398 data by summarizing the
#                    fixed-alpha2 sensitivity analysis under the selected tNLME
#                    model with Scenario (III), AR(1) errors, and MNAR dropout
#
#   Input data files  :  Data_and_Code/Data/fix_alpha/
#                         fit.t.III.ARp.MNAR_alpha*.RData
#
#   Intermediate file :  Data_and_Code/Data/fix_alpha/fixed_alpha.txt
#
#   Output data files :  Data_and_Code/Result/Figure2.pdf
#
#   R Version   :    R-4.6.0
#   Required R packages : ggplot2; dplyr; grid
#
################################################################################
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
DATA_DIR <- paste0(PATH, "/Data/fix_alpha")

ll <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 6)
ll <- c(ll, -c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 6))

alpha_tag <- function(x) {
  if (abs(abs(x) - 1e-4) < 1e-12) return(paste0(ifelse(x < 0, "-", ""), "1e-04"))
  format(x, scientific = FALSE, trim = TRUE, digits = 15)
}

file_names <- paste0("fit.t.III.ARp.MNAR_alpha", vapply(ll, alpha_tag, character(1)), ".RData")
file_paths <- file.path(DATA_DIR, file_names)

################################################################################
# 2. Check whether all 20 expected RData files are present
################################################################################
file_check <- data.frame(alpha = ll, file = file_names, exists = file.exists(file_paths), stringsAsFactors = FALSE)
# print(file_check, row.names = FALSE)

missing_files <- file_check$file[!file_check$exists]
actual_files <- list.files(DATA_DIR, pattern = "^fit\\.t\\.III\\.ARp\\.MNAR_alpha.*\\.RData$", full.names = FALSE)
extra_files <- setdiff(actual_files, file_names)

################################################################################
# 3. Load each fit and extract the eight beta estimates and SEs directly
#    from fit.t.III.ARp.MNAR$IM$out, then create fixed_alpha.txt
################################################################################
extract_one_fit <- function(rdata_file, alpha_value) {
  e <- new.env(parent = globalenv())
  loaded_names <- load(rdata_file, envir = e)

  expected_name <- "fit.t.III.ARp.MNAR"
  if (expected_name %in% loaded_names) {
    fit <- get(expected_name, envir = e)
  } else {
    candidate_names <- loaded_names[vapply(loaded_names, function(nm) {
      obj <- get(nm, envir = e)
      is.list(obj) && is.list(obj$IM) && !is.null(obj$IM$out)
    }, logical(1))]

    if (length(candidate_names) != 1) {
      stop("Cannot uniquely identify the fitted object in: ", basename(rdata_file),
           "\nObjects found: ", paste(loaded_names, collapse = ", "))
    }
    fit <- get(candidate_names, envir = e)
  }

  # IM$out is the displayed two-row table:
  # row 1 = estimates (EST), row 2 = standard errors.
  # Select only the eight columns whose column name is exactly "beta".
  im_out <- as.matrix(fit$IM$out)
  if (nrow(im_out) < 2 || is.null(colnames(im_out))) {
    stop("fit$IM$out must have at least two rows and column names in: ", basename(rdata_file))
  }

  beta_col <- which(colnames(im_out) == "beta")
  if (length(beta_col) != 8) {
    stop("Expected exactly 8 beta columns in fit$IM$out for: ", basename(rdata_file),
         " (found ", length(beta_col), ")")
  }

  beta_hat <- as.numeric(im_out[1, beta_col])
  beta_se <- as.numeric(im_out[2, beta_col])

  if (any(!is.finite(beta_hat)) || any(!is.finite(beta_se)) || any(beta_se < 0)) {
    stop("Invalid beta estimate or standard error in fit$IM$out for: ", basename(rdata_file))
  }

  z975 <- qnorm(0.975)
  data.frame(
    est.beta = beta_hat,
    est.upper = beta_hat + z975 * beta_se,
    est.lower = beta_hat - z975 * beta_se,
    alpha = alpha_value,
    beta = paste0("paste(beta[", seq_len(8), "])"),
    stringsAsFactors = FALSE
  )
}

fixed_alpha <- do.call(rbind, Map(extract_one_fit, file_paths, ll))
fixed_alpha$beta_id <- as.integer(sub(".*\\[([0-9]+)\\].*", "\\1", fixed_alpha$beta))
fixed_alpha <- fixed_alpha[order(fixed_alpha$alpha, fixed_alpha$beta_id), ]
fixed_alpha$beta_id <- NULL
rownames(fixed_alpha) <- NULL

fixed_alpha_file <- file.path(DATA_DIR, "fixed_alpha.txt")
write.table(fixed_alpha, fixed_alpha_file, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

################################################################################
# 4. Draw the fixed-alpha sensitivity figure
################################################################################
kk <- read.table(fixed_alpha_file, header = TRUE, na.strings = "NA")
kk <- as.data.frame(kk)
kk$beta <- gsub("^paste\\((beta\\[[0-9]+\\])\\)$", "\\1", as.character(kk$beta))
kk$alpha <- as.numeric(kk$alpha)

# Keep the same confidence-interval locations as the supplied Fig2.R.
# Because alpha = 0 is not fitted, the closest available value is selected there.
key_alpha <- c(-6, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 6)

kk_key <- kk %>%
  group_by(beta) %>%
  group_modify(~{
    do.call(rbind, lapply(key_alpha, function(a0) {
      .x %>% slice_min(order_by = abs(alpha - a0), n = 1, with_ties = FALSE)
    }))
  }) %>%
  ungroup() %>%
  distinct(beta, alpha, .keep_all = TRUE)

ya1 <- ggplot(kk, aes(x = alpha, y = est.beta)) +
  geom_line(aes(color = "Estimated value"), linewidth = 0.8) +
  geom_errorbar(
    data = kk_key,
    aes(ymin = est.lower, ymax = est.upper, color = "95% confidence interval"),
    width = 0.9, linetype = 1, linewidth = 0.65
  ) +
  geom_point(aes(color = "Estimated value"), size = 3) +
  facet_wrap(. ~ beta, labeller = label_parsed, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = c(-6, -4, -2, -1, 0, 1, 2, 4, 6)) +
  scale_y_continuous(labels = function(x) format(x, trim = TRUE, scientific = FALSE, digits = 4)) +
  scale_color_manual(
    values = c("Estimated value" = "black", "95% confidence interval" = "#084594"),
    breaks = c("Estimated value", "95% confidence interval"),
    labels = c("Estimated value", "95% Confidence interval")
  ) +
  labs(x = expression(alpha[2]), y = "Estimated fixed effects", color = NULL) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top",
    legend.key.width = unit(2.0, "cm"),
    legend.text = element_text(size = 22),
    strip.text.x = element_text(size = 24),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(
      reverse = FALSE,
      override.aes = list(linewidth = c(0.9, 0.8), shape = c(16, NA))
    )
  )

# print(ya1)

# postscript(paste0("D:/Data_and_Code", "/Result/Figur2.eps"), width = 15, height = 10, paper = "special")
# ya1
# dev.off()
PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

pdf(file.path(PATH, "Result", "Figure2.pdf"), width = 15, height = 10)
print(ya1)
dev.off()
