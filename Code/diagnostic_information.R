################################################################################
#
#   Filename    :    diagnostic_information.R
#   Project     :    Biometrics article "Modeling HIV Viral Dynamics Using a
#                    Nonlinear Mixed-Effects Framework for Heavy-Tailed Data
#                    with Informative Dropout"
#   Authors     :    Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang
#   Date        :    18.08.2026
#   Purpose     :    construct general convergence diagnostics for a fitted
#                    NLME or tNLME model, including the approximated observed-
#                    data log-likelihood trajectory, parameter trace plots, and
#                    Metropolis-Hastings acceptance-rate summaries
#
#   Input data files  :  none; the fitted model object must already exist in
#                        the current R session
#
#   Output data files :  none
#
#   R Version   :    R-4.6.0
#   Required R packages : ggplot2
#
################################################################################

library(ggplot2)

diagnostic_information <- function(fit) {
  
  fit.expression <- substitute(fit)
  fit.name <- deparse(fit.expression)
  
  if (is.name(fit.expression) && !exists(fit.name, envir = parent.frame(), inherits = TRUE)) {
    stop(
      paste0(
        "The fitted object '", fit.name, "' is not loaded in the current R session. ",
        "Load the RData file containing this object first."
      )
    )
  }
  
  fit <- eval(fit.expression, envir = parent.frame())
  
  ##############################################################################
  # 1. Check required components
  ##############################################################################
  
  if (is.null(fit$model.inf$iter.lnL)) stop("The fitted object does not contain model.inf$iter.lnL.")
  if (is.null(fit$model.inf$loglik)) stop("The fitted object does not contain model.inf$loglik.")
  if (is.null(fit$Tpara)) stop("The fitted object does not contain Tpara.")
  if (is.null(fit$Taccept.rate)) stop("The fitted object does not contain Taccept.rate.")
  if (is.null(fit$IM$out)) stop("The fitted object does not contain IM$out.")
  
  Tpara.all <- as.matrix(fit$Tpara)
  IM.out <- as.matrix(fit$IM$out)
  Taccept <- as.matrix(fit$Taccept.rate)
  
  if (nrow(Tpara.all) < 2) stop("Tpara does not contain enough stored iterations.")
  if (ncol(Tpara.all) < ncol(IM.out)) stop("Tpara contains fewer columns than the estimated parameters in IM$out.")
  
  ##############################################################################
  # 2. Construct parameter names from IM$out
  ##############################################################################
  
  raw.names <- colnames(IM.out)
  
  if (is.null(raw.names)) raw.names <- paste0("Parameter_", seq_len(ncol(IM.out)))
  
  raw.names[is.na(raw.names) | raw.names == ""] <- "Parameter"
  lower.names <- tolower(raw.names)
  parameter.names <- raw.names
  
  beta.idx <- which(lower.names == "beta")
  if (length(beta.idx) > 0) parameter.names[beta.idx] <- paste0("beta", seq_along(beta.idx))
  
  d.idx <- which(lower.names == "d")
  if (length(d.idx) == 1) parameter.names[d.idx] <- "d"
  if (length(d.idx) == 3) parameter.names[d.idx] <- c("d11", "d12", "d22")
  if (length(d.idx) == 10) parameter.names[d.idx] <- c("d11", "d12", "d22", "d13", "d23", "d33", "d14", "d24", "d34", "d44")
  if (length(d.idx) > 0 && !length(d.idx) %in% c(1, 3, 10)) parameter.names[d.idx] <- paste0("d", seq_along(d.idx))
  
  sigma.idx <- which(grepl("sigma|sig", lower.names))
  if (length(sigma.idx) == 1) parameter.names[sigma.idx] <- "sigma^2"
  if (length(sigma.idx) > 1) parameter.names[sigma.idx] <- paste0("sigma", seq_along(sigma.idx))
  
  phi.idx <- which(lower.names == "phi")
  if (length(phi.idx) == 1) parameter.names[phi.idx] <- "phi"
  if (length(phi.idx) > 1) parameter.names[phi.idx] <- paste0("phi", seq_along(phi.idx))
  
  nu.idx <- which(lower.names == "nu")
  if (length(nu.idx) == 1) parameter.names[nu.idx] <- "nu"
  if (length(nu.idx) > 1) parameter.names[nu.idx] <- paste0("nu", seq_along(nu.idx))
  
  alpha.idx <- which(lower.names == "alpha")
  if (length(alpha.idx) == 1) parameter.names[alpha.idx] <- "alpha00"
  if (length(alpha.idx) == 4) parameter.names[alpha.idx] <- c("alpha00", "alpha01", "alpha02", "alpha1")
  if (length(alpha.idx) == 5) parameter.names[alpha.idx] <- c("alpha00", "alpha01", "alpha02", "alpha1", "alpha2")
  if (length(alpha.idx) > 0 && !length(alpha.idx) %in% c(1, 4, 5)) parameter.names[alpha.idx] <- paste0("alpha", seq_along(alpha.idx))
  
  parameter.names <- make.unique(parameter.names)
  
  ##############################################################################
  # 3. Match Tpara columns to the estimated parameters in IM$out
  #
  # Tpara may contain inactive/fixed placeholder columns for some correlation
  # structures or missingness mechanisms. Therefore, Tpara and IM$out are not
  # required to have the same number of columns.
  ##############################################################################
  
  final.trace <- apply(
    Tpara.all,
    2,
    function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) return(NA_real_)
      tail(x, 1)
    }
  )
  
  trace.range <- apply(
    Tpara.all,
    2,
    function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) return(0)
      diff(range(x))
    }
  )
  
  target.est <- as.numeric(IM.out[1, ])
  p <- length(target.est)
  q <- ncol(Tpara.all)
  
  if (p == q) {
    
    trace.columns <- seq_len(q)
    
  } else {
    
    active.tol <- sqrt(.Machine$double.eps) * pmax(1, abs(final.trace))
    active.columns <- which(is.finite(final.trace) & trace.range > active.tol)
    
    if (length(active.columns) == p) {
      
      trace.columns <- active.columns
      
    } else {
      
      ############################################################################
      # Ordered matching:
      # choose p columns from Tpara, preserving their original order, so that
      # their final stored values are as close as possible to the estimates in
      # IM$out. This allows inactive placeholder columns to be skipped.
      ############################################################################
      
      scale.est <- pmax(abs(target.est), 1)
      cost <- matrix(Inf, nrow = p, ncol = q)
      
      for (i in seq_len(p)) {
        cost[i, ] <- abs(final.trace - target.est[i]) / scale.est[i]
        cost[i, !is.finite(cost[i, ])] <- 1e12
      }
      
      dp <- matrix(Inf, nrow = p + 1, ncol = q + 1)
      take <- matrix(FALSE, nrow = p + 1, ncol = q + 1)
      dp[1, ] <- 0
      
      for (i in 2:(p + 1)) {
        for (j in 2:(q + 1)) {
          
          skip.cost <- dp[i, j - 1]
          take.cost <- dp[i - 1, j - 1] + cost[i - 1, j - 1]
          
          if (take.cost <= skip.cost) {
            dp[i, j] <- take.cost
            take[i, j] <- TRUE
          } else {
            dp[i, j] <- skip.cost
          }
        }
      }
      
      if (!is.finite(dp[p + 1, q + 1])) stop("Unable to match Tpara columns to IM$out parameters.")
      
      trace.columns <- integer(p)
      i <- p + 1
      j <- q + 1
      
      while (i > 1 && j > 1) {
        if (take[i, j]) {
          trace.columns[i - 1] <- j - 1
          i <- i - 1
          j <- j - 1
        } else {
          j <- j - 1
        }
      }
      
      if (any(trace.columns == 0)) stop("Unable to identify all parameter traces in Tpara.")
    }
  }
  
  Tpara <- Tpara.all[, trace.columns, drop = FALSE]
  
  trace.column.map <- data.frame(
    Parameter = parameter.names,
    IM_Estimate = target.est,
    Tpara_column = trace.columns,
    Tpara_final_value = final.trace[trace.columns],
    Difference = final.trace[trace.columns] - target.est,
    stringsAsFactors = FALSE
  )
  
  ##############################################################################
  # 4. General axis-label formatter
  ##############################################################################
  
  axis_number <- function(x) {
    
    x <- as.numeric(x)
    out <- rep("", length(x))
    finite <- is.finite(x)
    
    if (!any(finite)) {
      out[is.infinite(x) & x > 0] <- "Inf"
      out[is.infinite(x) & x < 0] <- "-Inf"
      return(out)
    }
    
    xf <- x[finite]
    ux <- sort(unique(xf))
    dx <- abs(diff(ux))
    dx <- dx[is.finite(dx) & dx > 0]
    
    if (length(dx) == 0) {
      
      max.abs <- max(abs(xf), na.rm = TRUE)
      
      if (!is.finite(max.abs) || max.abs == 0) {
        digits <- 0L
      } else if (max.abs >= 1) {
        digits <- 2L
      } else {
        digits <- min(8L, max(2L, as.integer(ceiling(-log10(max.abs))) + 2L))
      }
      
    } else {
      
      step <- min(dx)
      
      if (!is.finite(step) || step <= 0) {
        digits <- 3L
      } else if (step >= 1) {
        digits <- 0L
      } else {
        digits <- min(8L, max(1L, as.integer(ceiling(-log10(step))) + 1L))
      }
    }
    
    scientific <- xf != 0 & (abs(xf) < 1e-6 | abs(xf) >= 1e6)
    
    plain.index <- which(!scientific)
    sci.index <- which(scientific)
    
    formatted <- character(length(xf))
    
    if (length(plain.index) > 0) {
      
      fmt <- paste0("%.", digits, "f")
      plain <- sprintf(fmt, xf[plain.index])
      
      if (digits > 0) {
        plain <- sub("(\\.[0-9]*[1-9])0+$", "\\1", plain)
        plain <- sub("\\.0+$", "", plain)
      }
      
      plain[plain == "-0"] <- "0"
      formatted[plain.index] <- plain
    }
    
    if (length(sci.index) > 0) {
      formatted[sci.index] <- sprintf("%.2e", xf[sci.index])
    }
    
    out[finite] <- formatted
    out[is.infinite(x) & x > 0] <- "Inf"
    out[is.infinite(x) & x < 0] <- "-Inf"
    
    out
  }
  
  ##############################################################################
  # 5. Approximated observed-data log-likelihood trajectory
  ##############################################################################
  
  lnL <- as.numeric(fit$model.inf$iter.lnL)
  final.loglik <- as.numeric(fit$model.inf$loglik)
  
  if (!any(is.finite(lnL))) stop("model.inf$iter.lnL does not contain finite values.")
  
  loglik.data <- data.frame(
    Iteration = 0:(length(lnL) - 1),
    LogLik = lnL
  )
  
  loglik.range <- range(loglik.data$LogLik, finite = TRUE)
  loglik.breaks <- pretty(loglik.range, n = 6)
  
  loglik.plot <- ggplot(loglik.data, aes(x = Iteration, y = LogLik)) +
    geom_line(linewidth = 0.45) +
    scale_x_continuous(breaks = pretty(loglik.data$Iteration, n = 6)) +
    scale_y_continuous(
      breaks = loglik.breaks,
      labels = axis_number
    ) +
    labs(
      title = paste0(
        fit.name,
        ": observed-data log-likelihood trajectory; final = ",
        sprintf("%.3f", final.loglik)
      ),
      x = "Iteration",
      y = "Observed log-likelihood value"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  ##############################################################################
  # 6. Parameter trace plots
  ##############################################################################
  
  trace.data <- data.frame(
    Iteration = rep(0:(nrow(Tpara) - 1), times = ncol(Tpara)),
    Parameter = factor(rep(parameter.names, each = nrow(Tpara)), levels = parameter.names),
    Estimate = as.vector(Tpara),
    stringsAsFactors = FALSE
  )
  
  trace.plot <- ggplot(trace.data, aes(x = Iteration, y = Estimate)) +
    geom_line(linewidth = 0.30) +
    facet_wrap(~Parameter, scales = "free_y", ncol = 4) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 4)) +
    scale_y_continuous(labels = axis_number, breaks = function(x) pretty(x, n = 4)) +
    labs(
      title = paste0(fit.name, ": parameter trace plots"),
      x = "Iteration",
      y = "Parameter estimate"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 14),
      strip.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(0.7, "lines")
    )
  
  ##############################################################################
  # 7. Parameter estimate table from IM$out
  ##############################################################################
  
  parameter.table <- data.frame(
    Parameter = parameter.names,
    Estimate = as.numeric(IM.out[1, ]),
    stringsAsFactors = FALSE
  )
  
  if (nrow(IM.out) >= 2) parameter.table$SE <- as.numeric(IM.out[2, ])
  
  ##############################################################################
  # 8. Metropolis-Hastings acceptance-rate summaries
  ##############################################################################
  
  acceptance.by.component <- colMeans(Taccept, na.rm = TRUE)
  
  component.names <- colnames(Taccept)
  
  if (is.null(component.names)) {
    component.names <- paste0("Component_", seq_along(acceptance.by.component))
  }
  
  missing.component.names <- is.na(component.names) | component.names == ""
  
  if (any(missing.component.names)) {
    component.names[missing.component.names] <- paste0("Component_", which(missing.component.names))
  }
  
  acceptance.table <- data.frame(
    Component = component.names,
    Acceptance_rate = as.numeric(acceptance.by.component),
    stringsAsFactors = FALSE
  )
  
  acceptance.table <- acceptance.table[is.finite(acceptance.table$Acceptance_rate), , drop = FALSE]
  
  acceptance.active <- acceptance.by.component[
    is.finite(acceptance.by.component) & acceptance.by.component != 0
  ]
  
  if (length(acceptance.active) == 0) stop("No non-zero Metropolis-Hastings acceptance rates were found.")
  
  overall.acceptance <- mean(acceptance.active)
  
  diagnostic.table <- data.frame(
    Diagnostic = c(
      "Stored log-likelihood values",
      "Stored parameter-trace iterations",
      "Number of estimated parameters",
      "Number of stored Tpara columns",
      "Final approximated observed log-likelihood",
      "Overall MH acceptance rate",
      "Minimum active MH acceptance rate",
      "Maximum active MH acceptance rate",
      "Overall acceptance within 0.2-0.4"
    ),
    Value = c(
      as.character(length(lnL)),
      as.character(nrow(Tpara.all)),
      as.character(ncol(IM.out)),
      as.character(ncol(Tpara.all)),
      formatC(final.loglik, format = "f", digits = 3),
      formatC(overall.acceptance, format = "f", digits = 3),
      formatC(min(acceptance.active), format = "f", digits = 3),
      formatC(max(acceptance.active), format = "f", digits = 3),
      ifelse(overall.acceptance >= 0.2 & overall.acceptance <= 0.4, "Yes", "No")
    ),
    stringsAsFactors = FALSE
  )
  
  ##############################################################################
  # 9. Return diagnostics; nothing is printed automatically
  ##############################################################################
  
  invisible(list(
    loglik.plot = loglik.plot,
    trace.plot = trace.plot,
    diagnostic.table = diagnostic.table,
    parameter.table = parameter.table,
    acceptance.table = acceptance.table,
    acceptance.by.component = acceptance.by.component,
    overall.acceptance = overall.acceptance,
    trace.column.map = trace.column.map,
    IM.out = IM.out,
    parameter.names = parameter.names,
    iter.loglik = lnL,
    parameter.trace = Tpara
  ))
}

################################################################################
# Examples
#
# Example 1: Scenario III
#
# load("D:/Data_and_Code/Data/fit_III_result.RData")
#
# gg = diagnostic_information(fit.t.III.ARp.MNAR)
# gg$loglik.plot
# gg$trace.plot
# gg$diagnostic.table
# gg$parameter.table
# gg$acceptance.table
# gg$trace.column.map
#
#
# Example 2: another correlation structure
#
# gg = diagnostic_information(fit.t.III.UNC.MNAR)
# gg$loglik.plot
# gg$trace.plot
# gg$diagnostic.table
# gg$parameter.table
# gg$acceptance.table
# gg$trace.column.map
#
#
# Example 3: Scenario II
#
# load("D:/Data_and_Code/Data/fit_II_result.RData")
#
# gg = diagnostic_information(fit.t.II.ARp.MNAR)
# gg$loglik.plot
# gg$trace.plot
# gg$diagnostic.table
# gg$parameter.table
# gg$acceptance.table
# gg$trace.column.map
#
# Example 4: Scenario III 5000 iteration
#
# load("D:/Data_and_Code/Data/fit.t.III.ARp.MNAR5000.RData")
#
# gg = diagnostic_information(fit.t.III.ARp.MNAR)
# gg$loglik.plot
# gg$trace.plot
# gg$diagnostic.table
# gg$parameter.table
# gg$acceptance.table
# gg$trace.column.map
#
################################################################################