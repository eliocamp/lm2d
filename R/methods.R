summary.lm2d <- function(object, ..., colors = TRUE) {
  cat("Call:")
  cat("\n ", deparse(object$call), sep = "")
  cat("\n")
  object$summary <- lapply(object$summary, function(d) as.character(signif(d, 2)))
  nchars <- vapply(object$summary, nchar, FUN.VALUE = 1)

  cat("\n")
  cat("R-squared   = ", object$summary$r2, sep = "")
  cat("\n")
  cat("F-statistic = ", object$summary$f.statistic, sep = "")
  cat("\n")
  cat("p-value     = ", object$summary$p.value, sep = "")
  cat("\n\n")

  cat("Fit used ", object$summary$non_zero,
      " principal components (of the possible ",
      length(object$coef_eof), ") with coefficients:", sep = "")
  cat("\n")
  trunc <- "...[truncated]"

  w <- min(length(object$coef_eof), options()$width - nchar(trunc)-8)
  w <- seq_len(w)
  print(sparkbars::sparkbars(object$coef_eof[w], colors = colors))
  if (length(object$coef_eof) > options()$width - nchar(trunc)) {
    cat(trunc)
  }
}

print.lm2d <- function(x, ..., colors = TRUE) {
  summary(x, colors = colors)
}



autoplot.lm2d <- function(object, ...) {
  variables <- formula.tools::get.vars(object$call$x)
  xy <- variables[seq(length(variables)-1, length(variables))]
  estimate <- as.character(object$call$y)

  ggplot(object$field, aes_string(xy[1], xy[2])) +
    geom_raster(aes_string(fill = estimate), ...)
}


plot.lm2d <- autoplot.lm2d
