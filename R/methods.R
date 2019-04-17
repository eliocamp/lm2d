#' @export
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
      length(object$fit), ")")

  if (check_package("sparkbars", NULL)) {
    cat(" with coefficients:", sep = "")
    cat("\n")
    trunc <- "...[truncated]"

    w <- min(length(object$fit), options()$width - nchar(trunc)-8)
    w <- seq_len(w)

    print(sparkbars::sparkbars(object$fit[w], colors = colors))
    if (length(object$fit) > options()$width - nchar(trunc)) {
      cat(trunc)
    }}

}

#' @export
print.lm2d <- function(x, ..., colors = TRUE) {
  summary(x, colors = colors)
}



#' @export
plot.lm2d <- function(x, ...) {
  variables <- formula.tools::get.vars(x$call$x)
  xy <- variables[seq(length(variables)-1, length(variables))]
  estimate <- as.character(x$call$y)

  check_package("ggplot2")
  ggplot2::ggplot(x$field, ggplot2::aes_string(xy[1], xy[2])) +
    ggplot2::geom_raster(ggplot2::aes_string(fill = estimate), ...)
}

