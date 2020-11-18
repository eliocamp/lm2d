#' Regress a 2d field
#'
#' Regresses a 2D scalar field into one time series using the method descibed by
#' DelSole and Yang (2011).
#'
#' @param x a formula to build the scalar field (see Details)
#' @param y vector
#' @param data data
#' @param method string with method
#' @param max_eof maximum number of principal components to use
#' @param k_fold number of folding for crossvalidation
#' @param alpha alpha number in lasso
#' @param seed seed used
#' @param verbose whether to talk a lot
#'
#' @details
#' The implementation from DelSole and Yang (2011) is `method = "cv"`. It uses
#' k-fold crossvalidation to select the first n principal components that
#' minimise the crossvalidation error when regressing `y` on `x`.
#'
#' The `lasso` method, uses LASSO/ridge rergession regularization to
#' automatically filter the principal components. The `alpha` parameter controls
#' the relative weight of the LASSO and ridge penalities (`alpha = 1` is the
#' LASSO penalty and `alpha = 0` is the ridge penalty).
#'
#' The `naive` method simply uses the first `max_eof` principal componets.
#'
#' @return
#' a list with elements:
#' \describe{
#'   \item{field}{the regression field}
#'   \item{fit}{the coefficients in EOF space}
#'   \item{summary}{a dataframe with r-squared, p.value and other aspects of the fit}
#'   \item{call}{the matched call}
#'   }
#'
#' @references
#' DelSole, T., & Yang, X. (2011). Field Significance of Regression Patterns. Journal of Climate, 24(19), 5094–5107. https://doi.org/10.1175/2011JCLI4105.1
#' @export
lm2d <- function(x, y, data = NULL,
                 method = c("cv", "lasso", "naive"),
                 max_eof = Inf,
                 k_fold = 10,
                 alpha = 1,
                 seed = 42,
                 verbose = FALSE) {
   fit_function <- match.fun(paste0("fit_", method[1]))

   maybe_message <- function(text) {
      if (verbose) {
         message(text)
      }
   }

   maybe_message("Parsing data")
   # Housekeeping. Getting data and transforming it to matrix
   x <- enrich_formula(x)

   y_name <- deparse(substitute(y))
   data <- data_from_formula(formula = x,
                             data = data,
                             extra.vars = y_name)

   g <- tidy2matrix(data, x$dims, c(x$value.var, y_name), fill = NULL)
   y <- g$matrix[[2]][, 1]

   g$matrix <- g$matrix[[1]]

   if (length(g$matrix) < nrow(data)) {
      stop(paste("The formula", as.character(x), "does not identify an unique observation for each cell."))
   }

   N <- length(g$rowdims[[1]])


   maybe_message("Computing EOF")
   g$matrix <- scale(g$matrix, scale = FALSE)

   if (max_eof < 1) {
      max_eof <- round(max_eof*min(nrow(g$matrix), ncol(g$matrix)))
   } else {
      max_eof <- min(max_eof, nrow(g$matrix), ncol(g$matrix))
   }

   eof <- smart_svd(g$matrix, max_eof)
   eof$d <- eof$d[seq_len(max_eof)]

   params <- list(max_eof = max_eof,
                  k_fold = k_fold,
                  alpha = alpha,
                  seed = seed)


   maybe_message("Fitting")
   fit <- fit_function(eof$u*sqrt(N), y, params)
   non_zero <- which(fit$coef != 0)

   E <- eof$v%*%diag(eof$d, nrow = length(eof$d))/sqrt(N)
   coef_real <- (fit$coef / var(y) * apply(eof$u*sqrt(N), 2, var)) %*% t(E)

   M <- length(non_zero)
   N <- length(y)
   f.statistic <-  fit$r2/(1 - fit$r2)*(N - M - 1)/M
   set(g$coldims, NULL, y_name, c(coef_real))
   model <- list(field = g$coldims,
                 fit = fit$coef,
                 summary = data.frame(r2 = fit$r2,
                                      f.statistic = f.statistic,
                                      p.value = pf(f.statistic, N, N - M - 1, lower.tail = FALSE),
                                      non_zero = M,
                                      N  = N),
                 call = match.call()
   )
   class(model) <-  "lm2d"
   return(model)
}

