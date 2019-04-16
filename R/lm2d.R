lm2d <- function(x, y, data = NULL,
                 method = c("cv", "lasso", "naive"),
                 max_eof = Inf,
                 k_fold = 10,
                 alpha = 1,
                 seed = 42,
                 verbose = FALSE) {
   fit_function <- match.fun(paste0("fit_", method[1]))

   maybe_message("Parsing data")
   # Housekeeping. Getting data and transforming it to matrix
   x <- enrich_formula(x)
   force(y)
   y_name <- deparse(substitute(y))
   data <- data_from_formula(formula = x,
                             data = data,
                             extra.vars = y_name)
   g <- tidy2matrix(data, x$dims, x$value.var, fill = NULL)

   if (length(g$matrix) < nrow(data)) {
      stop(paste("The formula", as.character(x), "does not identify an unique observation for each cell."))
   }

   N <- length(g$rowdims[[1]])

   y <- tidy2matrix(data, x$dims, y_name, fill = NULL)$matrix[, 1]

   maybe_message("Computing EOF")
   # g$matrix <- scale(g$matrix, scale = FALSE)

   if (max_eof < 1) {
      max_eof <- round(max_eof*min(nrow(g$matrix), ncol(g$matrix)))
   } else {
      max_eof <- min(max_eof, nrow(g$matrix) - 1, ncol(g$matrix) - 1)
   }

   eof <- smart_svd(g$matrix, max_eof)
   eof$d <- eof$d[seq_len(max_eof)]

   maybe_message("Fitting")
   fit <- fit_function(eof$u*sqrt(N), y)
   non_zero <- which(fit$coef != 0)

   E <- eof$v%*%diag(eof$d, nrow = length(eof$d))/sqrt(N)
   coef_eof <- fit$coef
   fit$coef <- (fit$coef / var(y) * apply(eof$u*sqrt(N), 2, var)) %*% t(E)

   M <- length(non_zero)
   N <- length(y)
   f.statistic <-  fit$r2/(1 - fit$r2)*(N - M - 1)/M
   set(g$coldims, NULL, y_name, c(fit$coef))
   model <- list(field = g$coldims,
                 coef_eof = coef_eof,
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

