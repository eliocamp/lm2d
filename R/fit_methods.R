

fit_naive <- function(x, y, params) {
  zeros <- rep(0, length = ncol(x) - params$max_eof)
  fit <- stats::.lm.fit(cbind(1, x), y)
  r2 <- 1 - var(fit$residuals)/var(y)

  return(list(coef = c(coef(fit)[-1], zeros),
              r2 = r2))
}



fit_lasso <- function(x, y, params) {
  check_package("glmnet")
  fit <- withr::with_seed(params$seed, glmnet::cv.glmnet(x, y,
                                                  nfolds = params$k_fold,
                                                  lambda.min.ratio = 0.01,
                                                  alpha = params$alpha,
                                                  standardize = FALSE))
  lambda <- fit$lambda.1se
  r2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == lambda)]
  coef <- c(as.matrix(coef(fit, s = lambda)))[-1]
  return(list(r2 = r2, coef = coef))
}

#' @importFrom stats as.formula coef pf sd var .lm.fit
fit_cv <- function(x, y, params) {
  N <- length(y)
  n_features <- ncol(x)
  k_fold <- params$k_fold
  if (k_fold > N) {
    warning("K-folding for crossvalidation greater than N.",
            " Reducing to leave-one-out crossvalidation")
    k_fold <- N
  }

  N <- length(y)
  start <- floor(seq(1, N, by = N/k_fold))
  end   <-  c((start - 1)[-1], N)

  mean_error <- rep(0, length = n_features)
  sd_error <-  rep(0, length = n_features)
  for (f in seq_len(n_features)) {
    error <- rep(0, N)
    for (k in seq_along(k_fold)) {

      test <- seq.int(start[k], end[k])

      fit <- stats::.lm.fit(cbind(1, x[-test, seq_len(f), drop = FALSE]), y[-test])
      error[test] <-  cbind(1, x[test, seq_len(f), drop = FALSE]) %*% coef(fit) - y[test]
    }
    mean_error[f] <- mean(error^2)
    sd_error[f] <- sd(error^2)/sqrt(length(error))
  }

  min_error <- which.min(mean_error)

  features <- which((mean_error - sd_error)[seq_len(min_error)] < (mean_error + sd_error)[min_error])[1]
  fit <- .lm.fit(cbind(1, x[, seq_len(features), drop = FALSE]), y)
  coef <- c(coef(fit)[-1], rep(0, length = ncol(x) - features))
  r2 <- 1 - var(fit$residuals)/var(y)

  return(list(coef = coef,
              r2 = r2))
}
