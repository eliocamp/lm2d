fit_naive <- function(X, y) {
  zeros <- rep(0, length = ncol(X) - max_eof)

  X <- X[, seq_len(max_eof), drop = FALSE]
  fit <- .lm.fit(cbind(1, X), y)
  r2 <- 1 - var(fit$residuals)/var(y)

  return(list(coef = c(coef(fit)[-1], zeros),
              r2 = r2))
}



fit_lasso <- function(X, y) {
  fit <- withr::with_seed(seed, glmnet::cv.glmnet(X, y, nfolds = k_fold,
                                                  lambda.min.ratio = 0.01,
                                                  alpha = alpha,
                                                  standardize = FALSE))
  lambda <- fit$lambda.1se
  r2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == lambda)]
  coef <- c(as.matrix(coef(fit, s = lambda)))[-1]
  return(list(r2 = r2, coef = coef))
}


fit_cv <- function(X, y) {
  n_features <- ncol(X)

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

      fit <- .lm.fit(cbind(1, X[-test, seq_len(f), drop = FALSE]), y[-test])
      error[test] <-  cbind(1, X[test, seq_len(f), drop = FALSE]) %*% coef(fit) - y[test]
    }
    mean_error[f] <- mean(error^2)
    sd_error[f] <- sd(error^2)/sqrt(length(error))
  }

  min_error <- which.min(mean_error)

  features <- which((mean_error - sd_error)[seq_len(min_error)] < (mean_error + sd_error)[min_error])[1]
  fit <- .lm.fit(cbind(1, X[, seq_len(features), drop = FALSE]), y)
  coef <- c(coef(fit)[-1], rep(0, length = ncol(X) - features))
  r2 <- 1 - var(fit$residuals)/var(y)

  return(list(coef = coef,
              r2 = r2))
}
