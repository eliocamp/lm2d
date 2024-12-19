#' Fit methods
#'
#' Functions that do the actual fit.
#'
#' @param k_fold number of folds for crossvalidation.
#' @param alpha elasticnet mixing parameter passed to [glmnet::glmnet].
#' `alpha = 1` is the lasso penalty, and `alpha = 0` the ridge penalty.
#' @param seed seed for crossvalidation.
#'
#' @export
#' @rdname fits
fit_cv <- function(k_fold = 10) {
  force(k_fold)
  function(x, y) {
    N <- length(y)
    n_features <- ncol(x)
    k_fold <- k_fold
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
        error[test] <-  cbind(1, x[test, seq_len(f), drop = FALSE]) %*% stats::coef(fit) - y[test]
      }
      mean_error[f] <- mean(error^2)
      sd_error[f] <- stats::sd(error^2)/sqrt(length(error))
    }

    min_error <- which.min(mean_error)

    features <- which((mean_error - sd_error)[seq_len(min_error)] < (mean_error + sd_error)[min_error])[1]
    fit <- stats::.lm.fit(cbind(1, x[, seq_len(features), drop = FALSE]), y)
    coef <- c(stats::coef(fit)[-1], rep(0, length = ncol(x) - features))
    r2 <- 1 - stats::var(fit$residuals)/stats::var(y)

    return(list(coef = coef,
                r2 = r2))
  }
}


#' @export
#' @rdname fits
fit_naive <- function() {
  function(x, y) {
    fit <- stats::.lm.fit(cbind(1, x), y)
    r2 <- 1 - stats::var(fit$residuals)/stats::var(y)

    return(list(coef = c(stats::coef(fit)[-1]),
                r2 = r2))
  }
}

#' @rdname fits
#' @export
fit_lasso <- function(k_fold = 10, alpha = 1, seed = 42) {
  force(k_fold)
  force(alpha)
  force(seed)
  function(x, y) {
    check_package("glmnet")
    fit <- withr::with_seed(seed, glmnet::cv.glmnet(x, y,
                                                    nfolds = k_fold,
                                                    lambda.min.ratio = 0.01,
                                                    alpha = alpha,
                                                    standardize = FALSE))
    lambda <- fit$lambda.1se
    r2 <- fit$glmnet.fit$dev.ratio[which(fit$glmnet.fit$lambda == lambda)]
    coef <- c(as.matrix(stats::coef(fit, s = lambda)))[-1]
    return(list(r2 = r2, coef = coef))
  }
}
