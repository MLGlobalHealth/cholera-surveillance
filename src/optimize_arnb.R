#' Autoregressive Negative Binomial likelihood function
#'
#' @param data The data to fit, numeric vector.
#' @param start_vals Starting values for the parameters
#' @param ... other parameters to be passed into the `optim` function
#'
#' @return the log likelihood of the AR negative binomial
#'
arnb_ll <- function(theta, data) {
  lls <- numeric(length(data))
  a <- theta[1]
  b <- theta[2]
  psi <- theta[3]
  
  lls[1] <- dnbinom(data[1], mu = exp(a) + 10e-25, size = exp(psi), log = TRUE)
  for (i in 2:(length(data))) {
    mu <- exp(a + b * data[i - 1])
    if (is.nan(mu)) {
      cat("a", a, "\n", "b", b, "\n", "mu", mu, "\n")
    }
    lls[i] <- dnbinom(data[i], mu = mu + 10e-25, size = exp(psi), log = TRUE)
  }
  -sum(lls)
}

arnb_fit <- function(
    data,
    start_vals = c(6, 0, 0),  # log-scale
    ...
) {
    out <- optim(
      start_vals, 
      arnb_ll, 
      data = data,
      ...
    )
    res <- list()
    res$theta <- out$par
    names(res$theta) <- c("Intercept", "AR1", "Size")
    res$cond_ll <- arnb_ll(theta = res$theta, data = data)
    res$ll <- -out$value
    class(res) <- 'AR1_NegBinom'
    res
}