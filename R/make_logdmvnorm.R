#' Build a fast multivariate normal log-density evaluator with a fixed covariance (internal use only)
#'
#' Factorises \code{sigma} once and returns a closure evaluating the multivariate normal
#' log-density for arbitrary mean vectors. Mathematically identical to
#' \code{mvtnorm::dmvnorm(x, mean, sigma, log = TRUE)}, but avoids repeating the
#' Cholesky decomposition on every call.
#'
#' @param sigma A symmetric positive-definite covariance matrix.
#' @return A function of two arguments, \code{x} and \code{mean}, returning a log-density.
#' @keywords internal
make_logdmvnorm <- function(sigma){
  R <- chol(sigma)
  const <- -0.5 * ncol(sigma) * log(2 * pi) - sum(log(diag(R)))
  function(x, mean){
    z <- backsolve(R, x - mean, transpose = TRUE)
    const - 0.5 * sum(z * z)
  }
}
