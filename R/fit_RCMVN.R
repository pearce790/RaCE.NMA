#' Fit a Bayesian Rank-Clustered Estimation model for Network Meta-Analysis (RaCE-NMA). Recommended for internal use only; use mcmc_RCMVN instead.
#'
#' This function fits a Bayesian RaCE-NMA model to data from a previous network meta-analysis.
#'
#' @param ybar A vector of estimated average relative treatment effects based on a previous NMA. The jth entry is the effect of treatment j.
#' @param s A vector of the assumed standard deviation of each treatment. The jth entry is the standard deviation of treatment j.
#' @param mu The hyperparameter mu_0, usually specified as the grand mean of the average treatment effects.
#' @param sigma0 The hyperparameter sigma_0, usually a large number as to be minimally informative.
#' @param tau The standard deviation of the Metropolis Hastings proposal distribution.
#' @param nu0 A numeric vector for the initialization of worth parameters, omega, in the MCMC algorithm. Default to \code{NULL}, indicating random initialization.
#' @param num_iters A numeric indicating the total number of outer MCMC iterations (i.e., the number of times the partition is updated in the Gibbs sampler).
#' @param nu_reps A numeric indicating the number of times each worth parameter is drawn per update of the parameter partition. There will be a total of \code{num_iters}x\code{nu_reps} samples from the posterior.
#'
#' @return A list of 4 elements: \code{omega}, a (\code{num_iters}x\code{nu_reps})x\code{J} matrix of approximate posterior draws of the object-specific worth parameters, omega; \code{nu} a (\code{num_iters}x\code{nu_reps})x\code{J} matrix of the unique parameter values corresponding to the jth partition cluster in posterior draw i, \code{g} a a (\code{num_iters}x\code{nu_reps})x\code{J} matrix indicating the cluster membership of object j in posterior draw i, and \code{K} a vector of the number of non-empty partition clusters in each posterior draw.
#'
#' @examples
#' fit_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu=0.5, sigma0=5, tau=0.5,
#'           nu0 = NULL, num_iters = 5000, nu_reps = 2)
#' @export
fit_RCMVN <-  function(ybar, s, mu, sigma0, tau, nu0 = NULL, num_iters = 5000, nu_reps = 2){
  J <- length(ybar)
  if(length(s)!=J){stop("ybar and s must be of the same length")}
  nu_samples <- matrix(NA, nrow = num_iters * nu_reps, ncol = J)
  g_samples <- matrix(NA, nrow = 0, ncol = J)
  K_samples <- c()
  curr <- 1
  if (is.null(nu0)) {
    nu <- rnorm(J,mean=mu,sd=sigma0)
    g <- 1:J
    K <- J
  }else {
    g <- as.numeric(factor(nu0))
    K <- length(unique(g))
    nu <- as.numeric(levels(factor(nu0)))
  }
  for (iter in 1:num_iters) {
    ## Update partitions, g (and correspondingly, K)
    g_curr <- sample_partition_normal(ybar=ybar, J=J, nu=nu, g=g, K=K, mu=mu,
                                      sigma0=sigma0, s=s, tau=tau)
    g <- g_curr$g
    K <- g_curr$K

    ## Update nu
    nu_curr <- matrix(data = NA, nrow = nu_reps, ncol = K)
    for(k in 1:K){
      C_k <- which(g==k)
      precision <- sum(c(1/sigma0^2,1/s[C_k]^2))
      post_mean <- sum(c(mu/sigma0^2,ybar[C_k]/s[C_k]^2))/precision
      post_sd <- sqrt(1/precision)
      nu_curr[,k] <- rnorm(nu_reps,mean=post_mean,sd=post_sd)
    }
    nu <- nu_curr[nu_reps,]

    ## Save results; update counter
    g_samples <- rbind(g_samples, matrix(rep(g, nu_reps), ncol = J, byrow = T))
    K_samples <- c(K_samples, rep(K, nu_reps))
    nu_samples[curr:(curr + nu_reps - 1), 1:K] <- nu_curr
    curr <- curr + nu_reps
  }
  omega_samples <- matrix(unlist(lapply(1:nrow(nu_samples),
                                        function(iter) {
                                          nu_samples[iter, g_samples[iter, ]]
                                        })), nrow = nrow(nu_samples), byrow = T)
  return(list(omega = omega_samples,
              nu = nu_samples,
              g = g_samples,
              K = K_samples))
}
