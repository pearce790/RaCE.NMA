#' Fit a Bayesian Rank-Clustered Estimation model for Network Meta-Analysis (RaCE-NMA). Recommended for internal use only; use mcmc_RCMVN instead.
#'
#' This function fits a Bayesian RaCE-NMA model to data from a previous network meta-analysis.
#'
#' @param mu_hat A vector of estimated average relative intervention effects based on a previous NMA. The jth entry is the effect of intervention j.
#' @param sigma_hat A vector of the estimated standard deviations of each intervention. The jth entry is the standard deviation of intervention j.
#' @param mu0 The hyperparameter mu0, usually specified as the grand mean of the average intervention effects.
#' @param sigma0 The hyperparameter sigma_0, usually a large number as to be minimally informative.
#' @param tau The standard deviation of the Metropolis Hastings proposal distribution.
#' @param nu0 A numeric vector for the initialization of intervention-specific mean parameters, mu, in the MCMC algorithm. Default to \code{NULL}, indicating random initialization.
#' @param num_iters A numeric indicating the total number of outer MCMC iterations (i.e., the number of times the partition is updated in the Gibbs sampler).
#' @param nu_reps A numeric indicating the number of times each worth parameter is drawn per update of the parameter partition. There will be a total of \code{num_iters}x\code{nu_reps} samples from the posterior.
#'
#' @return A list of 4 elements: \code{mu}, a (\code{num_iters}x\code{nu_reps})x\code{J} matrix of approximate posterior draws of the intervention-specific worth parameters, mu; \code{nu} a (\code{num_iters}x\code{nu_reps})x\code{J} matrix of the unique parameter values corresponding to the jth partition cluster in posterior draw i, \code{g} a a (\code{num_iters}x\code{nu_reps})x\code{J} matrix indicating the cluster membership of object j in posterior draw i, and \code{K} a vector of the number of non-empty partition clusters in each posterior draw.
#'
#' @examples
#' fit_RCMVN(mu_hat=c(0,0,1,1), sigma_hat=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,
#'           nu0 = NULL, num_iters = 5000, nu_reps = 2)
#' @export
fit_RCMVN <-  function(mu_hat, sigma_hat, mu0, sigma0, tau = min(diff(sort(mu_hat)))/2, nu0 = NULL, num_iters = 5000, nu_reps = 2){
  J <- length(mu_hat)
  if(length(sigma_hat)!=J){stop("mu_hat and sigma_hat must be of the same length")}
  nu_samples <- matrix(NA, nrow = num_iters * nu_reps, ncol = J)
  g_samples <- matrix(NA, nrow = 0, ncol = J)
  K_samples <- c()
  curr <- 1
  if (is.null(nu0)) {
    nu <- rnorm(J,mean=mu0,sd=sigma0)
    g <- 1:J
    K <- J
  }else {
    g <- as.numeric(factor(nu0))
    K <- length(unique(g))
    nu <- as.numeric(levels(factor(nu0)))
  }
  for (iter in 1:num_iters) {
    ## Update partitions, g (and correspondingly, K)
    g_curr <- sample_partition_normal(mu_hat=mu_hat, J=J, nu=nu, g=g, K=K, mu0=mu0,
                                      sigma0=sigma0, sigma_hat=sigma_hat, tau=tau)
    g <- g_curr$g
    K <- g_curr$K

    ## Update nu
    nu_curr <- matrix(data = NA, nrow = nu_reps, ncol = K)
    for(k in 1:K){
      C_k <- which(g==k)
      precision <- sum(c(1/sigma0^2,1/sigma_hat[C_k]^2))
      post_mean <- sum(c(mu0/sigma0^2,mu_hat[C_k]/sigma_hat[C_k]^2))/precision
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
  mu_samples <- matrix(
    unlist(lapply(1:nrow(nu_samples),function(iter) {nu_samples[iter, g_samples[iter, ]]})),
    nrow = nrow(nu_samples), byrow = T
  )

  return(list(mu = mu_samples,
              nu = nu_samples,
              g = g_samples,
              K = K_samples))
}
