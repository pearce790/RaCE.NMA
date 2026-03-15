#' Fit a Bayesian Rank-Clustered Estimation model for Network Meta-Analysis (RaCE-NMA). Recommended for internal use only; use mcmc_RCMVN instead.
#'
#' This function fits a Bayesian RaCE-NMA model to data from a previous network meta-analysis.
#'
#' @import invgamma
#' @import mvtnorm
#' @import reshape2
#' @importFrom stats dnorm
#' @importFrom stats na.exclude
#' @importFrom stats quantile
#' @importFrom stats rbinom
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats var
#'
#'
#' @param posterior A matrix of posterior draws of relative intervention effects based on a previous NMA. The (i,j) is the ith draw of the effect of intervention j.
#' @param mu_hat A vector of estimated average relative intervention effects based on a previous NMA. The jth entry is the effect of intervention j. Ignored if \code{posterior} is supplied.
#' @param cov A variance covariance matrix of relative intervention effects based on a previous NMA. The (i,j) entry is the covariane between intervention i and j's effects. Ignored if \code{posterior} is supplied.
#' @param s A vector of the estimated standard deviations of each intervention. The jth entry is the standard deviation of intervention j. Ignored if \code{posterior} is supplied.
#' @param mu0 The hyperparameter mu0. If \code{NULL}, set to mean(mu_hat).
#' @param sigma0 The hyperparameter sigma_0. If \code{NULL}, set to sqrt(10*var(mu_hat)) which aims to be minimally informative.
#' @param tau The standard deviation of the Metropolis Hastings proposal distribution. If \code{NULL}, set to min(|mu_hat_i-mu_hat_j|).
#' @param nu0 A numeric vector for the initialization of worth parameters, mu, in the MCMC algorithm. Default to \code{NULL}, indicating random initialization.
#' @param iter A numeric indicating the total number of outer MCMC iterations (i.e., the number of times the partition is updated in the Gibbs sampler).
#' @param nu_iter A numeric indicating the number of times each worth parameter is drawn per update of the parameter partition. There will be a total of \code{iter}x\code{nu_iter} samples from the posterior.
#' @return A list of 4 elements: \code{mu}, a (\code{iter}x\code{nu_iter})x\code{J} matrix of approximate posterior draws of the intervention-specific worth parameters, mu; \code{nu} a (\code{iter}x\code{nu_iter})x\code{J} matrix of the unique parameter values corresponding to the jth partition cluster in posterior draw i, \code{g} a a (\code{iter}x\code{nu_iter})x\code{J} matrix indicating the cluster membership of object j in posterior draw i, and \code{K} a vector of the number of non-empty partition clusters in each posterior draw.
#'
#' @examples
#' fit_raceNMA(mu_hat=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,
#'           nu0 = NULL, iter = 5000, nu_iter = 2)
#' @export
fit_raceNMA <-  function(posterior=NULL, mu_hat=NULL, cov=NULL, s=NULL, mu0=NULL, sigma0=NULL, tau, nu0, iter, nu_iter){

  ## Checks to ensure appropriate data provided
  if(!is.null(posterior)){
    corr <- TRUE
    if(!is.null(mu_hat)){
      rm(mu_hat)
      warning("Since `posterior` provided, ignoring mu_hat")
    }
    if(!is.null(s)){
      rm(s)
      warning("Since `posterior` provided, ignoring s")
    }
    if(!is.null(cov)){
      rm(cov)
      warning("Since `posterior` provided, ignoring cov")
    }
    mu_hat <- apply(posterior,2,mean)
    cov <- cov(posterior)
  }else{
    if(is.null(mu_hat)){stop("If posterior samples are not provided, must supply mu_hat")}
    if(!is.null(cov)){
      corr <- TRUE
      if(!is.null(s)){
        rm(s)
        warning("Since `posterior` provided, ignoring s")}
    }else{
      corr <- FALSE
      if(is.null(s)){stop("If posterior samples and cov are not provided, must supply s")}
      if(!is.null(cov)){
        rm(cov)
        warning("Since `posterior` provided, ignoring cov")
      }
    }
  }
  J <- length(mu_hat)
  if(corr){
    if(nrow(cov)!=J){stop("length(mu_hat) must equal nrow(cov)")}
  }else{
    if(length(s)!=J){stop("mu_hat and s must be of the same length")}
  }

  if(is.null(mu0)){mu0 <- mean(mu_hat)}
  if(is.null(sigma0)){sigma0 <- sqrt(10*var(mu_hat))}

  nu_samples <- matrix(NA, nrow = iter * nu_iter, ncol = J)
  g_samples <- matrix(NA, nrow = 0, ncol = J)
  K_samples <- c()
  curr <- 1
  if (is.null(tau)){
    tau <- 0.5*max(abs(diff(sort(mu_hat))))
  }
  if (is.null(nu0)) {
    if(corr){
      nu <- rnorm(J,mean=mu_hat,sd=tau)
    }else{
      nu <- rnorm(J,mean=mu_hat,sd=tau)
    }
    g <- 1:J
    K <- J
  }else {
    g <- as.numeric(factor(nu0))
    K <- length(unique(g))
    nu <- as.numeric(levels(factor(nu0)))
  }

  for (it in 1:iter) {
    ## Update partitions, g (and correspondingly, K)
    if(corr){
      g_curr <- sample_partition_correlation(mu_hat=mu_hat, J=J, nu=nu, g=g, K=K, mu0=mu0,
                                             sigma0=sigma0, cov=cov, tau=tau)
      g <- g_curr$g
      nu<- g_curr$nu
      K <- g_curr$K
    }else{
      g_curr <- sample_partition_independence(mu_hat=mu_hat, J=J, nu=nu, g=g, K=K, mu0=mu0,
                                              sigma0=sigma0, s=s, tau=tau)
      g <- g_curr$g
      K <- g_curr$K
    }

    ## Update nu
    nu_curr <- matrix(data = NA, nrow = nu_iter, ncol = K)
    if(corr){
      for(nu_it in 1:nu_iter){
        for(k in 1:K){
          nu_prop <- nu
          nu_prop[k] <- nu[k]+rnorm(1,mean=0,sd=tau)
          logtransition_prob <- (dmvnorm(mu_hat,mean=nu_prop[g],sigma=cov,log=T)+dnorm(nu_prop[k],mu0,sigma0,log=T))-
            (dmvnorm(mu_hat,mean=nu[g],sigma=cov,log=T)+dnorm(nu[k],mu0,sigma0,log=T))
          if(logtransition_prob > log(runif(1))){
            nu[k] <- nu_prop[k]
          }
        }
        nu_curr[nu_it,] <- nu
      }
    }else{
      for(k in 1:K){
        C_k <- which(g==k)
        precision <- sum(c(1/sigma0^2,1/s[C_k]^2))
        post_mean <- sum(c(mu0/sigma0^2,mu_hat[C_k]/s[C_k]^2))/precision
        post_sd <- sqrt(1/precision)
        nu_curr[,k] <- rnorm(nu_iter,mean=post_mean,sd=post_sd)
      }
      nu <- nu_curr[nu_iter,]
    }

    ## Save results; update counter
    g_samples <- rbind(g_samples, matrix(rep(g, nu_iter), ncol = J, byrow = T))
    K_samples <- c(K_samples, rep(K, nu_iter))
    nu_samples[curr:(curr + nu_iter - 1), 1:K] <- nu_curr
    curr <- curr + nu_iter
  }
  mu_samples <- matrix(sapply(1:nrow(nu_samples),function(it) {
    nu_samples[it, g_samples[it, ]]
  }), nrow = nrow(nu_samples), byrow = T)
  return(list(mu = mu_samples,
              nu = nu_samples,
              g = g_samples,
              K = K_samples))
}
