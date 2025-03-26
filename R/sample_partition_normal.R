#' Partition sampling in Bayesian RaCE-NMA models (internal use only)
#'
#' This function implements a reversible jump MCMC procedure for updating the parameter partition in Bayesian Rank-Clustered Estimation for Network Meta-Analysis models. For internal use only.
#'
#' @param mu_hat A vector of estimated average relative intervention effects based on a previous NMA. The jth entry is the effect of intervention j.
#' @param J A numeric indicating the total number of interventions being compared.
#' @param nu A vector indicating current values for nu in the Gibbs sampler.
#' @param g A vector indicating current values for g in the Gibbs sampler.
#' @param K A vector indicating current values for K in the Gibbs sampler.
#' @param mu0 The hyperparameter mu0, usually specified as the grand mean of the average intervention effects.
#' @param sigma0 The hyperparameter sigma_0, usually a large number as to be minimally informative.
#' @param sigma_hat A vector of the estimated standard deviations of each intervention. The jth entry is the standard deviation of intervention j.
#' @param tau The standard deviation of the Metropolis Hastings proposal distribution.
#' @param b_g The probability of "birth"ing a new partition cluster, if possible. Default is 0.5.
#' @param d_g The probability of "death"ing an existing partition cluster, if possible. Default is 0.5.
#'
#' @return A list containing updated values for g, nu, and K.
#'
#' @export
sample_partition_normal <- function (mu_hat, J, nu, g, K, mu0, sigma0, sigma_hat, tau = tau, b_g = 0.5, d_g = 0.5){
  logprior_partition <- log(rep(1,J))
  S_g <- unlist(lapply(1:K, function(k) {sum(g == k)}))
  if (rbinom(1, 1, b_g) == 1) {
    ## Birth!
    if (all(S_g == 1)) {
      logprob_accept <- -Inf
    }else {
      which_split <- which(S_g >= 2)
      if (length(which_split) == 1) {k <- which_split
      }else {k <- sample(which_split, 1)}
      new_clusts <- c(1, rbinom(S_g[k] - 1, 1, prob = 0.5))
      while (sum(new_clusts) == 0 | sum(new_clusts) == S_g[k]) {
        new_clusts <- c(1, rbinom(S_g[k] - 1, 1, prob = 0.5))
      }
      u <- rnorm(1,mean=0,sd=tau)
      nu_new <- nu
      nu_new[k] <- nu[k] - u
      nu_new[K + 1] <- nu[k] + u
      g_new <- g
      g_new[which(g == k)[new_clusts == 1]] <- K + 1
      K_new <- K+1
      if (any((min(nu_new[c(k, K + 1)]) <= nu[-k]) & (nu[-k] <= max(nu_new[c(k, K + 1)])))) {
        logprob_accept <- -Inf
      }else {
        logprob_accept <-
          sum(unlist(lapply(1:J,function(j){dnorm(mu_hat[j], mean = nu_new[g_new[j]], sd = sigma_hat[j],log=T)}))) -
          sum(unlist(lapply(1:J,function(j){dnorm(mu_hat[j], mean = nu[g[j]], sd = sigma_hat[j],log=T)}))) +
          sum(dnorm(nu_new[c(k, K+1)], mu0, sigma0, log=T)) - dnorm(nu[k], mu0, sigma0, log=T) +
          logprior_partition[K + 1] - logprior_partition[K] +
          log(d_g) + log(sum(S_g >= 2)) + log(2^(S_g[k]) - 2) - log(b_g) - log(K + 1 - 1) - log(2) -  dnorm(u,mean=0,sd=tau,log=T) +
          log(2)
      }
    }
    if (log(runif(1)) < logprob_accept) {
      g <- unlist(lapply(1:J, function(j) {
        which(sort(nu_new) == nu_new[g_new][j])
      }))
      nu <- sort(nu_new)
      K <- K_new
      if (any(nu_new[g_new] != nu[g])) {
        print("Something wrong!")
      }
    }
  } else {
    ## Death!
    if (K == 1) {
      logprob_accept <- -Inf
    }else {
      if (K == 2) {
        which_merge <- c(1, 2)
      }else {
        which_merge <- sample(1:(K - 1), 1)
        which_merge <- c(which_merge, which_merge + 1)
      }
      nu_new <- nu
      nu_new[which_merge[1]] <- mean(nu[which_merge])
      nu_new[which_merge[2]] <- NA
      g_new <- g
      g_new[g_new == which_merge[2]] <- which_merge[1]
      u <- (nu[which_merge[2]]-nu[which_merge[1]])/2
      S_gnew <- unlist(lapply(1:max(g_new), function(k) {
        sum(g_new == k)
      }))
      K_new <- K-1
      logprob_accept <-
        sum(unlist(lapply(1:J,function(j){dnorm(mu_hat[j], mean = nu_new[g_new[j]], sd = sigma_hat[j],log=T)}))) -
        sum(unlist(lapply(1:J,function(j){dnorm(mu_hat[j], mean = nu[g[j]], sd = sigma_hat[j],log=T)}))) +
        dnorm(nu_new[which_merge[1]], mu0, sigma0, log=T) -
        sum(dnorm(nu[which_merge], mu0, sigma0, log=T)) +
        logprior_partition[K - 1] - logprior_partition[K] +
        log(b_g) + log(K-1) + log(2) + dnorm(u,mean=0,sd=tau,log=T) - log(d_g) - log(sum(S_gnew >= 2)) - log(2^(S_gnew[which_merge[1]]) - 2) +
        log(1/2)
    }
    if (log(runif(1)) < logprob_accept) {
      g_new <- as.numeric(as.factor(g_new))
      nu_new <- c(na.exclude(nu_new))
      g <- unlist(lapply(1:J, function(j) {
        which(sort(nu_new) == nu_new[g_new][j])
      }))
      nu <- sort(nu_new)
      K <- K_new
      if (any(nu_new[g_new] != nu[g])) {
        print("Something wrong!")
      }
    }
  }
  g
  return(list(g = g, nu = nu, K = K))
}
