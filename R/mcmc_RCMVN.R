#' Fit a Bayesian Rank-Clustered Estimation model for Network Meta-Analysis (RaCE-NMA) using multiple MCMC chains
#'
#' This function fits a Bayesian RaCE-NMA model to data from a previous network meta-analysis. The function has input parameters to permit drawing multiple MCMC chains, as well as thinning and burn-in.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @import utils
#'
#' @param ybar A vector of estimated average relative treatment effects based on a previous NMA. The jth entry is the effect of treatment j.
#' @param s A vector of the assumed standard deviation of each treatment. The jth entry is the standard deviation of treatment j.
#' @param mu The hyperparameter mu_0, usually specified as the grand mean of the average treatment effects.
#' @param sigma0 The hyperparameter sigma_0, usually a large number as to be minimally informative.
#' @param tau The standard deviation of the Metropolis Hastings proposal distribution.
#' @param nu0 A numeric vector for the initialization of worth parameters, omega, in the MCMC algorithm. Default to \code{NULL}, indicating random initialization.
#' @param num_iters A numeric indicating the total number of outer MCMC iterations (i.e., the number of times the partition is updated in the Gibbs sampler).
#' @param nu_reps A numeric indicating the number of times each worth parameter is drawn per update of the parameter partition. There will be a total of \code{num_iters}x\code{nu_reps} samples from the posterior.
#' @param chains A numeric indicating the total number of independent MCMC chains to be run.
#' @param burn_prop A numeric between 0 and 1 indicating the proportion of MCMC samples in each chain to be removed as burn-in.
#' @param thin A numeric indicating that only every \code{thin}-th sample should be retained, to save computational memory.
#' @param seed A numeric indicating the random seed that should be set before running the first MCMC chain.
#' @param suppressPrint A boolean indicating if the function should not print progress updates as the MCMC chains run.
#'
#' @return A (\code{chains}x\code{num_iters}/\code{thin})x(3J+3) matrix of posterior draws, one row per posterior sample of omega, nu, and g, with additional columns indicating the MCMC chain index, iteration index, and number of non-empty partition clusters K of each posterior sample.
#'
#' @examples
#' mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' @export
mcmc_RCMVN <- function(ybar, s, mu, sigma0, tau, nu0 = NULL, num_iters = 100, nu_reps = 5, chains=4, burn_prop=0.5, thin=1, seed=NULL,suppressPrint=FALSE){
  J <- length(ybar)
  if (!is.null(seed)) {set.seed(seed)}
  counter <- 1
  mcmc <- replicate(n = chains,{
    if(!suppressPrint){print(paste0("Estimating chain ", counter, " of ", chains,"."))}
    counter <<- counter + 1
    res <- fit_RCMVN(ybar=ybar, s=s, mu = mu, sigma0 = sigma0, tau = tau, nu0 = nu0, num_iters = num_iters, nu_reps = nu_reps)
    nreps <- nrow(res$omega)
    keep_reps <- seq(ceiling(burn_prop * nreps) + 1, nreps,by = thin)
    tmp <- as.data.frame(cbind(res$omega[keep_reps, , drop=FALSE], res$nu[keep_reps, ,drop=FALSE],
                               res$g[keep_reps, ,drop=FALSE],res$K[keep_reps], keep_reps))
    names(tmp) <- c(paste0("omega", 1:J), paste0("nu", 1:ncol(res$nu)), paste0("G", 1:J), "K", "iteration")
    return(tmp)
  }, simplify = FALSE)
  mcmc <- do.call(rbind, mcmc)
  mcmc$chain <- factor(rep(1:chains, each = nrow(mcmc)/chains))
  mcmc <- mcmc %>% dplyr::select(chain, iteration, K, dplyr::everything())
  rm(counter)
  return(mcmc)
}