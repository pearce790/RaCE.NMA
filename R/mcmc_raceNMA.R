#' Fit a Bayesian Rank-Clustered Estimation model for Network Meta-Analysis (RaCE-NMA) using multiple MCMC chains
#'
#' This function fits a Bayesian RaCE-NMA model to data from a previous network meta-analysis. The function has input parameters to permit drawing multiple MCMC chains, as well as chain thinning and burn-in. Chains are independent and may be drawn in parallel across multiple CPU cores via the \code{cores} argument.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @import utils
#'
#' @param posterior A matrix of posterior draws of relative intervention effects based on a previous NMA. The (i,j) is the ith draw of the effect of intervention j.
#' @param mu_hat A vector of estimated average relative intervention effects based on a previous NMA. The jth entry is the effect of intervention j. Ignored if \code{posterior} is supplied.
#' @param cov A variance covariance matrix of relative intervention effects based on a previous NMA. The (i,j) entry is the covariance between intervention i and j's effects. Ignored if \code{posterior} is supplied.
#' @param s A vector of the estimated standard deviations of each intervention. The jth entry is the standard deviation of intervention j. Ignored if \code{posterior} is supplied.
#' @param mu0 The hyperparameter mu0. If \code{NULL}, set to mean(mu_hat).
#' @param sigma0 The hyperparameter sigma_0. If \code{NULL}, set to sqrt(10*var(mu_hat)) which aims to be minimally informative.
#' @param tau The standard deviation of the Metropolis Hastings proposal distribution. If \code{NULL}, set to min(|mu_hat_i-mu_hat_j|).
#' @param nu0 A numeric vector for the initialization of worth parameters, mu, in the MCMC algorithm. Default to \code{NULL}, indicating random initialization.
#' @param iter A numeric indicating the total number of outer MCMC iterations (i.e., the number of times the partition is updated in the Gibbs sampler).
#' @param nu_iter A numeric indicating the number of times each worth parameter is drawn per update of the parameter partition. There will be a total of \code{iter}x\code{nu_iter} samples from the posterior.
#' @param chains A numeric indicating the total number of independent MCMC chains to be run.
#' @param burn_prop A numeric between 0 and 1 indicating the proportion of MCMC samples in each chain to be removed as burn-in.
#' @param thin A numeric indicating that only every \code{thin}-th sample should be retained, to save computational memory.
#' @param seed A numeric indicating the random seed used to generate the per-chain random number streams.
#' @param cores A numeric indicating how many CPU cores to use when drawing chains. The default, \code{1}, draws chains sequentially. Values above 1 draw chains in parallel; \code{cores} is silently capped at \code{chains}. Results do not depend on \code{cores}.
#' @param verbose A boolean indicating if the function should print progress updates as the MCMC chains run. Default to \code{TRUE}.
#'
#' @details
#' Each chain is assigned its own L'Ecuyer-CMRG random number stream, derived deterministically
#' from \code{seed}. Chain \code{i} therefore produces the same draws whether it was run
#' sequentially or on a parallel worker, so \code{cores} affects only run time.
#'
#' @return A (\code{chains}x\code{iter}/\code{thin})x(3J+3) matrix of posterior draws, one row per posterior sample of mu, nu, and g, with additional columns indicating the MCMC chain index, iteration index, and number of non-empty partition clusters K of each posterior sample.
#'
#' @examples
#' mcmc <- mcmc_raceNMA(mu_hat=c(0,0,1,1), s=c(.1,.1,.1,.1), seed=1)
#' head(mcmc)
#'
#' @export
mcmc_raceNMA <- function(posterior = NULL, mu_hat = NULL, cov = NULL, s = NULL, mu0 = NULL, sigma0 = NULL, tau = NULL, nu0 = NULL,
                         iter = 4000, nu_iter = 5, chains = 2, burn_prop = 0.5, thin = 1, seed = NULL,
                         cores = 1, verbose = TRUE){

  if(!is.null(posterior)){
    J <- ncol(posterior)
  }else{
    J <- length(mu_hat)
  }

  if(!is.numeric(cores) || length(cores) != 1 || is.na(cores) || cores < 1){
    stop("`cores` must be a single number greater than or equal to 1.")
  }
  cores <- min(as.integer(cores), as.integer(chains), parallel::detectCores())

  ## Build one independent RNG stream per chain, derived from `seed`. This makes
  ## output invariant to `cores` and to the order in which workers pick up chains.
  oldkind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(oldkind[1]), add = TRUE)
  if(!is.null(seed)){set.seed(seed, kind = "L'Ecuyer-CMRG")}
  chain_seeds <- vector("list", chains)
  cs <- .Random.seed
  for(i in seq_len(chains)){
    chain_seeds[[i]] <- cs
    cs <- parallel::nextRNGStream(cs)
  }

  ## One chain, given its own RNG stream. Returns a tidy data frame of retained draws.
  run_chain <- function(i){
    assign(".Random.seed", chain_seeds[[i]], envir = globalenv())
    res <- fit_raceNMA(posterior = posterior, mu_hat = mu_hat, cov = cov, s = s, mu0 = mu0,
                       sigma0 = sigma0, tau = tau, nu0 = nu0, iter = iter, nu_iter = nu_iter)
    nreps <- nrow(res$mu)
    keep_reps <- seq(ceiling(burn_prop * nreps) + 1, nreps, by = thin)
    tmp <- as.data.frame(cbind(res$mu[keep_reps, , drop = FALSE], res$nu[keep_reps, , drop = FALSE],
                               res$g[keep_reps, , drop = FALSE], res$K[keep_reps], keep_reps))
    names(tmp) <- c(paste0("mu", 1:J), paste0("nu", 1:ncol(res$nu)), paste0("G", 1:J), "K", "iteration")
    tmp
  }

  if(cores == 1L){
    mcmc <- lapply(seq_len(chains), function(i){
      if(verbose){message("Estimating chain ", i, " of ", chains, ".")}
      run_chain(i)
    })
  }else{
    if(verbose){message("Estimating ", chains, " chains across ", cores, " cores.")}
    cl <- parallel::makeCluster(cores, type = if(.Platform$OS.type == "windows") "PSOCK" else "FORK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if(.Platform$OS.type == "windows"){
      ok <- unlist(parallel::clusterEvalQ(cl, requireNamespace("RaCE.NMA", quietly = TRUE)))
      if(!all(ok)){
        stop("Parallel workers could not load RaCE.NMA. Install the package, or use cores = 1.")
      }
      parallel::clusterExport(cl, varlist = c("chain_seeds", "posterior", "mu_hat", "cov", "s",
                                              "mu0", "sigma0", "tau", "nu0", "iter", "nu_iter",
                                              "burn_prop", "thin", "J"),
                              envir = environment())
      parallel::clusterEvalQ(cl, fit_raceNMA <- RaCE.NMA::fit_raceNMA)
    }
    ## Load balancing matters when chains outnumber cores and chain cost varies.
    mcmc <- parallel::parLapplyLB(cl, seq_len(chains), run_chain)
  }

  chain_rows <- vapply(mcmc, nrow, integer(1))
  mcmc <- do.call(rbind, mcmc)
  mcmc$chain <- factor(rep(seq_len(chains), times = chain_rows))
  mcmc <- mcmc %>% dplyr::select(chain, iteration, K, dplyr::everything())
  return(mcmc)
}
