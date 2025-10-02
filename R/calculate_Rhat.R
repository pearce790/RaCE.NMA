#' Calculate Gelman Diagnostics for Fitted RaCE NMA Models
#'
#' This function applies MCMC outputs from the mcmc_RCMVN package to the gelman.diag function in the coda package.
#'
#' @import coda
#'
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#' @param names A vector of intervention names (optional)
#' @param confidence The \code{confidence} parameter from the \code{gelman.diag} function in the \code{coda} package. Defaults to 0.95.
#' @param multivariate The \code{multivariate} parameter from the \code{gelman.diag} function in the \code{coda} package. Defaults to FALSE.
#'
#' @return Gelman diagnostics for the inputted MCMC chains, in the format of the output of the \code{gelman.diag} function in the \code{coda} package.
#'
#' @examples
#' mcmc <- mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' calculate_Rhat(mcmc,names=paste0("Treatment ",1:4))
#' @export
calculate_Rhat <- function(mcmc,names=NULL,confidence=0.95,multivariate=FALSE){
  start_iter <- min(mcmc$iteration)
  end_iter <- max(mcmc$iteration)
  thin <- diff(mcmc$iteration[1:2])
  mcmc_list <- list()
  for(chain_iter in 1:max(as.numeric(mcmc$chain))){
    tmp <- mcmc[as.numeric(mcmc$chain)==chain_iter,c(grep("omega",names(mcmc)))]
    if(!is.null(names)){names(tmp) <- names}
    mcmc_list[[chain_iter]] <- as.mcmc(tmp,start=start_iter,end=end_iter,thin=thin)
  }
  mcmc_list <- as.mcmc.list(mcmc_list)
  return(gelman.diag(mcmc_list,confidence=confidence,multivariate=multivariate))
}
