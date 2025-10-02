#' Create trace plots for the K parameter in RaCE NMA models
#'
#' This function creates trace plots for the K parameter in the form of a ggplot.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @import utils
#'
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#'
#' @return A ggplot of trace plots for the K parameter.
#'
#' @examples
#' mcmc <- mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' createtrace_K(mcmc)
#' @export
createtrace_K <- function(mcmc){
  g <- ggplot(data=mcmc,aes(x=iteration,y=K,color=chain))+
    geom_line()+theme_minimal()+
    theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
    labs(x="Iteration",y="Posterior of K",color="Chain")
  return(g)
}
