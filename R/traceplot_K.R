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
#' mcmc <- mcmc_raceNMA(mu_hat=c(0,0,1,1), s=c(.1,.1,.1,.1), seed=1)
#' traceplot_K(mcmc)
#' @export
traceplot_K <- function(mcmc){
  g <- ggplot(data=mcmc,aes(x=iteration,y=K,color=chain))+
    geom_line()+theme_minimal()+
    theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
    labs(x="Iteration",y="Posterior of K",color="Chain")
  return(g)
}
