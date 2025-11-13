#' Create trace plots for the mu parameter in RaCE NMA models
#'
#' This function creates trace plots for the mu parameter in the form of a ggplot.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @import utils
#'
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#' @param names A vector of intervention names (optional)
#'
#' @return A ggplot of trace plots for the mu parameter.
#'
#' @examples
#' mcmc <- mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' createtrace_mu(mcmc,names=paste0("Treatment ",1:4))
#' @export
createtrace_mu <- function(mcmc,names=NULL){
  J <- (ncol(mcmc)-3)/3
  mu_posterior <- mcmc[,c("chain","iteration",paste0("mu",1:J))]
  if(is.null(names)){
    names(mu_posterior)[3:ncol(mu_posterior)] <- paste0("Treatment ",1:J)
  }else{
    names(mu_posterior)[3:ncol(mu_posterior)] <- names
  }
  mu_posterior <- melt(mu_posterior,id.vars = 1:2)
  g <- ggplot(data=mu_posterior,aes(x=iteration,y=value,color=chain))+
    geom_line()+facet_wrap(~variable,scales = "free_y")+theme_minimal()+
    theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
    labs(x="Iteration",y=expression("Posterior of"~mu),color="Chain")
  return(g)
}
