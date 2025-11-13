#' Create posterior violin plots for the mu parameters in RaCE NMA models
#'
#' This function creates posterior violin plots for the mu parameters in the form of a ggplot.
#'
#' @import ggplot2
#'
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#' @param names A vector of intervention names (optional)
#'
#' @return A ggplot of posterior violin plots for the mu parameters.
#'
#' @examples
#' mcmc <- mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' create_violinplot(mcmc,names=paste0("Treatment ",1:4))
#' @export
create_violinplot <- function(mcmc,names=NULL){

  posterior_mu <- mcmc[,grep("mu",names(mcmc))]
  posterior_mu_order <- order(apply(posterior_mu,2,mean))
  suppressMessages(posterior_mu <- melt(posterior_mu))
  if(is.null(names)){
    posterior_mu$variable <- factor(posterior_mu$variable,
                                       levels=paste0("mu",posterior_mu_order),
                                       labels=paste0("Treatment",posterior_mu_order))
  }else{
    posterior_mu$variable <- factor(posterior_mu$variable,
                                       levels=paste0("mu",posterior_mu_order),
                                       labels=names[posterior_mu_order])
  }
  g <- ggplot(posterior_mu,aes(y=variable,x=value))+
    geom_violin()+theme_minimal()+
    scale_y_discrete(limits=rev)+
    theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
    labs(y=element_blank(),x="Relative Treatment Effect")
  g
}
