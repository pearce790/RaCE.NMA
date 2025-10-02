#' Create posterior violin plots for the omega parameters in RaCE NMA models
#'
#' This function creates posterior violin plots for the omega parameters in the form of a ggplot.
#'
#' @import ggplot2
#'
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#' @param names A vector of intervention names (optional)
#'
#' @return A ggplot of posterior violin plots for the omega parameters.
#'
#' @examples
#' mcmc <- mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' create_violinplot(mcmc,names=paste0("Treatment ",1:4))
#' @export
create_violinplot <- function(mcmc,names=NULL){

  posterior_omega <- mcmc[,grep("omega",names(mcmc))]
  posterior_omega_order <- order(apply(posterior_omega,2,mean))
  suppressMessages(posterior_omega <- melt(posterior_omega))
  if(is.null(names)){
    posterior_omega$variable <- factor(posterior_omega$variable,
                                       levels=paste0("omega",posterior_omega_order),
                                       labels=paste0("Treatment",posterior_omega_order))
  }else{
    posterior_omega$variable <- factor(posterior_omega$variable,
                                       levels=paste0("omega",posterior_omega_order),
                                       labels=names[posterior_omega_order])
  }
  g <- ggplot(posterior_omega,aes(y=variable,x=value))+
    geom_violin()+theme_minimal()+
    scale_y_discrete(limits=rev)+
    theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
    labs(y=element_blank(),x="Relative Treatment Effect")
  g
}
