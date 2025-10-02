#' Create posterior forest plots for the omega parameters in RaCE NMA models
#'
#' This function creates posterior forest plots for the omega parameters in the form of a ggplot.
#'
#' @import ggplot2
#'
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#' @param names A vector of intervention names (optional)
#' @param limits A numeric indicating the desired credible level to be displayed, as a proportion. Defaults to 0.95.
#'
#' @return A ggplot of posterior forest plots for the omega parameter.
#'
#' @examples
#' mcmc <- mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' create_forestplot(mcmc,names=paste0("Treatment ",1:4))
#' @export
create_forestplot <- function(mcmc,names=NULL,limits=0.95){

  posterior_omega <- mcmc[,grep("omega",names(mcmc))]
  posterior_omega_mean <- apply(posterior_omega,2,mean)
  posterior_omega_quantiles <- apply(posterior_omega,2,function(omega){
    quantile(omega,c((1-limits)/2,(1+limits)/2))
  })
  posterior_omega_order <- order(posterior_omega_mean)

  posterior_omega <- as.data.frame(t(rbind(posterior_omega_mean,posterior_omega_quantiles)))
  names(posterior_omega) <- c("mean","lower_CI","upper_CI")


  if(is.null(names)){
    posterior_omega$name <- factor(paste0("omega",1:nrow(posterior_omega)),
                                   levels=paste0("omega",posterior_omega_order))
  }else{
    posterior_omega$name <- factor(names,levels=names[posterior_omega_order])
  }
  g <- ggplot(posterior_omega,aes(y=name,x=mean,xmin=lower_CI,xmax=upper_CI))+
    geom_point()+geom_linerange()+theme_minimal()+
    scale_y_discrete(limits=rev)+
    theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
    labs(y=element_blank(),x="Relative Treatment Effect")
  g
}
