#' Create posterior forest plots for relative treatment effect estimates from a standard or RaCE NMA study
#'
#' This function creates posterior forest plots for relative treatment effect estimates from a standard or RaCE NMA study.
#'
#' @import ggplot2
#'
#' @param data A NxJ matrix of data to display as a forest plot, where N is the number of observations and J the number of treatments. This feature is designed for use to display a forest plot of results from a standard NMA study.
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#' @param names A vector of intervention names (optional)
#' @param limits A numeric indicating the desired credible level to be displayed, as a proportion. Defaults to 0.95.
#' @param order_by_average A boolean indicating if plot should order treatments by their average treatment effect.
#'
#' @return A ggplot of posterior forest plots for the mu parameter.
#'
#' @examples
#' data("toy_data")
#' create_forestplot(data=toy_data,names=paste0("Treatment ",1:4))
#' @export
create_forestplot <- function(data=NULL,mcmc=NULL,names=NULL,limits=0.95,order_by_average=TRUE){
  if(!is.null(data)){
    data_mean <- apply(data,2,mean)
    data_quantiles <- apply(data,2,function(mu){quantile(mu,c((1-limits)/2,(1+limits)/2))})
    data_summary <- as.data.frame(t(rbind(data_mean,data_quantiles)))
    names(data_summary) <- c("mean","lower_CI","upper_CI")

    if(order_by_average){
      mean_order <- order(data_mean)
    }else{
      mean_order <- 1:length(data_mean)
    }

    if(is.null(names)){
      data_summary$name <- factor(paste0(mean_order))
      g <- ggplot(data_summary,aes(y=name,x=mean,xmin=lower_CI,xmax=upper_CI))+
        geom_point()+geom_linerange()+theme_minimal()+
        scale_y_discrete(limits=rev)+
        theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
        labs(y="Treatment",x="Posterior")
    }else{
      data_summary$name <- factor(names,levels=names[mean_order])
      g <- ggplot(data_summary,aes(y=name,x=mean,xmin=lower_CI,xmax=upper_CI))+
        geom_point()+geom_linerange()+theme_minimal()+
        scale_y_discrete(limits=rev)+
        theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
        labs(y=NULL,x="Posterior")
    }
  }
  if(!is.null(mcmc)){
    posterior_mu <- mcmc[,grep("mu",names(mcmc))]
    posterior_mu_mean <- apply(posterior_mu,2,mean)
    posterior_mu_quantiles <- apply(posterior_mu,2,function(mu){
      quantile(mu,c((1-limits)/2,(1+limits)/2))
    })
    if(order_by_average){
      posterior_mu_order <- order(posterior_mu_mean)
    }else{
      posterior_mu_order <- 1:length(posterior_mu_mean)
    }


    posterior_mu <- as.data.frame(t(rbind(posterior_mu_mean,posterior_mu_quantiles)))
    names(posterior_mu) <- c("mean","lower_CI","upper_CI")


    if(is.null(names)){
      posterior_mu$name <- factor(paste0("mu",1:nrow(posterior_mu)),
                                  levels=paste0("mu",posterior_mu_order))
    }else{
      posterior_mu$name <- factor(names,levels=names[posterior_mu_order])
    }
    g <- ggplot(posterior_mu,aes(y=name,x=mean,xmin=lower_CI,xmax=upper_CI))+
      geom_point()+geom_linerange()+theme_minimal()+
      scale_y_discrete(limits=rev)+
      theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
      labs(y=NULL,x="Relative Treatment Effect")
  }
  return(g)
}
