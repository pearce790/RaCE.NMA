#' Create a cumulative ranking plot in RaCE NMA models
#'
#' This function creates a cumulative ranking plot in the form of a ggplot.
#'
#' @import ggplot2
#'
#' @param data A NxJ matrix of data, where N is the number of observations and J the number of treatments. This feature is designed to display results from a standard NMA study.
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#' @param names A vector of intervention names (optional)
#'
#' @return A ggplot of a cumulative ranking plot.
#'
#' @examples
#' mcmc <- mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' create_cumulativeranking(mcmc=mcmc,names=paste0("Treatment ",1:4))
#' @export
create_cumulativeranking <- function(data=NULL, mcmc=NULL, names=NULL){
  if(!is.null(data)){
    J <- ncol(data)
    if(is.null(names)){
      names <- paste(1:J)
    }
    data_ranks <- as.data.frame(t(apply(data,1,function(mu){rank(mu,ties.method="min")})))
    names(data_ranks) <- names
    data_meanorder <- order(apply(data_ranks,2,mean))
    data_rank_cumulative_probability <- melt(apply(data_ranks,2,function(ranks){
      unlist(lapply(1:J,function(j){mean(ranks<=j)}))}))

    g <- ggplot(data_rank_cumulative_probability,aes(x=Var1,y=value,group=Var2,color=Var2))+
      geom_line()+theme_minimal()+
      scale_x_continuous(breaks=1:J)+
      theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
      labs(x="Rank",y="Cumulative Probability",color=NULL)
  }
  if(!is.null(mcmc)){
    J <- (ncol(mcmc)-3)/3
    posterior_mu <- mcmc[,grep("mu",names(mcmc))]
    posterior_ranks <- t(apply(posterior_mu,1,function(mu){rank(mu,ties.method = "min")}))
    posterior_meanorder <- order(apply(posterior_ranks,2,mean))
    posterior_rank_cumulative_probability <- melt(apply(posterior_ranks,2,function(ranks){
      unlist(lapply(1:J,function(j){mean(ranks<=j)}))}))
    if(is.null(names)){
      posterior_rank_cumulative_probability$Var2 <- factor(posterior_rank_cumulative_probability$Var2,
                                                           levels=paste0("mu",posterior_meanorder),
                                                           labels=paste0("Treatment ",posterior_meanorder))
    }else{
      posterior_rank_cumulative_probability$Var2 <- factor(posterior_rank_cumulative_probability$Var2,
                                                           levels=paste0("mu",posterior_meanorder),
                                                           labels=names[posterior_meanorder])
    }
    g <- ggplot(posterior_rank_cumulative_probability,aes(x=Var1,y=value,group=Var2,color=Var2))+
      geom_line()+theme_minimal()+
      scale_x_continuous(breaks=1:J)+
      theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
      labs(x="Rank",y="Cumulative Probability",color=NULL)
  }
  return(g)
}
