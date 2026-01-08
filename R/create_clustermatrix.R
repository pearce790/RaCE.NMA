#' Create a posterior clustering matrix for the interventions based on RaCE NMA models
#'
#' This function creates a posterior clustering matrix for the interventions in the form of a ggplot.
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#'
#' @param data A NxJ matrix of data, where N is the number of observations and J the number of treatments. This feature is designed to display results from a standard NMA study.
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#' @param names A vector of intervention names (optional)
#' @param label_ranks A vector containing rank levels for which posterior rank probabilities should be displayed within the clustering matrix. Defaults to \code{NULL}, indicating no probabilities are displayed as text.
#'
#' @return A ggplot of a posterior clustering matrix for the interventions.
#'
#' @examples
#' mcmc <- mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' create_clustermatrix(mcmc=mcmc,names=paste0("Treatment ",1:4),label_ranks=c(1,3))
#' @export
create_clustermatrix <- function(data=NULL,mcmc=NULL,names=NULL,label_ranks=NULL){
  if(!is.null(data)){
    J <- ncol(data)
    if(is.null(names)){
      names <- paste(1:J)
    }
    data_ranks <- as.data.frame(t(apply(data,1,function(mu){rank(mu,ties.method="min")})))
    data_ranks_probs <- as.data.frame(apply(data_ranks,2,function(ranks){sapply(1:J,function(j){mean(ranks==j)})}))
    names(data_ranks_probs) <- names

    g <- ggplot(melt(cbind(rank=1:J,data_ranks_probs),id.vars = 1),
           aes(x=rank,y=factor(variable,levels=names[order(apply(data,2,mean))]),fill=value))+
      geom_tile()+scale_x_continuous(breaks=1:J,limits=c(0.5,J+0.5))+
      scale_y_discrete(limits=rev)+
      scale_fill_gradient(low="white",high="black",limits=c(0,1))+
      labs(x="Rank",y=NULL,fill="Probability ")+
      theme_minimal()+theme(panel.grid = element_blank(),legend.position = "right")
    if(!is.null(label_ranks)){
      rank_data <- melt(as.matrix(data_ranks_probs)) %>% filter(Var1 %in% label_ranks)
      g<- g + geom_text(data=rank_data,
                        aes(x=Var1,y=Var2,label=round(value,2)),
                        color=ifelse(rank_data$value>0.5,"white","black"))
    }
  }
  if(!is.null(mcmc)){
    J <- (ncol(mcmc)-3)/3
    posterior_mu <- mcmc[,grep("mu",names(mcmc))]
    posterior_ranks <- t(apply(posterior_mu,1,function(mu){rank(mu,ties.method = "min")}))
    posterior_meanorder <- order(apply(posterior_ranks,2,mean))
    posterior_rank_probability <- melt(apply(posterior_ranks,2,function(ranks){
      unlist(lapply(1:J,function(j){mean(ranks==j)}))
    }))
    if(is.null(names)){
      posterior_rank_probability$Var2 <- factor(posterior_rank_probability$Var2,
                                                levels=paste0("mu",posterior_meanorder),
                                                labels=paste0("Treatment ",posterior_meanorder))
    }else{
      posterior_rank_probability$Var2 <- factor(posterior_rank_probability$Var2,
                                                levels=paste0("mu",posterior_meanorder),
                                                labels=names[posterior_meanorder])
    }
    g<-ggplot(posterior_rank_probability,aes(x=Var1,y=Var2,fill=value))+
      geom_tile()+scale_x_continuous(breaks=1:J,limits=c(0.5,J+0.5))+
      scale_y_discrete(limits=rev)+
      scale_fill_gradient(low="white",high="black",limits=c(0,1))+
      labs(x="Rank",y=element_blank(),fill="Probability ")+
      theme_minimal()+theme(panel.grid = element_blank(),legend.position = "right")
    if(!is.null(label_ranks)){
      rank_data <- posterior_rank_probability %>% filter(Var1 %in% label_ranks)
      g<- g + geom_text(data=rank_data,
                        aes(x=Var1,y=Var2,label=round(value,2)),
                        color=ifelse(rank_data$value>0.5,"white","black"))
    }
  }
  return(g)
}
