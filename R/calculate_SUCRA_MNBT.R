#' Calculate SUCRA and MNBT based on MCMC draws in NMA or RaCE-NMA models
#'
#' This function creates a results table of SUCRA and Median Number of Better Treatments (MNBT) based on
#' MCMC draws from a standard NMA model or a fitted RaCE-NMA model.
#'
#' @import dplyr
#'
#' @param data A NxJ matrix of data, where N is the number of observations and J the number of treatments. This feature is designed to display results from a standard NMA study.
#' @param mcmc MCMC draws from the RaCE NMA model, in the form of the model output of the \code{mcmc_RCMVN} function.
#' @param confidence The \code{confidence} parameter from the \code{gelman.diag} function in the \code{coda} package. Defaults to 0.95.
#' @param names A vector of intervention names (optional)
#' @return A table containing SUCRA and MNBT values for each treatment, ordered by descending SUCRA values.
#'
#' @examples
#' mcmc <- mcmc_RCMVN(ybar=c(0,0,1,1), s=c(.1,.1,.1,.1), mu0=0.5, sigma0=5, tau=0.5,chains=2,seed=1)
#' calculate_SUCRA_MNBT(mcmc=mcmc)
#' @export
calculate_SUCRA_MNBT <- function(data=NULL,mcmc=NULL,confidence=0.50,names=NULL){
  if(!is.null(data)){
    J <- ncol(data)
    if(is.null(names)){
      names <- paste(1:J)
    }
    data_ranks <- as.data.frame(t(apply(data,1,function(mu){rank(mu,ties.method="min")})))
    data_ranks_probs <- apply(data_ranks,2,function(ranks){sapply(1:J,function(j){mean(ranks==j)})})
    data_ranks_probs_cumulative <- apply(data_ranks_probs,2,cumsum)


    SUCRA = apply(data_ranks_probs_cumulative[1:(J-1),],2,mean)
    MNBT = apply(data_ranks,2,function(ranks){
      values <- quantile(ranks-1,c(0.5,(1-confidence)/2,1-(1-confidence)/2))
      paste0(values[1]," (",values[2],", ",values[3],")")
    })
    result <- data.frame(Treatment=names,SUCRA=SUCRA,MNBT=MNBT) %>% arrange(desc(SUCRA),row.names=FALSE)
    rownames(result) <- NULL
    names(result)[3] <- paste0("MNBT (",confidence*100,"% CI)")
  }
  if(!is.null(mcmc)){
    mcmc_mu <- mcmc[,grep("mu",names(mcmc))]
    J <- ncol(mcmc_mu)
    if(is.null(names)){
      names <- paste(1:J)
    }
    mcmc_ranks <- as.data.frame(t(apply(mcmc_mu,1,function(mu){rank(mu,ties.method="min")})))
    mcmc_ranks_probs <- apply(mcmc_ranks,2,function(ranks){sapply(1:J,function(j){mean(ranks==j)})})
    mcmc_ranks_probs_cumulative <- apply(mcmc_ranks_probs,2,cumsum)


    SUCRA = apply(mcmc_ranks_probs_cumulative[1:(J-1),],2,mean)
    MNBT = apply(mcmc_ranks,2,function(ranks){
      values <- quantile(ranks-1,c(0.5,(1-confidence)/2,1-(1-confidence)/2))
      paste0(values[1]," (",values[2],", ",values[3],")")
    })
    result <- data.frame(Treatment=names,SUCRA=SUCRA,MNBT=MNBT) %>% arrange(desc(SUCRA),row.names=FALSE)
    rownames(result) <- NULL
    names(result)[3] <- paste0("MNBT (",confidence*100,"% CI)")
  }
  return(result)
}
