library(tidyverse)
library(invgamma)
library(mvtnorm)
library(reshape2)
library(coda)

## Provided Data (all optional; some needed)
##    posterior: a matrix of posterior samples, one column per treatment (if provided, ybar, s, cov ignored)
##    ybar: a vector of mean treatment-specific effects
##    cov: a covariance matrix of treatment-specific effects (if provided, s ignored)
##    s: a vector of standard deviations of treatment-specific effects
##    names: names of each treatment

## Parameters:
##    mu0: the prior mean of treatment-specific mean effects
##        (Empirical Bayes: Set as mean(ybar) )
##    sigma0: the prior standard deviation of treatment-specific mean effects
##        (Empirical Bayes: Set as sqrt(10*Variance(ybar)) )

## Algorithm Parameters:
##    tau: standard deviation of transition probability
##        (Default: min(|ybar_i-ybar_j|) )
##    nu0: starting values of nu
##        (Default: ybar)

sample_partition_independence <- function (ybar, J, nu, g, K, mu0, sigma0, s, tau = tau, b_g = 0.5, d_g = 0.5){
  logprior_partition <- log(rep(1,J))
  S_g <- unlist(lapply(1:K, function(k) {sum(g == k)}))
  if (rbinom(1, 1, b_g) == 1) {
    ## Birth!
    if (all(S_g == 1)) {
      logprob_accept <- -Inf
    }else {
      which_split <- which(S_g >= 2)
      if (length(which_split) == 1) {k <- which_split
      }else {k <- sample(which_split, 1)}
      new_clusts <- c(1, rbinom(S_g[k] - 1, 1, prob = 0.5))
      while (sum(new_clusts) == 0 | sum(new_clusts) == S_g[k]) {
        new_clusts <- c(1, rbinom(S_g[k] - 1, 1, prob = 0.5))
      }
      u <- rnorm(1,mean=0,sd=tau)
      nu_new <- nu
      nu_new[k] <- nu[k] - u
      nu_new[K + 1] <- nu[k] + u
      g_new <- g
      g_new[which(g == k)[new_clusts == 1]] <- K + 1
      K_new <- K+1
      if (any((min(nu_new[c(k, K + 1)]) <= nu[-k]) & (nu[-k] <= max(nu_new[c(k, K + 1)])))) {
        logprob_accept <- -Inf
      }else {
        logprob_accept <-
          sum(unlist(lapply(1:J,function(j){dnorm(ybar[j], mean = nu_new[g_new[j]], sd = s[j],log=T)}))) -
          sum(unlist(lapply(1:J,function(j){dnorm(ybar[j], mean = nu[g[j]], sd = s[j],log=T)}))) +
          sum(dnorm(nu_new[c(k, K+1)], mu0, sigma0, log=T)) - dnorm(nu[k], mu0, sigma0, log=T) +
          logprior_partition[K + 1] - logprior_partition[K] +
          log(d_g) + log(sum(S_g >= 2)) + log(2^(S_g[k]) - 2) - log(b_g) - log(K + 1 - 1) - log(2) -  dnorm(u,mean=0,sd=tau,log=T) +
          log(2)
      }
    }
    if (log(runif(1)) < logprob_accept) {
      g <- unlist(lapply(1:J, function(j) {
        which(sort(nu_new) == nu_new[g_new][j])
      }))
      nu <- sort(nu_new)
      K <- K_new
      if (any(nu_new[g_new] != nu[g])) {
        print("Something wrong!")
      }
    }
  } else {
    ## Death!
    if (K == 1) {
      logprob_accept <- -Inf
    }else {
      if (K == 2) {
        which_merge <- c(1, 2)
      }else {
        which_merge <- sample(1:(K - 1), 1)
        which_merge <- c(which_merge, which_merge + 1)
      }
      nu_new <- nu
      nu_new[which_merge[1]] <- mean(nu[which_merge])
      nu_new[which_merge[2]] <- NA
      g_new <- g
      g_new[g_new == which_merge[2]] <- which_merge[1]
      u <- (nu[which_merge[2]]-nu[which_merge[1]])/2
      S_gnew <- unlist(lapply(1:max(g_new), function(k) {
        sum(g_new == k)
      }))
      K_new <- K-1
      logprob_accept <-
        sum(unlist(lapply(1:J,function(j){dnorm(ybar[j], mean = nu_new[g_new[j]], sd = s[j],log=T)}))) -
        sum(unlist(lapply(1:J,function(j){dnorm(ybar[j], mean = nu[g[j]], sd = s[j],log=T)}))) +
        dnorm(nu_new[which_merge[1]], mu0, sigma0, log=T) -
        sum(dnorm(nu[which_merge], mu0, sigma0, log=T)) +
        logprior_partition[K - 1] - logprior_partition[K] + 
        log(b_g) + log(K-1) + log(2) + dnorm(u,mean=0,sd=tau,log=T) - log(d_g) - log(sum(S_gnew >= 2)) - log(2^(S_gnew[which_merge[1]]) - 2) +
        log(1/2)
    }
    if (log(runif(1)) < logprob_accept) {
      g_new <- as.numeric(as.factor(g_new))
      nu_new <- c(na.exclude(nu_new))
      g <- unlist(lapply(1:J, function(j) {
        which(sort(nu_new) == nu_new[g_new][j])
      }))
      nu <- sort(nu_new)
      K <- K_new
      if (any(nu_new[g_new] != nu[g])) {
        print("Something wrong!")
      }
    }
  }
  g
  return(list(g = g, nu = nu, K = K))
}
sample_partition_correlation <- function (ybar, J, nu, g, K, mu0, sigma0, cov, tau = tau, b_g = 0.5, d_g = 0.5){
  logprior_partition <- log(rep(1,J))
  S_g <- unlist(lapply(1:K, function(k) {sum(g == k)}))
  if (rbinom(1, 1, b_g) == 1) {
    ## Birth!
    if (all(S_g == 1)) {
      logprob_accept <- -Inf
    }else {
      which_split <- which(S_g >= 2)
      if (length(which_split) == 1) {k <- which_split
      }else {k <- sample(which_split, 1)}
      new_clusts <- c(1, rbinom(S_g[k] - 1, 1, prob = 0.5))
      while (sum(new_clusts) == 0 | sum(new_clusts) == S_g[k]) {
        new_clusts <- c(1, rbinom(S_g[k] - 1, 1, prob = 0.5))
      }
      u <- rnorm(1,mean=0,sd=tau)
      nu_new <- nu
      nu_new[k] <- nu[k] - u
      nu_new[K + 1] <- nu[k] + u
      g_new <- g
      g_new[which(g == k)[new_clusts == 1]] <- K + 1
      K_new <- K+1
      if (any((min(nu_new[c(k, K + 1)]) <= nu[-k]) & (nu[-k] <= max(nu_new[c(k, K + 1)])))) {
        logprob_accept <- -Inf
      }else {
        logprob_accept <-
          dmvnorm(x=ybar,mean=nu_new[g_new],sigma=cov,log=T) - dmvnorm(x=ybar,mean=nu[g],sigma=cov,log=T) +
          sum(dnorm(nu_new[c(k, K+1)], mu0, sigma0, log=T)) - dnorm(nu[k], mu0, sigma0, log=T) +
          logprior_partition[K + 1] - logprior_partition[K] +
          log(d_g) + log(sum(S_g >= 2)) + log(2^(S_g[k]) - 2) - log(b_g) - log(K + 1 - 1) - log(2) -  dnorm(u,mean=0,sd=tau,log=T) +
          log(2)
      }
    }
    if (log(runif(1)) < logprob_accept) {
      g <- unlist(lapply(1:J, function(j) {
        which(sort(nu_new) == nu_new[g_new][j])
      }))
      nu <- sort(nu_new)
      K <- K_new
      if (any(nu_new[g_new] != nu[g])) {
        print("Something wrong!")
      }
    }
  } else {
    ## Death!
    if (K == 1) {
      logprob_accept <- -Inf
    }else {
      if (K == 2) {
        which_merge <- c(1, 2)
      }else {
        which_merge <- sample(1:(K - 1), 1)
        which_merge <- c(which_merge, which_merge + 1)
      }
      nu_new <- nu
      nu_new[which_merge[1]] <- mean(nu[which_merge])
      nu_new[which_merge[2]] <- NA
      g_new <- g
      g_new[g_new == which_merge[2]] <- which_merge[1]
      u <- (nu[which_merge[2]]-nu[which_merge[1]])/2
      S_gnew <- unlist(lapply(1:max(g_new), function(k) {
        sum(g_new == k)
      }))
      K_new <- K-1
      logprob_accept <-
        dmvnorm(x=ybar,mean=nu_new[g_new],sigma=cov,log=T) - dmvnorm(x=ybar,mean=nu[g],sigma=cov,log=T) +
        dnorm(nu_new[which_merge[1]], mu0, sigma0, log=T) - sum(dnorm(nu[which_merge], mu0, sigma0, log=T)) +
        logprior_partition[K - 1] - logprior_partition[K] + 
        log(b_g) + log(K-1) + log(2) + dnorm(u,mean=0,sd=tau,log=T) - log(d_g) - log(sum(S_gnew >= 2)) - log(2^(S_gnew[which_merge[1]]) - 2) +
        log(1/2)
    }
    if (log(runif(1)) < logprob_accept) {
      g_new <- as.numeric(as.factor(g_new))
      nu_new <- c(na.exclude(nu_new))
      g <- unlist(lapply(1:J, function(j) {
        which(sort(nu_new) == nu_new[g_new][j])
      }))
      nu <- sort(nu_new)
      K <- K_new
      if (any(nu_new[g_new] != nu[g])) {
        print("Something wrong!")
      }
    }
  }
  return(list(g = g, nu = nu, K = K))
}
fit_RCMVN <-  function(posterior=NULL, ybar=NULL, cov=NULL, s=NULL, mu0=NULL, sigma0=NULL, tau, nu0, num_iters, nu_reps){
  
  ## Checks to ensure appropriate data provided
  if(!is.null(posterior)){
    corr <- TRUE
    if(!is.null(ybar)){
      rm(ybar)
      warning("Since `posterior` provided, ignoring ybar")
    }
    if(!is.null(s)){
      rm(s)
      warning("Since `posterior` provided, ignoring s")
    }
    if(!is.null(cov)){
      rm(cov)
      warning("Since `posterior` provided, ignoring cov")
    }
    ybar <- apply(posterior,2,mean)
    cov <- cov(posterior)
  }else{
    if(is.null(ybar)){stop("If posterior samples are not provided, must supply ybar")}
    if(!is.null(cov)){
      corr <- TRUE
      if(!is.null(s)){
        rm(s)
        warning("Since `posterior` provided, ignoring s")}
    }else{
      corr <- FALSE
      if(is.null(s)){stop("If posterior samples and cov are not provided, must supply s")}
      if(!is.null(cov)){
        rm(cov)
        warning("Since `posterior` provided, ignoring cov")
      }
    }
  }
  J <- length(ybar)
  if(corr){
    if(nrow(cov)!=J){stop("length(ybar) must equal nrow(cov)")}
  }else{
    if(length(s)!=J){stop("ybar and s must be of the same length")}
  }
  
  if(is.null(mu0)){mu0 <- mean(ybar)}
  if(is.null(sigma0)){sigma0 <- sqrt(10*var(ybar))}
  
  nu_samples <- matrix(NA, nrow = num_iters * nu_reps, ncol = J)
  g_samples <- matrix(NA, nrow = 0, ncol = J)
  K_samples <- c()
  curr <- 1
  if (is.null(tau)){
    tau <- 0.5*max(abs(diff(sort(ybar))))
  }
  if (is.null(nu0)) {
    if(corr){
      nu <- rnorm(J,mean=ybar,sd=tau)
    }else{
      nu <- rnorm(J,mean=ybar,sd=tau)
    }
    g <- 1:J
    K <- J
  }else {
    g <- as.numeric(factor(nu0))
    K <- length(unique(g))
    nu <- as.numeric(levels(factor(nu0)))
  }
  
  for (iter in 1:num_iters) {
    ## Update partitions, g (and correspondingly, K)
    if(corr){
      g_curr <- sample_partition_correlation(ybar=ybar, J=J, nu=nu, g=g, K=K, mu0=mu0,
                                             sigma0=sigma0, cov=cov, tau=tau)
      g <- g_curr$g
      nu<- g_curr$nu
      K <- g_curr$K
    }else{
      g_curr <- sample_partition_independence(ybar=ybar, J=J, nu=nu, g=g, K=K, mu0=mu0,
                                              sigma0=sigma0, s=s, tau=tau)
      g <- g_curr$g
      K <- g_curr$K
    }
    
    ## Update nu
    nu_curr <- matrix(data = NA, nrow = nu_reps, ncol = K)
    if(corr){
      for(nu_rep in 1:nu_reps){
        for(k in 1:K){
          nu_prop <- nu
          nu_prop[k] <- nu[k]+rnorm(1,mean=0,sd=tau)
          logtransition_prob <- (dmvnorm(ybar,mean=nu_prop[g],sigma=cov,log=T)+dnorm(nu_prop[k],mu0,sigma0,log=T))-
            (dmvnorm(ybar,mean=nu[g],sigma=cov,log=T)+dnorm(nu[k],mu0,sigma0,log=T))
          if(logtransition_prob > log(runif(1))){
            nu[k] <- nu_prop[k]
          }
        }
        nu_curr[nu_rep,] <- nu
      }
    }else{
      for(k in 1:K){
        C_k <- which(g==k)
        precision <- sum(c(1/sigma0^2,1/s[C_k]^2))
        post_mean <- sum(c(mu0/sigma0^2,ybar[C_k]/s[C_k]^2))/precision
        post_sd <- sqrt(1/precision)
        nu_curr[,k] <- rnorm(nu_reps,mean=post_mean,sd=post_sd)
      }
      nu <- nu_curr[nu_reps,]
    }
    
    ## Save results; update counter
    g_samples <- rbind(g_samples, matrix(rep(g, nu_reps), ncol = J, byrow = T))
    K_samples <- c(K_samples, rep(K, nu_reps))
    nu_samples[curr:(curr + nu_reps - 1), 1:K] <- nu_curr
    curr <- curr + nu_reps
  }
  omega_samples <- matrix(unlist(lapply(1:nrow(nu_samples),
                                        function(iter) {
                                          nu_samples[iter, g_samples[iter, ]]
                                        })), nrow = nrow(nu_samples), byrow = T)
  return(list(omega = omega_samples,
              nu = nu_samples,
              g = g_samples,
              K = K_samples))
}
mcmc_RCMVN <- function(posterior=NULL, ybar=NULL, cov=NULL, s=NULL, mu0=NULL, sigma0=NULL, 
                       tau = NULL, nu0 = NULL, 
                       num_iters = 5000, nu_reps = 3, chains=2, burn_prop=0.5, thin=1, 
                       seed=NULL,suppressPrint=FALSE){
  
  if(!is.null(posterior)){
    J <- ncol(posterior)
  }else{
    J <- length(ybar)
  }
  
  if (!is.null(seed)) {set.seed(seed)}
  counter <- 1
  mcmc <- replicate(n = chains,{
    if(!suppressPrint){print(paste0("Estimating chain ", counter, " of ", chains,"."))}
    counter <<- counter + 1
    res <- fit_RCMVN(posterior=posterior, ybar=ybar, cov=cov, s=s, mu0 = mu0, sigma0 = sigma0, 
                     tau = tau, nu0 = nu0, num_iters = num_iters, nu_reps = nu_reps)
    nreps <- nrow(res$omega)
    keep_reps <- seq(ceiling(burn_prop * nreps) + 1, nreps,by = thin)
    tmp <- as.data.frame(cbind(res$omega[keep_reps, , drop=FALSE], res$nu[keep_reps, ,drop=FALSE],
                               res$g[keep_reps, ,drop=FALSE],res$K[keep_reps], keep_reps))
    names(tmp) <- c(paste0("omega", 1:J), paste0("nu", 1:ncol(res$nu)), paste0("G", 1:J), "K", "iteration")
    return(tmp)
  }, simplify = FALSE)
  mcmc <- do.call(rbind, mcmc)
  mcmc$chain <- factor(rep(1:chains, each = nrow(mcmc)/chains))
  mcmc <- mcmc %>% dplyr::select(chain, iteration, K, dplyr::everything())
  rm(counter)
  return(mcmc)
}

calculate_Rhat <- function(mcmc,names=NULL,confidence=0.95,multivariate=FALSE){
  start_iter <- min(mcmc$iteration)
  end_iter <- max(mcmc$iteration)
  thin <- diff(mcmc$iteration[1:2])
  mcmc_list <- list()
  for(chain_iter in 1:max(as.numeric(mcmc$chain))){
    tmp <- mcmc[as.numeric(mcmc$chain)==chain_iter,c(grep("omega",names(mcmc)))]
    if(!is.null(names)){names(tmp) <- names}
    mcmc_list[[chain_iter]] <- as.mcmc(tmp,start=start_iter,end=end_iter,thin=thin)
  }
  mcmc_list <- as.mcmc.list(mcmc_list)
  return(gelman.diag(mcmc_list,confidence=confidence,multivariate=multivariate))
}
createtrace_omega <- function(mcmc,names=NULL){
  J <- (ncol(mcmc)-3)/3
  omega_posterior <- mcmc[,c("chain","iteration",paste0("omega",1:J))]
  if(is.null(names)){
    names(omega_posterior)[3:ncol(omega_posterior)] <- paste0("Treatment ",1:J)
  }else{
    names(omega_posterior)[3:ncol(omega_posterior)] <- names
  }
  omega_posterior <- melt(omega_posterior,id.vars = 1:2)
  g <- ggplot(data=omega_posterior,aes(x=iteration,y=value,color=chain))+
    geom_line()+facet_wrap(~variable,scales = "free_y")+theme_minimal()+
    theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
    labs(x="Iteration",y=expression("Posterior of"~omega),color="Chain")
  return(g)
}
createtrace_K <- function(mcmc){
  g <- ggplot(data=mcmc,aes(x=iteration,y=K,color=chain))+
    geom_line()+theme_minimal()+
    theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
    labs(x="Iteration",y="Posterior of K",color="Chain")
  return(g)
}
create_cumulativeranking <- function(mcmc,names=NULL){
  J <- (ncol(mcmc)-3)/3
  posterior_omega <- mcmc[,grep("omega",names(mcmc))]
  posterior_ranks <- t(apply(posterior_omega,1,function(omega){rank(omega,ties.method = "min")}))
  posterior_meanorder <- order(apply(posterior_ranks,2,mean))
  posterior_rank_cumulative_probability <- melt(apply(posterior_ranks,2,function(ranks){
    unlist(lapply(1:J,function(j){mean(ranks<=j)}))}))
  if(is.null(names)){
    posterior_rank_cumulative_probability$Var2 <- factor(posterior_rank_cumulative_probability$Var2,
                                                         levels=paste0("omega",posterior_meanorder),
                                                         labels=paste0("Treatment ",posterior_meanorder))
  }else{
    posterior_rank_cumulative_probability$Var2 <- factor(posterior_rank_cumulative_probability$Var2,
                                                         levels=paste0("omega",posterior_meanorder),
                                                         labels=names[posterior_meanorder])
  }
  g <- ggplot(posterior_rank_cumulative_probability,aes(x=Var1,y=value,group=Var2,color=Var2))+
    geom_line()+theme_minimal()+
    scale_x_continuous(breaks=1:J)+
    theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
    labs(x="Rank",y="Cumulative Probability",color=element_blank())
  
  return(g)
}
create_clustermatrix <- function(mcmc,names=NULL,label_ranks=NULL){
  J <- (ncol(mcmc)-3)/3
  posterior_omega <- mcmc[,grep("omega",names(mcmc))]
  posterior_ranks <- t(apply(posterior_omega,1,function(omega){rank(omega,ties.method = "min")}))
  posterior_meanorder <- order(apply(posterior_ranks,2,mean))
  posterior_rank_probability <- melt(apply(posterior_ranks,2,function(ranks){
    unlist(lapply(1:J,function(j){mean(ranks==j)}))
  }))
  if(is.null(names)){
    posterior_rank_probability$Var2 <- factor(posterior_rank_probability$Var2,
                                              levels=paste0("omega",posterior_meanorder),
                                              labels=paste0("Treatment ",posterior_meanorder))
  }else{
    posterior_rank_probability$Var2 <- factor(posterior_rank_probability$Var2,
                                              levels=paste0("omega",posterior_meanorder),
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
  return(g)
}
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

