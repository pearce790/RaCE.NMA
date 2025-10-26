source("source.R")
library(cowplot)
library(gridExtra)
library(kableExtra)
library(reshape2)

### Simulation Study 1 ####
s <- 0.1
g1a<-ggplot()+scale_x_continuous(limits=c(-3*s,1+3*s))+
  stat_function(fun = dnorm, args = list(mean = 0, sd =s))+
  stat_function(fun = dnorm, args = list(mean = 1, sd =s))+
  theme_bw()+labs(x=" ",y="Density",title=expression(hat(sigma)~"=0.1"))+
  theme(panel.grid.minor = element_blank(),panel.grid.major.y = element_blank())
s <- 0.3
g1b<-ggplot()+scale_x_continuous(limits=c(-3*s,1+3*s))+
  stat_function(fun = dnorm, args = list(mean = 0, sd =s))+
  stat_function(fun = dnorm, args = list(mean = 1, sd =s))+
  theme_bw()+labs(x="Posterior Intervention Effects",y=element_blank(),title=expression(hat(sigma)~"=0.3"))+
  theme(panel.grid.minor = element_blank(),panel.grid.major.y = element_blank())
s <- 0.5
g1c<-ggplot()+scale_x_continuous(limits=c(-3*s,1+3*s))+
  stat_function(fun = dnorm, args = list(mean = 0, sd =s))+
  stat_function(fun = dnorm, args = list(mean = 1, sd =s))+
  theme_bw()+labs(x=" ",y=element_blank(),title=expression(hat(sigma)~"=0.5"))+
  theme(panel.grid.minor = element_blank(),panel.grid.major.y = element_blank())
g1 <- grid.arrange(g1a,g1b,g1c,nrow=1)
ggsave("Simulation1_separation.pdf",g1,width=8,height=2.5)


set.seed(1)
results <- matrix(NA,nrow=0,ncol=6)
for(iter in 1:20){
  print(iter)
  for(J in c(6,12,18)){
    for(K in c(J/3,2*J/3,J)){
      for(s in c(0.1,0.3,0.5)){
        #if(all(iter==4,J==6,K==2,s==0.1)){stop()}
        if(K==J){
          ybar <- 1:K
        }else{
          ybar <- sample(1:K,J,replace=T)
          while(length(unique(ybar))<K){ybar <- sample(1:K,J,replace=T)}
        }
        
        mcmc <- mcmc_RCMVN(ybar=ybar,s=rep(s,J),mu=mean(ybar),sigma0=max(1,4*sd(ybar)),nu0=ybar,
                           tau=1,num_iters=5000,nu_reps=2,chains=4,burn_prop=0.5,thin=3,suppressPrint = TRUE)
        equal <- posterior_equal <- matrix(NA,nrow=J,ncol=J)
        for(i in 1:(J-1)){for(j in (i+1):J){equal[i,j] <- ifelse(ybar[i]==ybar[j],1,0)}}
        for(i in 1:(J-1)){for(j in (i+1):J){posterior_equal[i,j] <- mean(mcmc[,paste0("G",i)] == mcmc[,paste0("G",j)])}}
        
        results <- rbind(results,c(iter,J,K,s,mean(posterior_equal[equal==1],na.rm=T),mean(posterior_equal[equal==0],na.rm=T)))
      }
    }
  }
}
results <- as.data.frame(results)
names(results) <- c("Iteration","J","K","s","Prob_Clustered","Prob_Distinct")
results$s<-factor(results$s,levels=c(.1,.3,.5),labels=c(expression(hat(sigma)~"=0.1"),expression(hat(sigma)~"=0.3"),
                                                        expression(hat(sigma)~"=0.5")))
results$J <- factor(results$J,levels=c(6,12,18),labels=c(expression(J~"= 6"),expression(J~"= 12"),expression(J~"= 18")))
results_melt <- melt(results,id.vars=1:4)
g2<-ggplot(results_melt,aes(x=factor(K),y=value,color=factor(variable,labels=c("Rank-Clustered","Distinct"))))+
  facet_grid(s~J,scales="free_x",labeller=label_parsed)+
  geom_boxplot(outlier.alpha=0.5,position="identity")+theme_bw()+
  scale_color_manual(values=c("skyblue","darkblue"))+
  labs(x="Number of Rank-Clusters, K",y="Posterior Rank-Clustering Probability",
       color=element_blank())+
  theme(legend.position = "bottom",panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
g2
ggsave("Simulation1_results.pdf",g2,width=8,height=4)




### Simulation Study 2 ####
set.seed(2)

J <- 4
ybar <- c(-1,0,1,0)
s <- c(.1,.1,.1,1)
data <- t(matrix(rnorm(10000*4,mean=ybar,sd=s),nrow=4))
g3 <- ggplot(data.frame(name=1:J,mean=ybar,lower_CI=ybar-1.96*s,upper_CI=ybar+1.96*s),
       aes(y=factor(name),x=mean,xmin=lower_CI,xmax=upper_CI))+
  geom_point()+geom_linerange()+scale_y_discrete(limits=rev)+
  theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
  labs(y="Treatment",x="Posterior")
ggsave("Simulation2_data.pdf",g3,width=5,height=2)

mcmc <- mcmc_RCMVN(ybar=ybar,s=s,num_iters=50000,nu_reps=2,chains=4,burn_prop=0.5,thin=1)
calculate_Rhat(mcmc)

## table 1
posterior_equal <- matrix(NA,nrow=J,ncol=J)
for(i in 1:(J-1)){for(j in (i+1):J){posterior_equal[i,j] <- mean(mcmc[,paste0("G",i)] == mcmc[,paste0("G",j)])}}
round(posterior_equal,2)
# [,1] [,2] [,3] [,4]
# [1,]   NA    0    0 0.18
# [2,]   NA   NA    0 0.40
# [3,]   NA   NA   NA 0.19
# [4,]   NA   NA   NA   NA

post_ranks_traditional <- apply(t(apply(data,1,function(omega){rank(omega)})),2,
                                function(rank){unlist(lapply(1:J,function(j){mean(rank==j)}))})
post_ranks <- apply(t(apply(mcmc[,paste0("omega",1:J)],1,function(omega){rank(omega,ties.method="min")})),2,
                    function(rank){unlist(lapply(1:J,function(j){mean(rank==j)}))})
g4a<-ggplot(melt(post_ranks_traditional),aes(x=Var1,y=factor(Var2),fill=value))+
  geom_tile()+scale_y_discrete(limits=rev)+
  scale_x_continuous(breaks=1:4,limits=c(.5,4.5))+
  scale_fill_gradient(low="white",high="black",limits=c(0,1))+
  labs(x="Rank",y="Treatment",fill="Probability")+theme_minimal()+
  theme(panel.grid = element_blank(),legend.position = "bottom")+
  geom_text(aes(x=Var1,y=factor(Var2),label=round(value,2)),
            color=ifelse(melt(post_ranks_traditional)$value>0.4,"white","black"))
g4b<-ggplot(melt(post_ranks),aes(x=Var1,y=factor(Var2,levels=paste0("omega",1:J),labels=paste0(1:J)),fill=value))+
  geom_tile()+scale_y_discrete(limits=rev)+
  scale_x_continuous(breaks=1:4,limits=c(.5,4.5))+
  scale_fill_gradient(low="white",high="black",limits=c(0,1))+
  labs(x="Rank",y="Treatment",fill="Probability")+theme_minimal()+
  theme(panel.grid = element_blank(),legend.position = "bottom")+
  geom_text(aes(x=Var1,y=factor(Var2,levels=paste0("omega",1:J),labels=paste0(1:J)),label=round(value,2)),
            color=ifelse(melt(post_ranks_traditional)$value>0.4,"white","black"))
g4<-plot_grid(plot_grid(g4a+theme(legend.position = "none"), g4b+theme(legend.position = "none"), 
                        labels = c('A', 'B'), label_size = 12),
              get_plot_component(g4a, 'guide-box-bottom', return_all = TRUE),nrow=2,rel_heights = c(.9,.1))
save_plot("Simulation2_results.pdf",g4,base_width=8,base_height=4)



#### Case Study ####
load("WangPosteriors.RData")
wang_posterior <- mcmc.df[,4:13]
ybar <- c(0,apply(wang_posterior,2,mean))
cov <- as.matrix(cbind(c(min(apply(wang_posterior,2,var))/10,rep(0,10)),rbind(0,cov(wang_posterior))))
treatments <- c("R-CHOP","R-CHOP-R","R-Benda","R-Benda-R","R-Benda-R4",
              "R-CVP","R-CVP-R","R-2(Len)","G-CVP-G","G-CHOP-G","G-Benda-G")
rm(mcmc.df)

forestplot_data <- data.frame(
  name=treatments, mean=ybar,
  lower_CI=c(-1.96*sqrt(min(apply(wang_posterior,2,var))/10),apply(wang_posterior,2,function(x){quantile(x,0.025)})),
  upper_CI=c(1.96*sqrt(min(apply(wang_posterior,2,var))/10),apply(wang_posterior,2,function(x){quantile(x,0.975)}))
)
g5<-ggplot(forestplot_data[-1,],aes(y=factor(name,levels=name[order(mean)]),x=mean,xmin=lower_CI,xmax=upper_CI))+
  geom_point()+geom_linerange()+scale_y_discrete(limits=rev)+
  theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank())+
  labs(y="Treatment",x="Posterior")
ggsave("CaseStudy_WangForestPlot.pdf",g5,width=7,height=3)




mcmc <- mcmc_RCMVN(ybar=ybar,cov=cov,num_iters=50000,nu_reps=5,chains=4,seed=1)
calculate_Rhat(mcmc,names=treatments)

g6 <- createtrace_K(mcmc)+scale_x_continuous(labels=paste0(c(125, 150, 175,200,225,250),"k"))
ggsave("CaseStudy_TraceK.pdf",g6,width=12,height=6)

g7 <- createtrace_omega(mcmc,names=treatments)+scale_x_continuous(labels=paste0(c(125, 150, 175,200,225,250),"k"))
ggsave("CaseStudy_TraceOmega.pdf",g7,width=12,height=6)

g8 <- create_forestplot(mcmc,names=treatments)
ggsave("CaseStudy_RCForestPlot.pdf",g8,width=7,height=3)



mcmc_ranks <- as.data.frame(t(apply(mcmc[,4:14],1,function(omega){rank(omega,ties.method="min")})))
mcmc_ranks_probs <- as.data.frame(apply(mcmc_ranks,2,function(ranks){unlist(lapply(1:11,function(j){mean(ranks==j)}))}))
mcmc_ranks_cumprobs <- as.data.frame(apply(mcmc_ranks_probs,2,cumsum))
names(mcmc_ranks) <- names(mcmc_ranks_probs) <- names(mcmc_ranks_cumprobs) <- treatments


wang_ranks <- as.data.frame(t(apply(cbind(0,wang_posterior),1,function(omega){rank(omega,ties.method="min")})))
wang_ranks_probs <- as.data.frame(apply(wang_ranks,2,function(ranks){unlist(lapply(1:11,function(j){mean(ranks==j)}))}))
wang_ranks_cumprobs <- as.data.frame(apply(wang_ranks_probs,2,cumsum))
names(wang_ranks) <- names(wang_ranks_probs) <- names(wang_ranks_cumprobs) <- treatments


g9a <- ggplot(melt(cbind(rank=1:11,wang_ranks_probs),id.vars = 1),
              aes(x=rank,y=factor(variable,levels=treatments[order(ybar)]),fill=value))+
  geom_tile()+scale_x_continuous(breaks=1:11,limits=c(0.5,11.5))+
  scale_y_discrete(limits=rev)+
  scale_fill_gradient(low="white",high="black",limits=c(0,1))+
  labs(x="Rank",y=element_blank(),fill="Probability ")+
  theme_minimal()+theme(panel.grid = element_blank(),legend.position = "bottom")+
  geom_text(data = melt(cbind(rank=1:11,wang_ranks_probs),id.vars = 1) %>% filter(rank==1),
            aes(x=rank,y=variable,label=round(value,2)),
            color=c(rep("black",10),"white"))
g9b <- create_clustermatrix(mcmc,treatments,1)+theme(legend.position="bottom")
g9<-plot_grid(plot_grid(g9a+theme(legend.position = "none"), g9b+theme(legend.position = "none"), 
                        labels = c('A', 'B'), label_size = 12),
              get_plot_component(g9a, 'guide-box-bottom', return_all = TRUE),nrow=2,rel_heights = c(.9,.1))
save_plot("CaseStudy_ClusterMatrix.pdf",g9,base_width=11,base_height=5)


g10a <- ggplot(melt(cbind(rank=1:11,wang_ranks_cumprobs),id.vars = 1),
       aes(x=rank,y=value,color=factor(variable,levels=treatments[order(ybar)])))+
  geom_line()+labs(x="Rank",y="Cumulative Probability",color=element_blank())+
  theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
                   legend.position = "bottom")+guides(color = guide_legend(nrow = 2))
g10b <- create_cumulativeranking(mcmc,treatments)
g10<-plot_grid(plot_grid(g10a+theme(legend.position = "none"), g10b+theme(legend.position = "none"), 
                        labels = c('A', 'B'), label_size = 12),
              get_plot_component(g10a, 'guide-box-bottom', return_all = TRUE),nrow=2,rel_heights = c(.85,.15))
save_plot("Case_Study_CumRankCurve.pdf",g10,base_width=10,base_height=4)

table2 <- data.frame(
  SUCRA_Wang = apply(wang_ranks_cumprobs[1:10,],2,mean),
  SUCRA_RaCE = apply(mcmc_ranks_cumprobs[1:10,],2,mean),
  MedBetter_Wang = apply(wang_ranks,2,function(ranks){
    values <- quantile(ranks-1,c(0.5,0.25,0.75))
    paste0(values[1]," (",values[2],", ",values[3],")")
  }),
  MedBetter_RaCE = apply(mcmc_ranks,2,function(ranks){
    values <- quantile(ranks-1,c(0.5,0.25,0.75))
    paste0(values[1]," (",values[2],", ",values[3],")")
  })
)
kable(table2[order(ybar),],format = "latex",digits=3)
# \begin{tabular}{l|r|r|l|l}
# & SUCRA\_Wang & SUCRA\_RaCE & MedBetter\_Wang & MedBetter\_RaCE\\
# \hline
# G-Benda-G & 0.967 & 0.992 & 0 (0, 1) & 0 (0, 0)\\
# R-Benda-R4 & 0.884 & 0.943 & 1 (0, 1) & 0 (0, 1)\\
# R-Benda-R & 0.812 & 0.925 & 2 (2, 2) & 1 (0, 1)\\
# G-CHOP-G & 0.663 & 0.647 & 3 (3, 4) & 3 (3, 4)\\
# R-CHOP-R & 0.514 & 0.657 & 5 (4, 5) & 3 (3, 4)\\
# R-2(Len) & 0.501 & 0.604 & 5 (4, 6) & 3 (3, 5)\\
# R-Benda & 0.482 & 0.657 & 5 (4, 6) & 3 (3, 4)\\
# G-CVP-G & 0.278 & 0.276 & 7 (6, 8) & 7 (7, 8)\\
# R-CHOP & 0.186 & 0.279 & 8 (7, 9) & 7 (7, 8)\\
# R-CVP-R & 0.159 & 0.296 & 8 (8, 9) & 7 (7, 7)\\
# R-CVP & 0.054 & 0.264 & 10 (9, 10) & 7 (7, 8)\\
# \end{tabular}


## Sensitivity

cov_mod <- cov
cov_mod[1,1] <- min(apply(wang_posterior,2,var))/100
mcmc_mod <- mcmc_RCMVN(ybar=ybar,cov=cov_mod,num_iters=50000,nu_reps=5,chains=4,seed=1)
calculate_Rhat(mcmc_mod,names=treatments)

g11a <- create_clustermatrix(mcmc_mod,treatments,1)+theme(legend.position="bottom")
g11<-plot_grid(plot_grid(g11a+theme(legend.position = "none"), g9b+theme(legend.position = "none"), 
                        labels = c('A', 'B'), label_size = 12),
              get_plot_component(g11a, 'guide-box-bottom', return_all = TRUE),nrow=2,rel_heights = c(.9,.1))
save_plot("CaseStudy_ModVar_ClusterMatrix.pdf",g11,base_width=11,base_height=5)

mcmc_mod_ranks <- as.data.frame(t(apply(mcmc_mod[,4:14],1,function(omega){rank(omega,ties.method="min")})))
mcmc_mod_ranks_probs <- as.data.frame(apply(mcmc_mod_ranks,2,function(ranks){unlist(lapply(1:11,function(j){mean(ranks==j)}))}))
mcmc_mod_ranks_cumprobs <- as.data.frame(apply(mcmc_mod_ranks_probs,2,cumsum))
names(mcmc_mod_ranks) <- names(mcmc_mod_ranks_probs) <- names(mcmc_mod_ranks_cumprobs) <- treatments

kable(data.frame(
  SUCRA_RaCE_Mod = apply(mcmc_mod_ranks_cumprobs[1:10,],2,mean),
  MedBetter_RaCE_Mod = apply(mcmc_mod_ranks,2,function(ranks){
    values <- quantile(ranks-1,c(0.5,0.25,0.75))
    paste0(values[1]," (",values[2],", ",values[3],")")})
)[order(ybar),],format = "latex",digits=3)
# \begin{tabular}{l|r|l}
# \hline
# & SUCRA\_RaCE\_Mod & MedBetter\_RaCE\_Mod\\
# \hline
# G-Benda-G & 0.992 & 0 (0, 0)\\
# R-Benda-R4 & 0.941 & 0 (0, 1)\\
# R-Benda-R & 0.922 & 1 (0, 1)\\
# G-CHOP-G & 0.652 & 3 (3, 4)\\
# R-CHOP-R & 0.646 & 3 (3, 4)\\
# R-2(Len) & 0.632 & 3 (3, 4)\\
# R-Benda & 0.648 & 3 (3, 4)\\
# G-CVP-G & 0.292 & 7 (7, 8)\\
# R-CHOP & 0.264 & 7 (7, 8)\\
# R-CVP-R & 0.290 & 7 (7, 8)\\
# R-CVP & 0.256 & 7 (7, 8)\\
# \end{tabular}

