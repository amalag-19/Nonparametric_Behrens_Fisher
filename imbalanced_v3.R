#########################################################################################################
## Loading the required packages
library(perm)
library(lawstat)
library(foreach)
library(doParallel)

#########################################################################################################
## Writing a function to conduct Bootstrap t test and Wilcoxon test
# pooled sample, assumes equal variance
boot_test<-function(sample_1,sample_2,test_stat_vec,alternative_val,boot_size){
  sample_combined<-c(sample_1,sample_2)
  sample_1_normalized<-sample_1-mean(sample_1)+mean(sample_combined)
  sample_2_normalized<-sample_2-mean(sample_2)+mean(sample_combined)
  t_val_vec<-rep(NA_real_,boot_size)
  wilcox_val_vec<-rep(NA_real_,boot_size)
  BM_val_vec<-rep(NA_real_,boot_size)
  #ptm<-proc.time()
  statistic_mat<-cbind(t_val_vec,wilcox_val_vec,BM_val_vec)
  statistic_mat_output<-t(apply(X = statistic_mat,MARGIN = 1,FUN = function(x){
    sample_1_boot_individual<-sample(sample_1_normalized,replace = T)
    sample_2_boot_individual<-sample(sample_2_normalized,replace = T)
    t_val_individual<-t.test(x = sample_1_boot_individual,y = sample_2_boot_individual,alternative = alternative_val)$statistic
    wilcox_val_individual<-wilcox.test(x = sample_1_boot_individual,y = sample_2_boot_individual,alternative = alternative_val)$statistic
    BM_val_individual<-brunner.munzel.test(x = sample_1_boot_individual,y = sample_2_boot_individual,alternative = alternative_val)$statistic
    ## Getting the pooled bootstrap sample
    sample_boot_pooled_indices<-sample(1:length(sample_combined),replace=TRUE)
    sample_boot_pooled<-sample_combined[sample_boot_pooled_indices]
    sample_1_boot_pooled<-sample_boot_pooled[1:length(sample_1)]
    sample_2_boot_pooled<-sample_boot_pooled[(length(sample_1)+1):(length(sample_1)+length(sample_2))]
    t_val_pooled<-t.test(x = sample_1_boot_pooled,y = sample_2_boot_pooled,alternative = alternative_val)$statistic
    wilcox_val_pooled<-wilcox.test(x = sample_1_boot_pooled,y = sample_2_boot_pooled,alternative = alternative_val)$statistic
    BM_val_pooled<-brunner.munzel.test(x = sample_1_boot_pooled,y = sample_2_boot_pooled,alternative = alternative_val)$statistic
    return(c(t_val_individual,wilcox_val_individual,BM_val_individual,t_val_pooled,wilcox_val_pooled,BM_val_pooled))
  }))
  #print(proc.time()-ptm)
  ## Getting the p values
  if(alternative_val=="two.sided"){
    t_individual_p_val<-(sum(abs(statistic_mat_output[,1])>abs(test_stat_vec[1]))+1)/(boot_size+1)
    t_pooled_p_val<-(sum(abs(statistic_mat_output[,4])>abs(test_stat_vec[1]))+1)/(boot_size+1)
    wilcox_individual_p_val<-(sum(abs(statistic_mat_output[,2])>abs(test_stat_vec[2]))+1)/(boot_size+1)
    wilcox_pooled_p_val<-(sum(abs(statistic_mat_output[,5])>abs(test_stat_vec[2]))+1)/(boot_size+1)
    BM_individual_p_val<-(sum(abs(statistic_mat_output[,3])>abs(test_stat_vec[3]))+1)/(boot_size+1)
    BM_pooled_p_val<-(sum(abs(statistic_mat_output[,6])>abs(test_stat_vec[3]))+1)/(boot_size+1)
  }else if(alternative_val=="less"){
    t_individual_p_val<-(sum(statistic_mat_output[,1]<test_stat_vec[1])+1)/(boot_size+1)
    t_pooled_p_val<-(sum(statistic_mat_output[,4]<test_stat_vec[1])+1)/(boot_size+1)
    wilcox_individual_p_val<-(sum(statistic_mat_output[,2]<test_stat_vec[2])+1)/(boot_size+1)
    wilcox_pooled_p_val<-(sum(statistic_mat_output[,5]<test_stat_vec[2])+1)/(boot_size+1)
    BM_individual_p_val<-(sum(statistic_mat_output[,3]<test_stat_vec[3])+1)/(boot_size+1)
    BM_pooled_p_val<-(sum(statistic_mat_output[,6]<test_stat_vec[3])+1)/(boot_size+1)
  }else if(alternative_val=="greater"){
    t_individual_p_val<-(sum(statistic_mat_output[,1]>test_stat_vec[1])+1)/(boot_size+1)
    t_pooled_p_val<-(sum(statistic_mat_output[,4]>test_stat_vec[1])+1)/(boot_size+1)
    wilcox_individual_p_val<-(sum(statistic_mat_output[,2]>test_stat_vec[2])+1)/(boot_size+1)
    wilcox_pooled_p_val<-(sum(statistic_mat_output[,5]>test_stat_vec[2])+1)/(boot_size+1)
    BM_individual_p_val<-(sum(statistic_mat_output[,3]>test_stat_vec[3])+1)/(boot_size+1)
    BM_pooled_p_val<-(sum(statistic_mat_output[,6]>test_stat_vec[3])+1)/(boot_size+1)
  }
  ## Combining the p values
  boot_p_val_vec<-c(t_individual_p_val,t_pooled_p_val,wilcox_individual_p_val,wilcox_pooled_p_val,BM_individual_p_val,BM_pooled_p_val)
  return(boot_p_val_vec)
}

#########################################################################################################
## Writing the function to generate table for test results
test_results_generator<-function(n_simulations,sample_1_size,sample_2_size,dist_type,sample_1_mean=NA,sample_2_mean=NA,sample_1_var=NA,sample_2_var=NA,sample_1_df=NA,sample_2_df=NA,sample_1_ncp=NA,sample_2_ncp=NA,sample_1_shape=NA, sample_2_shape=NA, sample_1_rate=NA, sample_2_rate=NA, sample_1_mix_p = NA, sample_2_mix_p = NA, sample_1_mixComp1_rate = NA, sample_1_mixComp2_rate = NA, sample_2_mixComp1_rate = NA, sample_2_mixComp2_rate = NA, sample_1_mixComp1_shape = NA, sample_1_mixComp2_shape = NA, sample_2_mixComp1_shape = NA, sample_2_mixComp2_shape = NA, alternative_val,boot_size){
  ## Initializing the p-value dataframes for t and wilcoxon test for simple, permutation and bootstrap cases
  df_t_p_val<-data.frame(matrix(NA_real_,n_simulations,4))
  names(df_t_p_val)<-c("Simple","Permutation","Bootstrap_UnPooled","Bootstrap_Pooled")
  df_wilcox_p_val<-data.frame(matrix(NA_real_,n_simulations,3))
  names(df_wilcox_p_val)<-c("Simple","Bootstrap_UnPooled","Bootstrap_Pooled")
  ## Loop to conduct simulations
  for (j in 1:n_simulations){
    set.seed(j)
    if(dist_type=="normal"){
      sample_1<-rnorm(n = sample_1_size,mean = sample_1_mean,sd = sqrt(sample_1_var))
      sample_2<-rnorm(n = sample_2_size,mean = sample_2_mean,sd = sqrt(sample_2_var))
    }else if(dist_type=="t"){
      sample_1<-rt(n = sample_1_size,df = sample_1_df,ncp = sample_1_ncp)
      sample_2<-rt(n = sample_2_size,df = sample_2_df,ncp = sample_2_ncp)
    } else if(dist_type=="exp"){
      sample_1<-rexp(n = sample_1_size,rate = sample_1_rate)
      sample_2<-rexp(n = sample_2_size,rate = sample_2_rate)
    }else if(dist_type=="gamma"){
      sample_1<-rgamma(n = sample_1_size,rate = sample_1_rate,shape=sample_1_shape)
      sample_2<-rgamma(n = sample_2_size,rate = sample_2_rate,shape=sample_2_shape)
    }else if(dist_type=="mixgamma"){
      sample_1_components <- sample(1:2,prob=c(sample_1_mix_p,1-sample_1_mix_p),size=sample_1_size,replace=TRUE)
      sample_1_rates<-c(sample_1_mixComp1_rate,sample_1_mixComp2_rate)
      sample_1_shapes<-c(sample_1_mixComp1_shape,sample_1_mixComp2_shape)
      sample_1<-rgamma(n = sample_1_size,rate = sample_1_rates[sample_1_components],shape=sample_1_shapes[sample_1_components])
      
      sample_2_components <- sample(1:2,prob=c(sample_2_mix_p,1-sample_2_mix_p),size=sample_2_size,replace=TRUE)
      sample_2_rates<-c(sample_2_mixComp1_rate,sample_2_mixComp2_rate)
      sample_2_shapes<-c(sample_2_mixComp1_shape,sample_2_mixComp2_shape)
      sample_2<-rgamma(n = sample_2_size,rate = sample_2_rates[sample_2_components],shape=sample_2_shapes[sample_2_components])
    }
    ## Simple tests
    t_test_obj<-t.test(x = sample_1,y = sample_2,alternative = alternative_val)
    wilcox_test_obj<-wilcox.test(x = sample_1,y = sample_2,alternative = alternative_val)
    df_t_p_val[j,1]<-t_test_obj$p.value
    df_wilcox_p_val[j,1]<-wilcox_test_obj$p.value
    # Getting the general and rank score vectors and group ID vector for permutation t and wilcoxon tests
    gen_score_vec<-c(sample_1,sample_2)
    group_ID_vec<-factor(rep(c("A", "B"), c(length(sample_1), length(sample_2))))
    ## Permutation tests
    df_t_p_val[j,2]<-permTS(gen_score_vec~group_ID_vec,alternative=alternative_val, method="exact.mc",control=permControl(nmc=10^4-1))$p.value
    ## Bootstrap tests
    boot_test_output<-boot_test(sample_1 = sample_1,sample_2 = sample_2,test_stat = c(t_test_obj$statistic,wilcox_test_obj$statistic),alternative_val=alternative_val,boot_size = boot_size)
    df_t_p_val[j,3]<-boot_test_output[1]
    df_t_p_val[j,4]<-boot_test_output[2]
    df_wilcox_p_val[j,2]<-boot_test_output[3]
    df_wilcox_p_val[j,3]<-boot_test_output[4]
    #print(j)
  }
  power_size_t_vec<-colSums(df_t_p_val<0.05)/n_simulations
  power_size_wilcox_vec<-colSums(df_wilcox_p_val<0.05)/n_simulations
  
  if(dist_type=="normal"){
    if(((sample_1_mean==sample_2_mean)&(sample_1_var==sample_2_var))){
      df_output<-data.frame("Size of Simple Test"=numeric(),"Size of Permutation Test"=numeric(),"Size of Bootstrap Test"=numeric())
    }else{
      df_output<-data.frame("Power of Simple Test"=numeric(),"Power of Permutation Test"=numeric(),"Power of Bootstrap Test"=numeric())
    }
  }else if(dist_type=="t"){
    if(((sample_1_df==sample_2_df)&(sample_1_ncp==sample_2_ncp))){
      df_output<-data.frame("Size of Simple Test"=numeric(),"Size of Permutation Test"=numeric(),"Size of Bootstrap Test"=numeric())
    }else{
      df_output<-data.frame("Power of Simple Test"=numeric(),"Power of Permutation Test"=numeric(),"Power of Bootstrap Test"=numeric())
    }
  }else if(dist_type=="exp"){
    if(sample_1_rate==sample_2_rate){
      df_output<-data.frame("Size of Simple Test"=numeric(),"Size of Permutation Test"=numeric(),"Size of Bootstrap Test"=numeric())
    }else{
      df_output<-data.frame("Power of Simple Test"=numeric(),"Power of Permutation Test"=numeric(),"Power of Bootstrap Test"=numeric())
    }
  }else if(dist_type=="gamma"){
    if(((sample_1_rate==sample_2_rate)&(sample_1_shape==sample_2_shape))){
      df_output<-data.frame("Size of Simple Test"=numeric(),"Size of Permutation Test"=numeric(),"Size of Bootstrap Test"=numeric())
    }else{
      df_output<-data.frame("Power of Simple Test"=numeric(),"Power of Permutation Test"=numeric(),"Power of Bootstrap Test"=numeric())
    }
  }else if(dist_type=="mixgamma"){
    if(((sample_1_mix_p==sample_2_mix_p)&(sample_1_mixComp1_rate==sample_2_mixComp1_rate)&(sample_1_mixComp2_rate==sample_2_mixComp2_rate)&(sample_1_mixComp1_shape==sample_2_mixComp1_shape)&(sample_1_mixComp2_shape==sample_2_mixComp2_shape))){
      df_output<-data.frame("Size of Simple Test"=numeric(),"Size of Permutation Test"=numeric(),"Size of Bootstrap Test"=numeric())
    }else{
      df_output<-data.frame("Power of Simple Test"=numeric(),"Power of Permutation Test"=numeric(),"Power of Bootstrap Test"=numeric())
    }
  }
  df_output[1,]<-power_size_t_vec
  df_output[2,]<-power_size_wilcox_vec
  row.names(df_output)<-c("t Test","Wilcoxon Test")
  return(df_output)
}

#########################################################################################################
#undebug(test_results_generator)
#test_results_generator(n_simulations = 100,dist_type = "normal",sample_1_mean = 0,sample_2_mean = 0.5,sample_1_size = 300,sample_2_size = 1200,boot_size = 2)

df_normal_cases<-data.frame("dist"=rep("normal",50),"var"=rep(c(1,5,10,15,20),each=10),"mean"=rep(seq(0,0.45,by=0.05),5),stringsAsFactors = F)

df_normal_cases$simple_t<-rep(NA_real_,nrow(df_normal_cases))
df_normal_cases$permutation_t<-rep(NA_real_,nrow(df_normal_cases))
df_normal_cases$bootstrap_t<-rep(NA_real_,nrow(df_normal_cases))
df_normal_cases$simple_wilcoxon<-rep(NA_real_,nrow(df_normal_cases))
df_normal_cases$permutation_wilcoxon<-rep(NA_real_,nrow(df_normal_cases))

for(i in 1:nrow(df_cases)){
  df_output<-test_results_generator(n_simulations = 100, sample_1_size = 300,sample_2_size = 1600, dist_type = df_normal_cases$dist[i], sample_1_mean = 0, sample_2_mean = df_normal_cases$mean[i], sample_1_var = 1, sample_2_var = df_normal_cases$var[i], boot_size = 100)
  df_normal_cases$simple_t[i]<-df_output[1,1]
  df_normal_cases$permutation_t[i]<-df_output[1,2]
  df_normal_cases$bootstrap_t[i]<-df_output[1,3]
  df_normal_cases$simple_wilcoxon[i]<-df_output[2,1]
  df_normal_cases$permutation_wilcoxon[i]<-df_output[2,2]
  print(i)
}

#########################################################################################################
#########################################################################################################
## Reading the dataset old and New
df_old<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/Old_NW_PA_csv.csv",stringsAsFactors = F)

str(df_old)

prod(df_old[,"CensorCode"]=="nc")
prod(df_old[,"CensorCode.1"]=="nc")
prod(df_old[,"CensorCode.2"]=="nc")
prod(df_old[,"CensorCode.3"]=="nc")
prod(df_old[,"CensorCode.4"]=="nc")
prod(df_old[,"CensorCode.5"]=="nc")
prod(df_old[,"CensorCode.6"]=="nc")
prod(df_old[,"CensorCode.7"]=="nc")
prod(df_old[,"CensorCode.8"]=="nc")
prod(df_old[,"CensorCode.9"]=="nc")
prod(df_old[,"CensorCode.10"]=="nc")
prod(df_old[,"CensorCode.11"]=="nc")
prod(df_old[,"CensorCode.12"]=="nc")
prod(df_old[,"CensorCode.13"]=="nc")
prod(df_old[,"CensorCode.14"]=="nc")

sum(df_old[,"CensorCode"]=="nc")
sum(df_old[,"CensorCode.1"]=="nc")
sum(df_old[,"CensorCode.2"]=="nc")
sum(df_old[,"CensorCode.3"]=="nc")
sum(df_old[,"CensorCode.4"]=="nc")
sum(df_old[,"CensorCode.5"]=="nc")
sum(df_old[,"CensorCode.6"]=="nc")
sum(df_old[,"CensorCode.7"]=="nc")
sum(df_old[,"CensorCode.8"]=="nc")
sum(df_old[,"CensorCode.9"]=="nc")
sum(df_old[,"CensorCode.10"]=="nc")
sum(df_old[,"CensorCode.11"]=="nc")
sum(df_old[,"CensorCode.12"]=="nc")
sum(df_old[,"CensorCode.13"]=="nc")
sum(df_old[,"CensorCode.14"]=="nc")

## Based on highest number of non-censored values, we first analyze top seven analytes as follows
analytes<-c("Chloride","pH","TDS","Specific_conductance","Magnesium_total","Sulfate_total")

######################################################
## Analysis fpr Chloride
df_old$Chloride_total

hist(df_old$Chloride_total,breaks = 50)

library("stats")    
Cl_vec <- df_old$Chloride_total
Cl_pdf <- density(Cl_vec, from= 0, to=1, bw=0.1)

N<-length(df_old$Chloride_total)
Cl_sampled <- rnorm(N, sample(Cl_vec, size = N, replace = TRUE), Cl_pdf$bw)
Cl_sampled<-Cl_sampled[which(Cl_sampled>0)]

# Histogram of the draws with the distribution superimposed.
hist(Cl_sampled, freq = FALSE)
lines(Cl_pdf)

library(MASS)
lambda<-fitdistr(df_old$Chloride_total, "exponential")$estimate
lambda

gamma_param_estimate<-fitdistr(df_old$Chloride_total, "gamma")$estimate
gamma_param_estimate

par(mfrow=c(1,2))
curve(expr = dexp(x = x,rate=lambda), from = 0, to = 1000, n = 101, add = FALSE, type = "l",ylab="probab",main="Exponential (rate=0.041)")
hist(df_old$Chloride_total,breaks = 50,prob=T,add=T)

curve(expr = dgamma(x = x,shape = gamma_param_estimate[1], rate=gamma_param_estimate[2]), from = 0, to = 1000, n = 101, add = FALSE, type = "l",ylab="probab",main="Gamma (shape = 0.629,rate=0.0255)",ylim=c(0,0.04))
hist(df_old$Chloride_total,breaks = 50,prob=T,add=T)

######################################################
## Analysis fpr Magnesium
df_old$Magnesium_total

par(mfrow=c(1,1))
hist(df_old$Magnesium_total,breaks = 50)

Mg_vec<-df_old$Magnesium_total[which(!is.na(df_old$Magnesium_total))]

library(MASS)
gamma_param_estimate<-fitdistr(Mg_vec, "gamma")$estimate
gamma_param_estimate

curve(expr = dgamma(x = x,shape = gamma_param_estimate[1], rate=gamma_param_estimate[2]), from = 0, to = 200, n = 1000, add = FALSE, type = "l",ylab="probab",main="Gamma (shape = 1.344,rate=0.108)")
hist(Mg_vec,breaks = 50,prob=T,add=T)

#########################################################################################################
#########################################################################################################
#########################################################################################################
## New NW PA Dataset

df_new<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/New_NW_PA_csv.csv",stringsAsFactors = F)

str(df_new)

######################################################
## Analysis fpr Chloride
library(MASS)
lambda<-fitdistr(df_new$Chloride_total, "exponential")$estimate
lambda

gamma_param_estimate<-fitdistr(df_new$Chloride_total, "gamma")$estimate
gamma_param_estimate

par(mfrow=c(1,2))
curve(expr = dexp(x = x,rate=lambda), from = 0, to = 1000, n = 101, add = FALSE, type = "l",ylab="probab",main="Exponential (rate=0.0256)")
hist(df_new$Chloride_total,breaks = 50,prob=T,add=T)

curve(expr = dgamma(x = x,shape = gamma_param_estimate[1], rate=gamma_param_estimate[2]), from = 0, to = 1000, n = 101, add = FALSE, type = "l", ylab="probab",main="Gamma (shape = 0.669,rate=0.017)",ylim=c(0,0.025))
hist(df_new$Chloride_total,breaks = 50,prob=T,add=T)

######################################################
## Analysis fpr Magnesium
df_new$Magnesium_total

par(mfrow=c(1,1))
hist(df_new$Magnesium_total,breaks = 50)

Mg_vec<-df_new$Magnesium_total[which(!is.na(df_new$Magnesium_total))]
hist(Mg_vec,breaks = 30)

library(MASS)
gamma_param_estimate<-fitdistr(Mg_vec, "gamma")$estimate
gamma_param_estimate

require(mixtools)
Mg_gmix_output<-gammamixEM(Mg_vec)
Mg_gmix_output$lambda
Mg_gmix_output$gamma.pars

curve(expr = Mg_gmix_output$lambda[1]*dgamma(x = x,shape = Mg_gmix_output$gamma.pars[1,1], rate=1/(Mg_gmix_output$gamma.pars[2,1]))+Mg_gmix_output$lambda[2]*dgamma(x = x,shape = Mg_gmix_output$gamma.pars[1,2], rate=1/(Mg_gmix_output$gamma.pars[2,2])), from = 0, to = 200, n = 1000, add = FALSE, type = "l",ylab = "probab")
hist(Mg_vec,breaks = 50,prob=T,add=T)

#########################################################################################################
#########################################################################################################
## Simulations for size based on real data densties

## Old NW data, Exponential distribution with rate lambda=0.041
df_old_Cl<-data.frame("sample_size_1"=c(300,100,10,10,10), "sample_size_2"= c(1600,11000,50,200,1000),stringsAsFactors = F)

df_old_Cl$simple_t<-rep(NA_real_,nrow(df_old_Cl))
df_old_Cl$permutation_t<-rep(NA_real_,nrow(df_old_Cl))
df_old_Cl$bootstrap_t<-rep(NA_real_,nrow(df_old_Cl))
df_old_Cl$simple_wilcoxon<-rep(NA_real_,nrow(df_old_Cl))
df_old_Cl$permutation_wilcoxon<-rep(NA_real_,nrow(df_old_Cl))

for(i in 1:nrow(df_old_Cl)){
  df_output<-test_results_generator(n_simulations = 10,sample_1_size = df_old_Cl[i,"sample_size_1"],sample_2_size = df_old_Cl[i,"sample_size_2"],dist_type = "exp", sample_1_rate=0.041, sample_2_rate=0.041,boot_size = 1)
  df_old_Cl$simple_t[i]<-df_output[1,1]
  df_old_Cl$permutation_t[i]<-df_output[1,2]
  df_old_Cl$bootstrap_t[i]<-df_output[1,3]
  df_old_Cl$simple_wilcoxon[i]<-df_output[2,1]
  df_old_Cl$permutation_wilcoxon[i]<-df_output[2,2]
  print(i)
}

#########################################################################################################
## New NW data, Exponential distribution with rate lambda=0.026
df_new_Cl<-data.frame("sample_size_1"=c(300,100,10,10,10), "sample_size_2"= c(1600,11000,50,200,1000),stringsAsFactors = F)

df_new_Cl$simple_t<-rep(NA_real_,nrow(df_new_Cl))
df_new_Cl$permutation_t<-rep(NA_real_,nrow(df_new_Cl))
df_new_Cl$bootstrap_t<-rep(NA_real_,nrow(df_new_Cl))
df_new_Cl$simple_wilcoxon<-rep(NA_real_,nrow(df_new_Cl))
df_new_Cl$permutation_wilcoxon<-rep(NA_real_,nrow(df_new_Cl))

for(i in 1:nrow(df_new_Cl)){
  df_output<-test_results_generator(n_simulations = 1000,sample_1_size = df_new_Cl[i,"sample_size_1"],sample_2_size = df_new_Cl[i,"sample_size_2"],dist_type = "exp", sample_1_rate=0.026, sample_2_rate=0.026,boot_size = 1000)
  df_new_Cl$simple_t[i]<-df_output[1,1]
  df_new_Cl$permutation_t[i]<-df_output[1,2]
  df_new_Cl$bootstrap_t[i]<-df_output[1,3]
  df_new_Cl$simple_wilcoxon[i]<-df_output[2,1]
  df_new_Cl$permutation_wilcoxon[i]<-df_output[2,2]
  print(i)
}

#########################################################################################################
## Simulations for power based on real data densties
df_new_old_Cl<-data.frame("sample_size_1"=c(300,100,10,10,10), "sample_size_2"= c(1600,11000,50,200,1000),stringsAsFactors = F)

df_new_old_Cl$simple_t<-rep(NA_real_,nrow(df_new_old_Cl))
df_new_old_Cl$permutation_t<-rep(NA_real_,nrow(df_new_old_Cl))
df_new_old_Cl$bootstrap_t<-rep(NA_real_,nrow(df_new_old_Cl))
df_new_old_Cl$simple_wilcoxon<-rep(NA_real_,nrow(df_new_old_Cl))
df_new_old_Cl$permutation_wilcoxon<-rep(NA_real_,nrow(df_new_old_Cl))

for(i in 1:nrow(df_new_old_Cl)){
  df_output<-test_results_generator(n_simulations = 1000,sample_1_size = df_new_old_Cl[i,"sample_size_1"],sample_2_size = df_new_old_Cl[i,"sample_size_2"],dist_type = "exp", sample_1_rate=0.026, sample_2_rate=0.041,boot_size = 1000)
  df_new_old_Cl$simple_t[i]<-df_output[1,1]
  df_new_old_Cl$permutation_t[i]<-df_output[1,2]
  df_new_old_Cl$bootstrap_t[i]<-df_output[1,3]
  df_new_old_Cl$simple_wilcoxon[i]<-df_output[2,1]
  df_new_old_Cl$permutation_wilcoxon[i]<-df_output[2,2]
  print(i)
}

#########################################################################################################
## Exporting results to Word
load(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/results_normal_10_1000.RData")

write.table(df_normal_cases[,-1], file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/results_txt/results_normal_10_1000.txt", row.names=FALSE, sep = ",")


load(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/results_t_10_1000.RData")

write.table(df_t_cases[,-1], file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/results_txt/results_t_10_1000.txt", row.names=FALSE, sep = ",")

load(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Cluster_Code/real_data/Cl/results/results_new_old_Cl.RData")

write.table(df_new_old_Cl, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Cluster_Code/real_data/Cl/results/results_txt/results_new_old_Cl.txt", row.names=FALSE, sep = ",")

df_normal_cases
xtable(df_normal_cases[,-1])

#########################################################################################################
#########################################################################################################
## Tests for real data (NW PA)
df_old<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/Old_NW_PA_csv.csv",stringsAsFactors = F)

df_new<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/New_NW_PA_csv.csv",stringsAsFactors = F)

analytes<-c("Chloride_total","pH","TDS","Specific_conductance","Magnesium_total","Sulfate_total")

Cl_vec_old<-df_old$Chloride_total[which(!is.na(df_old$Chloride_total))]
Mg_vec_old<-df_old$Magnesium_total[which(!is.na(df_old$Magnesium_total))]
pH_vec_old<-df_old$pH[which(!is.na(df_old$pH))]
TDS_vec_old<-df_old$TDS[which(!is.na(df_old$TDS))]
SC_vec_old<-df_old$Specific_conductance[which(!is.na(df_old$Specific_conductance))]
Su_vec_old<-df_old$Sulfate_total[which(!is.na(df_old$Sulfate_total))]

Cl_vec_new<-df_new$Chloride_total[which(!is.na(df_new$Chloride_total))]
Mg_vec_new<-df_new$Magnesium_total[which(!is.na(df_new$Magnesium_total))]
pH_vec_new<-df_new$pH[which(!is.na(df_new$pH))]
TDS_vec_new<-df_new$TDS[which(!is.na(df_new$TDS))]
SC_vec_new<-df_new$Specific_conductance[which(!is.na(df_new$Specific_conductance))]
Su_vec_new<-df_new$Sulfate_total[which(!is.na(df_new$Sulfate_total))]

######################################################
four_test_p_val_generator<-function(sample_1,sample_2,alternative_val){
  p_val_vec<-rep(NA_real_,4)
  t_test_obj<-t.test(x = sample_1,y = sample_2,alternative = alternative_val)
  wilcox_test_obj<-wilcox.test(x = sample_1,y = sample_2,alternative = alternative_val)
  p_val_vec[1]<-t_test_obj$p.value
  p_val_vec[3]<-wilcox_test_obj$p.value
  # Getting the general and rank score vectors and group ID vector for permutation t and wilcoxon tests
  gen_score_vec<-c(sample_1,sample_2)
  rank_score_vec<-rank(c(sample_1,sample_2))
  group_ID_vec<-factor(rep(c("A", "B"), c(length(sample_1), length(sample_2))))
  ## Permutation tests
  p_val_vec[2]<-permTS(gen_score_vec~group_ID_vec,alternative=alternative_val, method="exact.mc",control=permControl(nmc=10^4-1))$p.value
  p_val_vec[4]<-permTS(rank_score_vec~group_ID_vec,alternative=alternative_val, method="exact.mc",control=permControl(nmc=10^4-1))$p.value
  return(p_val_vec)
}

######################################################
df_analytes_tests_two_sided<-data.frame("t"=numeric(),"tperm"=numeric(),"Wilcoxon"=numeric(),"Wilcoxonperm"=numeric())

df_analytes_tests_two_sided[1,]<-formatC(four_test_p_val_generator(sample_1 = Cl_vec_old,sample_2 = Cl_vec_new,alternative_val = "two.sided"), format = "e", digits = 2)
df_analytes_tests_two_sided[2,]<-formatC(four_test_p_val_generator(sample_1 = Mg_vec_old,sample_2 = Mg_vec_new,alternative_val = "two.sided"), format = "e", digits = 2)
df_analytes_tests_two_sided[3,]<-formatC(four_test_p_val_generator(sample_1 = pH_vec_old,sample_2 = pH_vec_new,alternative_val = "two.sided"), format = "e", digits = 2)
df_analytes_tests_two_sided[4,]<-formatC(four_test_p_val_generator(sample_1 = TDS_vec_old,sample_2 = TDS_vec_new,alternative_val = "two.sided"), format = "e", digits = 2)
df_analytes_tests_two_sided[5,]<-formatC(four_test_p_val_generator(sample_1 = SC_vec_old,sample_2 = SC_vec_new,alternative_val = "two.sided"), format = "e", digits = 2)
df_analytes_tests_two_sided[6,]<-formatC(four_test_p_val_generator(sample_1 = Su_vec_old,sample_2 = Su_vec_new,alternative_val = "two.sided"), format = "e", digits = 2)

write.table(df_analytes_tests_two_sided, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/df_analytes_tests_two_sided.txt", row.names=FALSE, sep = ",")

######################################################
df_analytes_tests_one_sided<-data.frame("t"=numeric(),"tperm"=numeric(),"Wilcoxon"=numeric(),"Wilcoxonperm"=numeric())

df_analytes_tests_one_sided[1,]<-formatC(four_test_p_val_generator(sample_1 = Cl_vec_old,sample_2 = Cl_vec_new,alternative_val = "less"), format = "e", digits = 2)
df_analytes_tests_one_sided[2,]<-formatC(four_test_p_val_generator(sample_1 = Mg_vec_old,sample_2 = Mg_vec_new,alternative_val = "less"), format = "e", digits = 2)
df_analytes_tests_one_sided[3,]<-formatC(four_test_p_val_generator(sample_1 = pH_vec_old,sample_2 = pH_vec_new,alternative_val = "less"), format = "e", digits = 2)
df_analytes_tests_one_sided[4,]<-formatC(four_test_p_val_generator(sample_1 = TDS_vec_old,sample_2 = TDS_vec_new,alternative_val = "less"), format = "e", digits = 2)
df_analytes_tests_one_sided[5,]<-formatC(four_test_p_val_generator(sample_1 = SC_vec_old,sample_2 = SC_vec_new,alternative_val = "less"), format = "e", digits = 2)
df_analytes_tests_one_sided[6,]<-formatC(four_test_p_val_generator(sample_1 = Su_vec_old,sample_2 = Su_vec_new,alternative_val = "less"), format = "e", digits = 2)

write.table(df_analytes_tests_one_sided, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/df_analytes_tests_one_sided.txt", row.names=FALSE, sep = ",")

#########################################################################################################
summary_generator<-function(sample_1, sample_2){
  return(c(round(mean(sample_1),2),round(mean(sample_2),2),round(median(sample_1),2),round(median(sample_2),2),round(sqrt(var(sample_1)),2),round(sqrt(var(sample_2)),2),length(sample_1), length(sample_2)))
}


df_analytes_summary<-data.frame("Mean Old"=numeric(),"Mean New"=numeric(),"Median Old"=numeric(),"Median New"=numeric(),"Std. Dev. Old"=numeric(),"Std. Dev. New"=numeric(), "No. of Obs. Old"=numeric(), "No. of Obs. New"=numeric())

df_analytes_summary[1,]<-summary_generator(sample_1 = Cl_vec_old,sample_2 = Cl_vec_new)
df_analytes_summary[2,]<-summary_generator(sample_1 = Mg_vec_old,sample_2 = Mg_vec_new)
df_analytes_summary[3,]<-summary_generator(sample_1 = pH_vec_old,sample_2 = pH_vec_new)
df_analytes_summary[4,]<-summary_generator(sample_1 = TDS_vec_old,sample_2 = TDS_vec_new)
df_analytes_summary[5,]<-summary_generator(sample_1 = SC_vec_old,sample_2 = SC_vec_new)
df_analytes_summary[6,]<-summary_generator(sample_1 = Su_vec_old,sample_2 = Su_vec_new)

write.table(df_analytes_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/df_analytes_summary.txt", row.names=FALSE, sep = ",")

#########################################################################################################
#########################################################################################################
## Tests for real data (Bradford)

## Reading and cleaning the old dataset
df_old<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/Bradford_Data/v2/Bradford_old_csv.csv",stringsAsFactors = F)
str(df_old)
list_analytes_old<-list()
k<-1
for(j in seq(3,50,by = 2)){
  non_censored_ids<-which((df_old[,j+1]=="nc")|(is.na(df_old[,j+1])))
  if(length(non_censored_ids)>0){
    non_censored_vec<-df_old[non_censored_ids,j]
    list_analytes_old[[k]]<-non_censored_vec[which(!is.na(non_censored_vec))]
    k<-k+1
  }
}
names(list_analytes_old)<-colnames(df_old[,seq(3,50,by = 2)])

## Reading and cleaning the new dataset
df_new<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/Bradford_Data/v1/Bradford_New_csv.csv",stringsAsFactors = F)
str(df_new)
list_analytes_new<-list()
k<-1
for(j in seq(6,41,by = 2)){
  non_censored_ids<-which((df_new[,j+1]=="nc")|(is.na(df_new[,j+1])))
  if(length(non_censored_ids)>0){
    non_censored_vec<-df_new[non_censored_ids,j]
    list_analytes_new[[k]]<-non_censored_vec[which(!is.na(non_censored_vec))]
    k<-k+1
  }
}
names(list_analytes_new)<-colnames(df_new[,seq(6,41,by = 2)])

names(list_analytes_old)
## Reassigning the names
names(list_analytes_old)<-c("Specific_conductance","Hardness","pH_acidity","Calcium","Magnesium","Sodium","Potassium","Alkalinity","Sulfate","Chloride","Fluoride","TDS","Nitrate","Aluminium","Arsenic","Barium","Cadmium","Cromium","Iron","Lead","Manganese","Nickel", "Strontium","Zinc")

names(list_analytes_new)<-c("Alkalinity","Arsenic","Barium","Calcium","Chloride","Hardness","Iron","Lead","Magnesium","Manganese","Methane","pH_acidity","Potassium","Sodium","TDS","Specific_conductance","Sulfate","Turbidity")

######################################################
## Function to get p-values from seven tests 
test_p_val_generator<-function(sample_1,sample_2,alternative_val,boot_size){
  p_val_vec<-rep(NA_real_,11)
  names(p_val_vec)<-c("t_pool","t_Welch","t_perm","t_boot_unpooled","t_boot_pooled","W","W_boot_unpooled","W_boot_pooled","BM","BM_boot_unpooled","BM_boot_pooled")
  p_val_vec["t_pool"]<-t.test(x = sample_1,y = sample_2,alternative = alternative_val,var.equal=TRUE)$p.value
  t_test_obj<-t.test(x = sample_1,y = sample_2,alternative = alternative_val)
  wilcox_test_obj<-wilcox.test(x = sample_1,y = sample_2,alternative = alternative_val)
  BM_test_obj<-brunner.munzel.test(x = sample_1, y = sample_2, alternative = alternative_val)
  p_val_vec["t_Welch"]<-t_test_obj$p.value
  p_val_vec["W"]<-wilcox_test_obj$p.value
  p_val_vec["BM"]<-BM_test_obj$p.value
  # Getting the general and rank score vectors and group ID vector for permutation t and wilcoxon tests
  gen_score_vec<-c(sample_1,sample_2)
  group_ID_vec<-factor(rep(c("A", "B"), c(length(sample_1), length(sample_2))))
  ## Permutation test
  p_val_vec["t_perm"]<-permTS(gen_score_vec~group_ID_vec,alternative=alternative_val, method="exact.mc",control=permControl(nmc=10^4-1))$p.value
  ## Bootstrap tests
  boot_output<-boot_test(sample_1 = sample_1,sample_2 = sample_2,test_stat_vec = c(t_test_obj$statistic,wilcox_test_obj$statistic,BM_test_obj$statistic),alternative_val = alternative_val,boot_size = boot_size)
  p_val_vec["t_boot_unpooled"]<-boot_output[1]
  p_val_vec["t_boot_pooled"]<-boot_output[2]
  p_val_vec["W_boot_unpooled"]<-boot_output[3]
  p_val_vec["W_boot_pooled"]<-boot_output[4]
  p_val_vec["BM_boot_unpooled"]<-boot_output[5]
  p_val_vec["BM_boot_pooled"]<-boot_output[6]
  p_val_vec<-formatC(p_val_vec, format = "e", digits = 2)
  return(p_val_vec)
}

######################################################
## Conducting two sided tests and one_sided tests
df_analytes_tests_two_sided<-data.frame("Analyte"=character(),"t_Pool"=numeric(),"t"=numeric(),"tperm"=numeric(),"tboot_UnPooled"=numeric(),"tboot_Pooled"=numeric(),"Wilcoxon"=numeric(),"Wilcoxonboot_Unpooled"=numeric(),"Wilcoxonboot_Pooled"=numeric(),"BM"=numeric(),"BMboot_UnPooled"=numeric(),"BMboot_Pooled"=numeric(),stringsAsFactors = F)
df_analytes_tests_less<-data.frame("Analyte"=character(),"t_Pool"=numeric(),"t"=numeric(),"tperm"=numeric(),"tboot_UnPooled"=numeric(),"tboot_Pooled"=numeric(),"Wilcoxon"=numeric(),"Wilcoxonboot_Unpooled"=numeric(),"Wilcoxonboot_Pooled"=numeric(),"BM"=numeric(),"BMboot_UnPooled"=numeric(),"BMboot_Pooled"=numeric(),stringsAsFactors = F)
df_analytes_tests_greater<-data.frame("Analyte"=character(),"t_Pool"=numeric(),"t"=numeric(),"tperm"=numeric(),"tboot_UnPooled"=numeric(),"tboot_Pooled"=numeric(),"Wilcoxon"=numeric(),"Wilcoxonboot_Unpooled"=numeric(),"Wilcoxonboot_Pooled"=numeric(),"BM"=numeric(),"BMboot_UnPooled"=numeric(),"BMboot_Pooled"=numeric(),stringsAsFactors = F)
# for(i in 1:length(list_analytes_old)){
#   print(c(names(list_analytes_old)[i],agrep(names(list_analytes_old)[i],names(list_analytes_new),value=T)))
# }
name_IDs_new<-rep(NA_integer_,length(list_analytes_old))
k<-1
for(i in 1:length(list_analytes_old)){
  sample_1<-list_analytes_old[[i]]
  if(length(which(names(list_analytes_new)==names(list_analytes_old)[i]))>0){
    name_IDs_new[i]<-which(names(list_analytes_new)==names(list_analytes_old)[i])
    sample_2<-list_analytes_new[[name_IDs_new[i]]]
    ## Filling out the analyte names
    df_analytes_tests_two_sided[k,1]<-names(list_analytes_old)[i]
    df_analytes_tests_less[k,1]<-names(list_analytes_old)[i]
    df_analytes_tests_greater[k,1]<-names(list_analytes_old)[i]
    ## Filling out the p values
    df_analytes_tests_two_sided[k,-1]<-test_p_val_generator(sample_1 = sample_1,sample_2 = sample_2,alternative_val = "two.sided",boot_size = 1000)
    df_analytes_tests_less[k,-1]<-test_p_val_generator(sample_1 = sample_1,sample_2 = sample_2,alternative_val = "less",boot_size = 1000)
    df_analytes_tests_greater[k,-1]<-test_p_val_generator(sample_1 = sample_1,sample_2 = sample_2,alternative_val = "greater",boot_size = 1000)
    k<-k+1
  }
  print(i)
}

write.table(df_analytes_tests_two_sided, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v3/df_Bradford_two_sided.txt", row.names=FALSE, sep = ",")

write.table(df_analytes_tests_less, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v3/df_Bradford_less.txt", row.names=FALSE, sep = ",")

write.table(df_analytes_tests_greater, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v3/df_Bradford_greater.txt", row.names=FALSE, sep = ",")

#########################################################################################################
## Writing a function to extract the summary statustics
summary_generator<-function(sample_1, sample_2){
  return(c(round(mean(sample_1),2),round(mean(sample_2),2),round(median(sample_1),2),round(median(sample_2),2),round(sqrt(var(sample_1)),2),round(sqrt(var(sample_2)),2),length(sample_1), length(sample_2)))
}

df_analytes_summary<-data.frame("Analyte"=character(),"Mean Old"=numeric(),"Mean New"=numeric(),"Median Old"=numeric(),"Median New"=numeric(),"Std. Dev. Old"=numeric(),"Std. Dev. New"=numeric(), "No. of Obs. Old"=numeric(), "No. of Obs. New"=numeric(),stringsAsFactors = F)

name_IDs_new<-rep(NA_integer_,length(list_analytes_old))
k<-1
for(i in 1:length(list_analytes_old)){
  sample_1<-list_analytes_old[[i]]
  if(length(which(names(list_analytes_new)==names(list_analytes_old)[i]))>0){
    name_IDs_new[i]<-which(names(list_analytes_new)==names(list_analytes_old)[i])
    sample_2<-list_analytes_new[[name_IDs_new[i]]]
    df_analytes_summary[k,1]<-names(list_analytes_old)[i]
    df_analytes_summary[k,-1]<-summary_generator(sample_1 = sample_1,sample_2 = sample_2)
    k<-k+1
  }
}

write.table(df_analytes_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v3/df_Bradford_summary.txt", row.names=FALSE, sep = ",")

#########################################################################################################
#########################################################################################################
#########################################################################################################
## Tests for real data (NW PA revised)
## Reading and cleaning the old dataset
df_old<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/NW_PA_Data/v1/Old_NW_PA_csv.csv",stringsAsFactors = F)
str(df_old)
list_analytes_old<-list()
k<-1
for(j in seq(6,35,by = 2)){
  non_censored_ids<-which((df_old[,j+1]=="nc")|(is.na(df_old[,j+1])))
  if(length(non_censored_ids)>0){
    non_censored_vec<-df_old[non_censored_ids,j]
    list_analytes_old[[k]]<-non_censored_vec[which(!is.na(non_censored_vec))]
    k<-k+1
  }
}
names(list_analytes_old)<-colnames(df_old[,seq(6,35,by = 2)])

## Reading and cleaning the new dataset
df_new<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/NW_PA_Data/v3/New NW PA_csv.csv",stringsAsFactors = F)
str(df_new)
list_analytes_new<-list()
k<-1
for(j in seq(5,34,by = 2)){
  non_censored_ids<-which((df_new[,j+1]=="nc")|(is.na(df_new[,j+1])))
  if(length(non_censored_ids)>0){
    non_censored_vec<-df_new[non_censored_ids,j]
    list_analytes_new[[k]]<-non_censored_vec[which(!is.na(non_censored_vec))]
    k<-k+1
  }
}
names(list_analytes_new)<-colnames(df_new[,seq(5,34,by = 2)])

names(list_analytes_old)
## Reassigning the names
names(list_analytes_old)<-c("pH_acidity","Hardness","TDS","Alkalinity","Specific_conductance","Potassium","Magnesium","Calcium","Chloride","Sodium","Sulfate","Methane","Turbidity","Iron","Manganese")

names(list_analytes_new)
names(list_analytes_new)<-c("Alkalinity","Calcium","Chloride","Hardness","Iron","Manganese","Methane","pH_acidity","Potassium","Sodium","Specific_conductance","Sulfate","TDS","Turbidity","Magnesium")

######################################################
## Function to get p-values from seven tests 
test_p_val_generator<-function(sample_1,sample_2,alternative_val,boot_size){
  p_val_vec<-rep(NA_real_,11)
  names(p_val_vec)<-c("t_pool","t_Welch","t_perm","t_boot_unpooled","t_boot_pooled","W","W_boot_unpooled","W_boot_pooled","BM","BM_boot_unpooled","BM_boot_pooled")
  p_val_vec["t_pool"]<-t.test(x = sample_1,y = sample_2,alternative = alternative_val,var.equal=TRUE)$p.value
  t_test_obj<-t.test(x = sample_1,y = sample_2,alternative = alternative_val)
  wilcox_test_obj<-wilcox.test(x = sample_1,y = sample_2,alternative = alternative_val)
  BM_test_obj<-brunner.munzel.test(x = sample_1, y = sample_2, alternative = alternative_val)
  p_val_vec["t_Welch"]<-t_test_obj$p.value
  p_val_vec["W"]<-wilcox_test_obj$p.value
  p_val_vec["BM"]<-BM_test_obj$p.value
  # Getting the general and rank score vectors and group ID vector for permutation t and wilcoxon tests
  gen_score_vec<-c(sample_1,sample_2)
  group_ID_vec<-factor(rep(c("A", "B"), c(length(sample_1), length(sample_2))))
  ## Permutation test
  p_val_vec["t_perm"]<-permTS(gen_score_vec~group_ID_vec,alternative=alternative_val, method="exact.mc",control=permControl(nmc=10^4-1))$p.value
  ## Bootstrap tests
  boot_output<-boot_test(sample_1 = sample_1,sample_2 = sample_2,test_stat_vec = c(t_test_obj$statistic,wilcox_test_obj$statistic,BM_test_obj$statistic),alternative_val = alternative_val,boot_size = boot_size)
  p_val_vec["t_boot_unpooled"]<-boot_output[1]
  p_val_vec["t_boot_pooled"]<-boot_output[2]
  p_val_vec["W_boot_unpooled"]<-boot_output[3]
  p_val_vec["W_boot_pooled"]<-boot_output[4]
  p_val_vec["BM_boot_unpooled"]<-boot_output[5]
  p_val_vec["BM_boot_pooled"]<-boot_output[6]
  p_val_vec<-formatC(p_val_vec, format = "e", digits = 2)
  return(p_val_vec)
}

######################################################
## Conducting two sided tests and one_sided tests
df_analytes_tests_two_sided<-data.frame("Analyte"=character(),"t_Pool"=numeric(),"t"=numeric(),"tperm"=numeric(),"tboot_UnPooled"=numeric(),"tboot_Pooled"=numeric(),"Wilcoxon"=numeric(),"Wilcoxonboot_Unpooled"=numeric(),"Wilcoxonboot_Pooled"=numeric(),"BM"=numeric(),"BMboot_UnPooled"=numeric(),"BMboot_Pooled"=numeric(),stringsAsFactors = F)
df_analytes_tests_less<-data.frame("Analyte"=character(),"t_Pool"=numeric(),"t"=numeric(),"tperm"=numeric(),"tboot_UnPooled"=numeric(),"tboot_Pooled"=numeric(),"Wilcoxon"=numeric(),"Wilcoxonboot_Unpooled"=numeric(),"Wilcoxonboot_Pooled"=numeric(),"BM"=numeric(),"BMboot_UnPooled"=numeric(),"BMboot_Pooled"=numeric(),stringsAsFactors = F)
df_analytes_tests_greater<-data.frame("Analyte"=character(),"t_Pool"=numeric(),"t"=numeric(),"tperm"=numeric(),"tboot_UnPooled"=numeric(),"tboot_Pooled"=numeric(),"Wilcoxon"=numeric(),"Wilcoxonboot_Unpooled"=numeric(),"Wilcoxonboot_Pooled"=numeric(),"BM"=numeric(),"BMboot_UnPooled"=numeric(),"BMboot_Pooled"=numeric(),stringsAsFactors = F)
# for(i in 1:length(list_analytes_old)){
#   print(c(names(list_analytes_old)[i],agrep(names(list_analytes_old)[i],names(list_analytes_new),value=T)))
# }
name_IDs_new<-rep(NA_integer_,length(list_analytes_old))
k<-1
for(i in 1:length(list_analytes_old)){
  sample_1<-list_analytes_old[[i]]
  if(length(which(names(list_analytes_new)==names(list_analytes_old)[i]))>0){
    name_IDs_new[i]<-which(names(list_analytes_new)==names(list_analytes_old)[i])
    sample_2<-list_analytes_new[[name_IDs_new[i]]]
    ## Filling out the analyte names
    df_analytes_tests_two_sided[k,1]<-names(list_analytes_old)[i]
    df_analytes_tests_less[k,1]<-names(list_analytes_old)[i]
    df_analytes_tests_greater[k,1]<-names(list_analytes_old)[i]
    ## Filling out the p values
    df_analytes_tests_two_sided[k,-1]<-test_p_val_generator(sample_1 = sample_1,sample_2 = sample_2,alternative_val = "two.sided",boot_size = 1000)
    df_analytes_tests_less[k,-1]<-test_p_val_generator(sample_1 = sample_1,sample_2 = sample_2,alternative_val = "less",boot_size = 1000)
    df_analytes_tests_greater[k,-1]<-test_p_val_generator(sample_1 = sample_1,sample_2 = sample_2,alternative_val = "greater",boot_size = 1000)
    k<-k+1
  }
  print(i)
}

write.table(df_analytes_tests_two_sided, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v3/df_NW_PA_two_sided.txt", row.names=FALSE, sep = ",")

write.table(df_analytes_tests_less, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v3/df_NW_PA_less.txt", row.names=FALSE, sep = ",")

write.table(df_analytes_tests_greater, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v3/df_NW_PA_greater.txt", row.names=FALSE, sep = ",")

#########################################################################################################
## Writing a function to extract the summary statustics
## Writing a function to extract the summary statustics
summary_generator<-function(sample_1, sample_2){
  return(c(round(mean(sample_1),2),round(mean(sample_2),2),round(median(sample_1),2),round(median(sample_2),2),round(sqrt(var(sample_1)),2),round(sqrt(var(sample_2)),2),length(sample_1), length(sample_2)))
}

df_analytes_summary<-data.frame("Analyte"=character(),"Mean Old"=numeric(),"Mean New"=numeric(),"Median Old"=numeric(),"Median New"=numeric(),"Std. Dev. Old"=numeric(),"Std. Dev. New"=numeric(), "No. of Obs. Old"=numeric(), "No. of Obs. New"=numeric(),stringsAsFactors = F)

name_IDs_new<-rep(NA_integer_,length(list_analytes_old))
k<-1
for(i in 1:length(list_analytes_old)){
  sample_1<-list_analytes_old[[i]]
  if(length(which(names(list_analytes_new)==names(list_analytes_old)[i]))>0){
    name_IDs_new[i]<-which(names(list_analytes_new)==names(list_analytes_old)[i])
    sample_2<-list_analytes_new[[name_IDs_new[i]]]
    df_analytes_summary[k,1]<-names(list_analytes_old)[i]
    df_analytes_summary[k,-1]<-summary_generator(sample_1 = sample_1,sample_2 = sample_2)
    k<-k+1
  }
}

write.table(df_analytes_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v3/df_NW_PA_summary.txt", row.names=FALSE, sep = ",")

#########################################################################################################


