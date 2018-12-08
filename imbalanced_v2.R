#########################################################################################################
## Loading the required packages
library(perm)
library(foreach)
library(doParallel)

#########################################################################################################
## Writing a function to conduct Bootstrap t test and Wilcoxon test
boot_test<-function(sample_1,sample_2,test_stat_vec,boot_size){
  sample_combined<-c(sample_1,sample_2)
  sample_1_normalized<-sample_1-mean(sample_1)+mean(sample_combined)
  sample_2_normalized<-sample_2-mean(sample_2)+mean(sample_combined)
  t_val_vec<-rep(NA_real_,boot_size)
  wilcox_val_vec<-rep(NA_real_,boot_size)
  #ptm<-proc.time()
  statistic_mat<-cbind(t_val_vec,wilcox_val_vec)
  statistic_mat_output<-t(apply(X = statistic_mat,MARGIN = 1,FUN = function(x){
    sample_boot_1_t<-sample(sample_1_normalized,replace = T)
    sample_boot_2_t<-sample(sample_2_normalized,replace = T)
    sample_boot_combined<-sample(x = sample_combined,replace = T)
    sample_boot_1_wilcox<-sample_boot_combined[1:length(sample_1)]
    sample_boot_2_wilcox<-sample_boot_combined[(length(sample_1)+1):(length(sample_1)+length(sample_2))]
    t_val<-t.test(x = sample_boot_1_t,y = sample_boot_2_t,alternative = "two.sided")$statistic
    wilcox_val<-wilcox.test(x = sample_boot_1_wilcox,y = sample_boot_2_wilcox,alternative = "two.sided")$statistic
    return(c(t_val,wilcox_val))
  }))
  #print(proc.time()-ptm)
  boot_p_val_vec<-c((sum(abs(statistic_mat_output[,1])>=abs(test_stat_vec[1]))/boot_size),(sum((statistic_mat_output[,2])>=(test_stat_vec[2]))/boot_size))
  return(boot_p_val_vec)
}

#########################################################################################################
## Writing the function to generate table for test results
test_results_generator<-function(n_simulations,sample_1_size,sample_2_size,dist_type,sample_1_mean=NA,sample_2_mean=NA,sample_1_var=NA,sample_2_var=NA,sample_1_df=NA,sample_2_df=NA,sample_1_ncp=NA,sample_2_ncp=NA,sample_1_shape=NA, sample_2_shape=NA, sample_1_rate=NA, sample_2_rate=NA, sample_1_mix_p = NA, sample_2_mix_p = NA, sample_1_mixComp1_rate = NA, sample_1_mixComp2_rate = NA, sample_2_mixComp1_rate = NA, sample_2_mixComp2_rate = NA, sample_1_mixComp1_shape = NA, sample_1_mixComp2_shape = NA, sample_2_mixComp1_shape = NA, sample_2_mixComp2_shape = NA, boot_size){
  ## Initializing the p-value dataframes for t and wilcoxon test for simple, permutation and bootstrap cases
  df_t_p_val<-data.frame(matrix(NA_real_,n_simulations,3))
  names(df_t_p_val)<-c("Simple","Permutation","Bootstrap")
  df_wilcox_p_val<-data.frame(matrix(NA_real_,n_simulations,3))
  names(df_wilcox_p_val)<-c("Simple","Permutation","Bootstrap")
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
      sample_1<-rep(NA_real_,sample_1_size)
      sample_2<-rep(NA_real_,sample_2_size)
      N <- 100000
      
      components <- sample(1:2,prob=c(0.3,0.5,0.2),size=N,replace=TRUE)
      mus <- c(0,10,3)
      sds <- sqrt(c(1,1,0.1))
      
      samples <- rnorm(n=N,mean=mus[components],sd=sds[components])
      
      for(k in 1:sample_1_size){
        p1<-rbinom(n = 1,size = 1,prob = sample_1_mix_p)
        if(p1==1){
          sample_1[k]<-rgamma(n = 1,rate = sample_1_mixComp1_rate,shape=sample_1_mixComp1_shape)
        }else{
          sample_1[k]<-rgamma(n = 1,rate = sample_1_mixComp2_rate,shape=sample_1_mixComp2_shape)
        }
      }
      for(k in 1:sample_2_size){
        p1<-rbinom(n = 1,size = 1,prob = sample_2_mix_p)
        if(p1==1){
          sample_2[k]<-rgamma(n = 1,rate = sample_2_mixComp1_rate,shape=sample_2_mixComp1_shape)
        }else{
          sample_2[k]<-rgamma(n = 1,rate = sample_2_mixComp2_rate,shape=sample_2_mixComp2_shape)
        }
      }
    }
    ## Simple tests
    t_test_obj<-t.test(x = sample_1,y = sample_2)
    wilcox_test_obj<-wilcox.test(x = sample_1,y = sample_2)
    df_t_p_val[j,1]<-t_test_obj$p.value
    df_wilcox_p_val[j,1]<-wilcox_test_obj$p.value
    # Getting the general and rank score vectors and group ID vector for permutation t and wilcoxon tests
    gen_score_vec<-c(sample_1,sample_2)
    rank_score_vec<-rank(c(sample_1,sample_2))
    group_ID_vec<-factor(rep(c("A", "B"), c(length(sample_1), length(sample_2))))
    ## Permutation tests
    df_t_p_val[j,2]<-permTS(gen_score_vec~group_ID_vec,alternative="two.sided", method="exact.mc",control=permControl(nmc=10^4-1))$p.value
    df_wilcox_p_val[j,2]<-permTS(rank_score_vec~group_ID_vec,alternative="two.sided", method="exact.mc",control=permControl(nmc=10^4-1))$p.value
    ## Bootstrap tests
    boot_test_output<-boot_test(sample_1 = sample_1,sample_2 = sample_2,test_stat = c(t_test_obj$statistic,wilcox_test_obj$statistic),boot_size = boot_size)
    df_t_p_val[j,3]<-boot_test_output[1]
    df_wilcox_p_val[j,3]<-boot_test_output[2]
    print(j)
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
## Tests for real data
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