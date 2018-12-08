
## Defining the sample size vectors for the two samples
sample_1_size_vec<-c(20,50,100,500,1000,10000)
sample_2_size_vec<-c(20,50,100,500,1000,10000)

## Creating all combinations of sample sizes
df_sample_size<-expand.grid(x= sample_1_size_vec,y=sample_2_size_vec)

## Initializing the power
df_sample_size$power<-NA_real_

ditribution_vec<-c("normal","t")

## test type could be simple, permutation and bootstrap

## test name could be t, wilcoxon 

test_results_generator<-function(sample_1_size,sample_2_size,dist_type,mean_sample_1,mean_sample_2,n_simulations){
  ## Initializing the p-value vector
  p_val_t_vec<-rep(NA_real_,n_simulations)
  p_val_wilcox_vec<-rep(NA_real_,n_simulations)
  ## Loop to conduct simulations
  for (j in 1:n_simulations){
    sample_1<-rnorm(n = sample_1_size,mean = mean_sample_1,sd = 1)
    sample_2<-rnorm(n = sample_1_size,mean = mean_sample_2,sd = 1)
    p_val_t_vec[j]<-t.test(x = sample_1,y = sample_2)$p.value
    p_val_wilcox_vec[j]<-wilcox.test(x = sample_1,y = sample_2)$p.value
  }
  power_size_t_val<-sum(p_val_t_vec<0.05)/n_simulations
  power_size_wilcox_val<-sum(p_val_wilcox_vec<0.05)/n_simulations
  return(c(power_size_t_val,power_size_wilcox_val))
}



## Loop over different sample size combinations
for (i in 1:nrow(df_sample_size)){
  ## Defining the number of simulations
  n_simulations<-1000
  ## Initializing the p-value vector
  p_val_vec<-rep(NA_real_,n_simulations)
  ## Loop to conduct simulations
  for (j in 1:n_simulations){
    sample_1<-rnorm(n = df_sample_size[i,1],mean = 0,sd = 1)
    sample_2<-rnorm(n = df_sample_size[i,2],mean = 0.5,sd = 1)
    p_val_vec[j]<-wilcox.test(x = sample_1,y = sample_2)$p.value
  }
  df_sample_size$power[i]<-sum(p_val_vec<0.05)/n_simulations
  print(i)
}

names(df_sample_size)[c(1,2)]<-c("Sample_size_1","Sample_size_2")
df_sample_size





