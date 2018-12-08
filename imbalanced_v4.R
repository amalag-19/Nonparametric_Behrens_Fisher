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
#########################################################################################################
## Tests for real data (Bradford)

## Reading and cleaning the old dataset
df_old<-read.csv(file = "/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/Data/Bradford_Data/v2/Bradford_old_csv.csv",stringsAsFactors = F)
str(df_old)
list_analytes_old<-list()
k<-1
for(j in seq(3,50,by = 2)){
  obs_vec<-df_old[which(!is.na(df_old[,j])),j]
  list_analytes_old[[k]]<-obs_vec
  k<-k+1
}
names(list_analytes_old)<-colnames(df_old[,seq(3,50,by = 2)])

## Reading and cleaning the new dataset
df_new<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/Bradford_Data/v1/Bradford_New_csv.csv",stringsAsFactors = F)
str(df_new)
list_analytes_new<-list()
k<-1
for(j in seq(6,41,by = 2)){
  obs_vec<-df_new[which(!is.na(df_new[,j])),j]
  list_analytes_new[[k]]<-obs_vec
  k<-k+1
}
names(list_analytes_new)<-colnames(df_new[,seq(6,41,by = 2)])

names(list_analytes_old)
## Reassigning the names
names(list_analytes_old)<-c("Specific_conductance","Hardness","pH_acidity","Calcium","Magnesium","Sodium","Potassium","Alkalinity","Sulfate","Chloride","Fluoride","TDS","Nitrate","Aluminium","Arsenic","Barium","Cadmium","Cromium","Iron","Lead","Manganese","Nickel", "Strontium","Zinc")

names(list_analytes_new)
names(list_analytes_new)<-c("Alkalinity","Arsenic","Barium","Calcium","Chloride","Hardness","Iron","Lead","Magnesium","Manganese","Methane","pH_acidity","Potassium","Sodium","TDS","Specific_conductance","Sulfate","Turbidity")

#########################################################################################################
## Writing a function to get cdf plots

library(ggplot2)

plot_generator<-function(analyte_name,xmin=-1,xmax){
  sample_1<-list_analytes_old[[which(names(list_analytes_old)==analyte_name)]]
  sample_2<-list_analytes_new[[which(names(list_analytes_new)==analyte_name)]]
  df_sample<-data.frame("ID"=c(rep("Before 2010",length(sample_1)),rep("After 2010",length(sample_2))),"Value"=c(sample_1,sample_2),stringsAsFactors = F)
  df_sample$ID<-factor(df_sample$ID,levels=c("Before 2010","After 2010"))
  cdf_plot_out<-ggplot(df_sample, aes(x = Value, colour = factor(ID))) + stat_ecdf() +coord_cartesian(xlim = c(xmin,xmax))+theme_bw()+labs(x="Observed Value", y="Fraction of Data")+theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25),legend.text=element_text(size=20),strip.text.x = element_text(size = 30),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),legend.title = element_text(hjust = 0.5,size = 20))+scale_color_discrete(guide = guide_legend(title = "Temporal Range",labels=c("Before 2010","After 2010"),keywidth=0.2,keyheight=0.4,default.unit = "inch"))
  
  df_summary<-rbind(round(summary(sample_1),2),round(summary(sample_2),2))
  
  sample_1_modified<-sample_1[which(sample_1<1000)]
  df_sample<-data.frame("ID"=c(rep("Before 2010",length(sample_1_modified)),rep("After 2010",length(sample_2))),"Value"=c(sample_1_modified,sample_2),stringsAsFactors = F)
  df_sample$ID<-factor(df_sample$ID,levels=c("Before 2010","After 2010"))
  
  violin_plot_out<-ggplot(df_sample, aes(x = factor(ID), y = Value)) + geom_violin()+coord_cartesian(ylim = c(xmin,xmax))+theme_bw()+labs(x="Temporal Range", y="Observed Values")+theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25),legend.text=element_text(size=20),strip.text.x = element_text(size = 30),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),legend.title = element_text(hjust = 0.5,size = 20))# +coord_cartesian(xlim = c(-1,xmax))+theme_bw()+labs(x="Observed Value", y="Fraction of Data")+theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25),legend.text=element_text(size=20),strip.text.x = element_text(size = 30),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),legend.title = element_text(hjust = 0.5,size = 20))+scale_color_discrete(guide = guide_legend(title = "Temporal Range",labels=c("Before 2010","After 2010"),keywidth=0.2,keyheight=0.4,default.unit = "inch"))
  return(list(cdf_plot_out,violin_plot_out,sample_1,sample_2,df_summary))
}

#undebug(plot_generator)
out<-plot_generator(analyte_name="Potassium",xmax = 10)
out[[2]]
df_summary<-out[[5]]
write.table(df_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/summary_stats_Bradford/K_summary.txt", row.names=FALSE, sep = ",")

out<-plot_generator(analyte_name="Sulfate",xmax = 200)
out[[2]]
df_summary<-out[[5]]
write.table(df_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/summary_stats_Bradford/SO4_summary.txt", row.names=FALSE, sep = ",")


out<-plot_generator(analyte_name="Arsenic",xmin=0,xmax = 0.1)
out[[2]]
df_summary<-out[[5]]
write.table(df_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/summary_stats_Bradford/As_summary.txt", row.names=FALSE, sep = ",")


out<-plot_generator(analyte_name="Barium",xmin=-0.25,xmax = 4)
out[[2]]
df_summary<-out[[5]]
write.table(df_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/summary_stats_Bradford/Ba_summary.txt", row.names=FALSE, sep = ",")

out<-plot_generator(analyte_name="Iron",xmin=-0.25,xmax = 5)
out[[2]]
df_summary<-out[[5]]
write.table(df_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/summary_stats_Bradford/Fe_summary.txt", row.names=FALSE, sep = ",")

out<-plot_generator(analyte_name="Lead",xmin=0,xmax = 0.05)
out[[2]]
df_summary<-out[[5]]
write.table(df_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/summary_stats_Bradford/Pb_summary.txt", row.names=FALSE, sep = ",")


out<-plot_generator(analyte_name="Manganese",xmin=-0.25,xmax = 2)
out[[2]]
df_summary<-out[[5]]
write.table(df_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/summary_stats_Bradford/Mn_summary.txt", row.names=FALSE, sep = ",")


#out<-cdf_plot_generator(analyte_name="pH_acidity",xmin=5,xmax = 10)

#out[[1]]


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

write.table(df_analytes_tests_two_sided, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/df_Bradford_two_sided.txt", row.names=FALSE, sep = ",")

write.table(df_analytes_tests_less, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/df_Bradford_less.txt", row.names=FALSE, sep = ",")

write.table(df_analytes_tests_greater, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/df_Bradford_greater.txt", row.names=FALSE, sep = ",")

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

write.table(df_analytes_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/df_Bradford_summary.txt", row.names=FALSE, sep = ",")

#########################################################################################################
#########################################################################################################
#########################################################################################################
## Tests for real data (NW PA revised)
## Reading and cleaning the old dataset
df_old<-read.csv(file = "/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/Data/NW_PA_Data/v1/Old_NW_PA_csv.csv",stringsAsFactors = F)
str(df_old)
list_analytes_old<-list()
k<-1
for(j in seq(6,35,by = 2)){
  obs_vec<-df_old[which(!is.na(df_old[,j])),j]
  list_analytes_old[[k]]<-obs_vec
  k<-k+1
}
names(list_analytes_old)<-colnames(df_old[,seq(6,35,by = 2)])

## Reading and cleaning the new dataset
df_new<-read.csv(file = "/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/Data/NW_PA_Data/v3/New NW PA_csv.csv",stringsAsFactors = F)
str(df_new)
list_analytes_new<-list()
k<-1
for(j in seq(5,34,by = 2)){
  obs_vec<-df_new[which(!is.na(df_new[,j])),j]
  list_analytes_new[[k]]<-obs_vec
  k<-k+1
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

write.table(df_analytes_tests_two_sided, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/df_NW_PA_two_sided.txt", row.names=FALSE, sep = ",")

write.table(df_analytes_tests_less, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/df_NW_PA_less.txt", row.names=FALSE, sep = ",")

write.table(df_analytes_tests_greater, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/df_NW_PA_greater.txt", row.names=FALSE, sep = ",")

#########################################################################################################
## Writing a function to extract the summary statistics
## Writing a function to extract the summary statistics
summary_generator<-function(sample_1, sample_2){
  summary_old<-c(round(min(sample_1),2),round(mean(sample_1),2),round(median(sample_1),2),round(quantile(x = sample_1, probs = 0.95),2),round(max(sample_1),2),round(sqrt(var(sample_1)),2),length(sample_1))
  summary_new<-c(round(min(sample_2),2),round(mean(sample_2),2),round(median(sample_2),2),round(quantile(x = sample_2, probs = 0.95),2),round(max(sample_2),2),round(sqrt(var(sample_2)),2),length(sample_2))
  return(list(summary_old,summary_new))
}

df_analytes_summary_old<-data.frame("Analyte"=character(),"Minimum"=numeric(),"Mean"=numeric(),"Median"=numeric(),"95th Percentile"=numeric(),"Maximum"=numeric(),"Std. Dev."=numeric(), "No. of Obs."=numeric(),stringsAsFactors = F)

df_analytes_summary_new<-data.frame("Analyte"=character(),"Minimum"=numeric(),"Mean"=numeric(),"Median"=numeric(),"95th Percentile"=numeric(),"Maximum"=numeric(),"Std. Dev."=numeric(), "No. of Obs."=numeric(),stringsAsFactors = F)

name_IDs_new<-rep(NA_integer_,length(list_analytes_old))
k<-1
for(i in 1:length(list_analytes_old)){
  sample_1<-list_analytes_old[[i]]
  if(length(which(names(list_analytes_new)==names(list_analytes_old)[i]))>0){
    name_IDs_new[i]<-which(names(list_analytes_new)==names(list_analytes_old)[i])
    sample_2<-list_analytes_new[[name_IDs_new[i]]]
    df_analytes_summary_old[k,1]<-names(list_analytes_old)[i]
    df_analytes_summary_new[k,1]<-names(list_analytes_old)[i]
    df_analytes_summary_old[k,-1]<-summary_generator(sample_1 = sample_1,sample_2 = sample_2)[[1]]
    df_analytes_summary_new[k,-1]<-summary_generator(sample_1 = sample_1,sample_2 = sample_2)[[2]]
    k<-k+1
  }
}

write.table(df_analytes_summary_old, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/df_NW_PA_summary_old.txt", row.names=FALSE, sep = ",")

write.table(df_analytes_summary_new, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/df_NW_PA_summary_new.txt", row.names=FALSE, sep = ",")

#########################################################################################################
## Writing a function to get cdf plots

library(ggplot2)

plot_generator<-function(analyte_name,xmin=-1,xmax){
  sample_1<-list_analytes_old[[which(names(list_analytes_old)==analyte_name)]]
  sample_2<-list_analytes_new[[which(names(list_analytes_new)==analyte_name)]]
  df_sample<-data.frame("ID"=c(rep("Before 2010",length(sample_1)),rep("After 2010",length(sample_2))),"Value"=c(sample_1,sample_2),stringsAsFactors = F)
  df_sample$ID<-factor(df_sample$ID,levels=c("Before 2010","After 2010"))
  cdf_plot_out<-ggplot(df_sample, aes(x = Value, colour = factor(ID))) + stat_ecdf() +coord_cartesian(xlim = c(xmin,xmax))+theme_bw()+labs(x="Observed Value", y="Fraction of Data")+theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25),legend.text=element_text(size=20),strip.text.x = element_text(size = 30),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),legend.title = element_text(hjust = 0.5,size = 20))+scale_color_discrete(guide = guide_legend(title = "Temporal Range",labels=c("Before 2010","After 2010"),keywidth=0.2,keyheight=0.4,default.unit = "inch"))
  
  df_summary<-rbind(round(summary(sample_1),2),round(summary(sample_2),2))
  
  sample_1_modified<-sample_1[which(sample_1<1000)]
  df_sample<-data.frame("ID"=c(rep("Before 2010",length(sample_1_modified)),rep("After 2010",length(sample_2))),"Value"=c(sample_1_modified,sample_2),stringsAsFactors = F)
  df_sample$ID<-factor(df_sample$ID,levels=c("Before 2010","After 2010"))
  
  violin_plot_out<-ggplot(df_sample, aes(x = factor(ID), y = Value)) + geom_violin()+coord_cartesian(ylim = c(xmin,xmax))+theme_bw()+labs(x="Temporal Range", y="Observed Values")+theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25),legend.text=element_text(size=20),strip.text.x = element_text(size = 30),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),legend.title = element_text(hjust = 0.5,size = 20))# +coord_cartesian(xlim = c(-1,xmax))+theme_bw()+labs(x="Observed Value", y="Fraction of Data")+theme(axis.text.x = element_text(size=25),axis.text.y = element_text(size=25),legend.text=element_text(size=20),strip.text.x = element_text(size = 30),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),legend.title = element_text(hjust = 0.5,size = 20))+scale_color_discrete(guide = guide_legend(title = "Temporal Range",labels=c("Before 2010","After 2010"),keywidth=0.2,keyheight=0.4,default.unit = "inch"))
  return(list(cdf_plot_out,violin_plot_out,sample_1,sample_2,df_summary))
}

debug(plot_generator)
out<-plot_generator(analyte_name="Turbidity",xmax = 100)

out[[1]]

out[[2]]

max(out[[3]])
summary(out[[3]])
IQR<-summary(out[[3]])[3]-summary(out[[3]])[1]
outlier_cutoff<-summary(out[[3]])[3]+1.5*IQR

out[[3]][which(out[[3]]>outlier_cutoff)]

hist(out[[3]][which(out[[3]]>outlier_cutoff)],breaks = 100)

summary(out[[4]])

out<-plot_generator(analyte_name="Manganese",xmin=0,xmax = 2.5)
out[[2]]
out[[1]]
df_summary<-out[[5]]
write.table(df_summary, file="/Users/Amal/Box Sync/PSU/Summer 2018/Geoscience_Research/Imbalanced Project/results/real_data/v4_all/summary_stats_NWPA/Mn_summary.txt", row.names=FALSE, sep = ",")


####################################################
out<-plot_generator(analyte_name="Methane",xmin=-0.1,xmax = 2)

out[[1]]

summary(out[[3]])
summary(out[[4]])

