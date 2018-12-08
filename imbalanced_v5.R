#########################################################################################################
## Loading the required packages
library(perm)
library(lawstat)
library(foreach)
library(doParallel)
library(zoo)
library(data.table)
library(pastecs)
library(ggplot2)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#library(xts)

#########################################################################################################
#########################################################################################################
#########################################################################################################
## Tests for real data (NW PA revised)
## Reading and cleaning the old dataset
df_old<-read.csv(file = "/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/Data/NW_PA_Data/v4/Old_NW_PA_Amal.csv",stringsAsFactors = F)
str(df_old)

list_analytes_old<-list()
k<-1
for(j in seq(6,35,by = 2)){
  date_vec<-as.Date(gsub("(.*)/(..)$", "\\1/19\\2", df_old$Date_Sampled),format = "%m/%d/%Y",by = "month")
  
  ## Extracting non NA IDs as boolean vec
  date_IDs<-(!is.na(date_vec))
  sample_IDs<-(!is.na(df_old[,j]))
  common_IDs<-date_IDs&sample_IDs
  
  ## Creating dataframe for jth analyte
  df_temp_j<-data.frame("Date"=date_vec[common_IDs],"Analyte"=df_old[common_IDs,j])
  #names(df_temp_j)[2]<-colnames(df_old)[j]
  
  ## Considering only samples after 1980
  df_temp_j<-df_temp_j[which(df_temp_j$Date>as.Date("1/1/1980",format="%m/%d/%Y")),]
  
  #temp_j_mat<-as.matrix(df_temp_j$pH)
  #attributes(temp_j_mat)<-list("Date"=df_temp_j$Date)
  
  ## Converting the data frame into the time series object
  df_old_ts <- xts(x = df_temp_j$Analyte, order.by = df_temp_j$Date)
  
  #xts::plot.xts(x = df_old_ts,type = "p")
  df_temp_j_ordered <- df_temp_j[order(df_temp_j$Date),]
  #year_month_means <- setDT(df_temp_j_ordered)[, .(MontlyMeans = mean(pH)), by = .(year(Date), month(Date))]
  year_means_j <- setDT(df_temp_j_ordered)[, .(YearlyMeans = mean(Analyte)), by = .(year(Date))]
  
  library(pastecs)
  trend_test_j <- trend.test(ts(year_means_j$YearlyMeans), R=1)
  
  #xts(x = year_means_j$YearlyMeans,order.by = as.Date(x = year_means_j$year,format="%Y"))
  #plot(xts(x = year_means_j$YearlyMeans,order.by = as.Date(as.character(year_means_j$year),format="%Y")))
  # temp<-ts(data = df_temp_j$pH,start = min(df_temp_j$Date),end = max(df_temp_j$Date),names=df_temp_j$Date)
  # rownames(temp)<-matrix(df_temp_j$Date)
  list_analytes_old[[k]] <- list()
  
  list_analytes_old[[k]][[1]] <- year_means_j
  list_analytes_old[[k]][[2]] <- trend_test_j
  list_analytes_old[[k]][[3]] <- df_old_ts
  list_analytes_old[[k]][[4]] <- df_temp_j_ordered
  
  k<-k+1
}
names(list_analytes_old)<-colnames(df_old[,seq(6,35,by = 2)])

df_trend_summary<-data.frame("Analyte"=names(list_analytes_old),"Spearman rho"=rep(NA_real_,length(names(list_analytes_old))),"Spearman p value"=rep(NA_real_,length(names(list_analytes_old))))

for (k in 1: length(list_analytes_old)){
  df_trend_summary[k,2] <- list_analytes_old[[k]][[2]]$estimate
  df_trend_summary[k,3] <- list_analytes_old[[k]][[2]]$p.value
}

df_trend_summary_rounded<-data.frame("Analyte"=df_trend_summary[,1],round(x = as.matrix(df_trend_summary[,c(2,3)]),digits = 3))

df_trend_summary_rounded[which(df_trend_summary_rounded[,3]<0.05),3]<-paste0(df_trend_summary_rounded[which(df_trend_summary_rounded[,3]<0.05),3],"*")

#write.csv(df_trend_summary_rounded,file="/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/trend_summary_old.csv")

## Checking the plots, change k from 1 to 15 for different analytes.
plots_scatter_old<-list()
plots_year_old<-list()
df_scatter

log_indices<-c(6,7,13,14,15)

for (k in 1:15){
  temp<-as.data.frame(list_analytes_old[[k]][[1]])
  p<-ggplot(data = temp)
  plots_year_old[[k]] <- p+geom_line(aes(x = year, y = YearlyMeans))+labs(y = paste0("Yearly Averages for ",names(list_analytes_old)[k]), x = "Year") + theme_bw(base_size = 15)+coord_cartesian(xlim = c(1985,2000))+ theme(plot.title = element_blank(),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(face="bold",size=20),axis.text.y = element_text(size=20))+scale_x_continuous(minor_breaks = seq(from = 1985, to = 2000,length.out = 16))
  #plots_year_old[[k]]
  #, axis.title.x = element_blank(), axis.text.x = element_blank()
  #plots_scatter_old[[k]]<-xts::plot.xts(x = list_analytes_old[[k]][[3]],type = "p",main=paste0("Scatterplot for ",names(list_analytes_old)[k]))
  p<-ggplot(data = list_analytes_old[[k]][[4]])
  if(k %in% log_indices){
    plots_scatter_old[[k]] <- p+geom_point(aes(x = Date, y = log(Analyte)))+labs(y = paste0("Log of Scatterplot for ",names(list_analytes_old)[k]), x = "Year")+theme_bw(base_size = 15)+coord_cartesian(xlim = c(as.Date("1985-01-01", "%Y-%m-%d"),as.Date("2000-12-31", "%Y-%m-%d")))+ theme(plot.title = element_blank(),axis.title.x = element_blank(), axis.text.x = element_blank(),panel.grid.minor = element_line(),axis.title.y = element_text(face="bold",size=20),axis.text.y = element_text(size=20))+scale_x_date(minor_breaks = seq.Date(from = as.Date("1985-01-01", "%Y-%m-%d"), to = as.Date("2000-01-01", "%Y-%m-%d"),length.out = 16))
  }
  else{
    plots_scatter_old[[k]] <- p+geom_point(aes(x = Date, y = Analyte))+labs(y = paste0("Scatterplot for ",names(list_analytes_old)[k]), x = "Year")+theme_bw(base_size = 15)+coord_cartesian(xlim = c(as.Date("1985-01-01", "%Y-%m-%d"),as.Date("2000-12-31", "%Y-%m-%d")))+ theme(plot.title = element_blank(),axis.title.x = element_blank(), axis.text.x = element_blank(),panel.grid.minor = element_line(),axis.title.y = element_text(face="bold",size=20),axis.text.y = element_text(size=20))+scale_x_date(minor_breaks = seq.Date(from = as.Date("1985-01-01", "%Y-%m-%d"), to = as.Date("2000-01-01", "%Y-%m-%d"),length.out = 16))
  }
}
plots_year_old[[1]]
plots_scatter_old[[1]]


for (k in 1:15) {
  pdf(paste0("/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/plots/NW_PA_old/yearMeans/",names(list_analytes_old)[k],".pdf"))
  print(plots_year_old[[k]])
  dev.off()
}

for (k in 1:15) {
  pdf(paste0("/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/plots/NW_PA_old/scatter/",names(list_analytes_old)[k],".pdf"))
  print(plots_scatter_old[[k]])
  dev.off()
}

# #k<-7
# multiplot(plots_year_old[[1]],plots_year_old[[2]],plots_year_old[[3]],plots_year_old[[4]],plots_year_old[[5]],plots_year_old[[6]],plots_year_old[[7]],plots_year_old[[8]],plots_year_old[[9]],plots_year_old[[10]],plots_year_old[[11]],plots_year_old[[12]],plots_year_old[[13]],plots_year_old[[14]],plots_year_old[[15]],cols = 1)
# 
# multiplot(plots_scatter_old[[1]],plots_scatter_old[[2]],plots_scatter_old[[3]],plots_scatter_old[[4]],cols = 2)
# 

#########################################################################################################
## Reading and cleaning the new dataset
df_new<-read.csv(file = "/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/Data/NW_PA_Data/v4/New_NW PA_Amal.csv",stringsAsFactors = F)
str(df_new)

list_analytes_new<-list()
k<-1
for(j in seq(5,34,by = 2)){
  date_vec<-as.Date(gsub("(.*)/(..)$", "\\1/20\\2", df_new$Date.Sampled),format = "%m/%d/%Y",by = "month")
  
  ## Extracting non NA IDs as boolean vec
  date_IDs<-(!is.na(date_vec))
  sample_IDs<-(!is.na(df_new[,j]))
  common_IDs<-date_IDs&sample_IDs
  
  ## Creating dataframe for jth analyte
  df_temp_j<-data.frame("Date"=date_vec[common_IDs],"Analyte"=df_new[common_IDs,j])
  #names(df_temp_j)[2]<-colnames(df_old)[j]
  
  #temp_j_mat<-as.matrix(df_temp_j$pH)
  #attributes(temp_j_mat)<-list("Date"=df_temp_j$Date)
  
  ## Converting the data frame into the time series object
  df_new_ts <- xts(x = df_temp_j$pH, order.by = df_temp_j$Date)
  
  #xts::plot.xts(x = df_old_ts)
  df_temp_j_ordered <- df_temp_j[order(df_temp_j$Date),]
  #year_month_means <- setDT(df_temp_j_ordered)[, .(MontlyMeans = mean(pH)), by = .(year(Date), month(Date))]
  year_means_j <- setDT(df_temp_j_ordered)[, .(YearlyMeans = mean(Analyte)), by = .(year(Date))]
  
  trend_test_j <- trend.test(ts(year_means_j$YearlyMeans), R=1)
  
  #xts(x = year_means_j$YearlyMeans,order.by = as.Date(x = year_means_j$year,format="%Y"))
  #plot(xts(x = year_means_j$YearlyMeans,order.by = as.Date(as.character(year_means_j$year),format="%Y")))
  # temp<-ts(data = df_temp_j$pH,start = min(df_temp_j$Date),end = max(df_temp_j$Date),names=df_temp_j$Date)
  # rownames(temp)<-matrix(df_temp_j$Date)
  list_analytes_new[[k]] <- list()
  
  list_analytes_new[[k]][[1]] <- year_means_j
  list_analytes_new[[k]][[2]] <- trend_test_j
  list_analytes_new[[k]][[3]] <- df_new_ts
  list_analytes_new[[k]][[4]] <- df_temp_j_ordered
  k<-k+1
}
names(list_analytes_new)<-colnames(df_new[,seq(5,34,by = 2)])

df_trend_summary<-data.frame("Analyte"=names(list_analytes_new),"Spearman rho"=rep(NA_real_,length(names(list_analytes_new))),"Spearman p value"=rep(NA_real_,length(names(list_analytes_new))))

for (k in 1: length(list_analytes_new)){
  df_trend_summary[k,2] <- list_analytes_new[[k]][[2]]$estimate
  df_trend_summary[k,3] <- list_analytes_new[[k]][[2]]$p.value
}
df_trend_summary

df_trend_summary_rounded<-data.frame("Analyte"=df_trend_summary[,1],round(x = as.matrix(df_trend_summary[,c(2,3)]),digits = 3))

df_trend_summary_rounded[which(df_trend_summary_rounded[,3]<0.05),3]<-paste0(df_trend_summary_rounded[which(df_trend_summary_rounded[,3]<0.05),3],"*")

#write.csv(df_trend_summary_rounded,file="/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/trend_summary_new.csv")

## Checking the plots, change k from 1 to 15 for different analytes.
plots_scatter_new<-list()
plots_year_new<-list()

log_indices<-c(3,5,6,14)

for (k in 1:15){
  temp<-as.data.frame(list_analytes_new[[k]][[1]])
  p<-ggplot(data = temp)
  plots_year_new[[k]] <- p+geom_line(aes(x = year, y = YearlyMeans))+labs(y = paste0("Yearly Averages for ",names(list_analytes_new)[k]), x = "Year") + theme_bw(base_size = 15)+coord_cartesian(xlim = c(2012,2016))+ theme(plot.title = element_blank(),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(face="bold",size=20),axis.text.y = element_text(size=20))+scale_x_continuous(minor_breaks = seq(from = 2012, to = 2016,length.out = 5))
  #plots_year_new[[k]]
  #, axis.title.x = element_blank(), axis.text.x = element_blank()
  #plots_scatter_new[[k]]<-xts::plot.xts(x = list_analytes_new[[k]][[3]],type = "p",main=paste0("Scatterplot for ",names(list_analytes_new)[k]))
  p<-ggplot(data = list_analytes_new[[k]][[4]])
  if(k %in% log_indices){
    plots_scatter_new[[k]] <- p+geom_point(aes(x = Date, y = log(Analyte)))+labs(y = paste0("Log of Scatterplot for ",names(list_analytes_new)[k]), x = "Year")+theme_bw(base_size = 15)+coord_cartesian(xlim = c(as.Date("2012-01-01", "%Y-%m-%d"),as.Date("2015-12-31", "%Y-%m-%d")))+ theme(plot.title = element_blank(),axis.title.x = element_blank(), axis.text.x = element_blank(),panel.grid.minor = element_line(),axis.title.y = element_text(face="bold",size=20),axis.text.y = element_text(size=20))+scale_x_date(minor_breaks = seq.Date(from = as.Date("2012-01-01", "%Y-%m-%d"), to = as.Date("2016-12-31", "%Y-%m-%d"),length.out = 5))
  }
  else{
    plots_scatter_new[[k]] <- p+geom_point(aes(x = Date, y = Analyte))+labs(y = paste0("Scatterplot for ",names(list_analytes_new)[k]), x = "Year")+theme_bw(base_size = 15)+coord_cartesian(xlim = c(as.Date("2012-01-01", "%Y-%m-%d"),as.Date("2015-12-31", "%Y-%m-%d")))+ theme(plot.title = element_blank(),panel.grid.minor = element_line(),axis.title.x = element_blank(), axis.text.x = element_blank(),axis.title.y = element_text(face="bold",size=20),axis.text.y = element_text(size=20))+scale_x_date(minor_breaks = seq.Date(from = as.Date("2012-01-01", "%Y-%m-%d"), to = as.Date("2016-01-01", "%Y-%m-%d"),length.out = 5))
  }
}
plots_year_new[[1]]
plots_scatter_new[[1]]


for (k in 1:15) {
  pdf(paste0("/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/plots/NW_PA_new/yearMeans/",names(list_analytes_new)[k],".pdf"))
  print(plots_year_new[[k]])
  dev.off()
}

for (k in 1:15) {
  pdf(paste0("/Users/Amal/Box Sync/PSU/Fall 2018/Geoscience_Research/Imbalanced Project/plots/NW_PA_new/scatter/",names(list_analytes_new)[k],".pdf"))
  print(plots_scatter_new[[k]])
  dev.off()
}

#########################################################################################################
