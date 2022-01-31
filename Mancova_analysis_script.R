#install.packages("car")
library("car")
library("ggplot2")
library("mvoutlier")
library("MASS")
library("mvnormtest")
library("pastecs")
library("reshape")
library("WRS") # install.packages("WRS", repos="http://R-Forge.R-project.org")

rm(list = ls())

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

curr_path <- "C:/Users/Ilona/Dropbox/Postdoc/Papers/DRTT/Stats"
tractlist <- c("L_DRT_chopped", "R_DRT_chopped", "L_SPCT_chopped", "R_SPCT_chopped", "L_MCP", "L_ICP", "R_ICP", "L_SCP", "R_SCP")

for (tract in tractlist){

  print(tract)
  sprtshtname <- paste(curr_path,"/",tract,".csv",sep="")
  CM <- read.csv(sprtshtname,sep=";")
  ### how is ("Cerebellar_mancova.csv") set up?
  
  idx <-CM$GROUP < 3 & CM$FA > 0
  idx[ is.na(idx) ] <- FALSE
  CM <- CM[idx,2:(ncol(CM)-1)]  ### group 3 would be HC60, we don't need them
  
  CM$GROUP<-as.factor(CM$GROUP)
  CM$SEX<-as.factor(CM$SEX)
  
  CM <- na.omit(CM) ### get rid of NAs
  
  ### explore multivariance
  Y<-CM[,c(4:5,12)]
  figname <- paste(curr_path,"/Figures/DV_scatters_",tract,".jpg",sep="")
  jpeg(figname, width=800, height=300, bg="white")
  pairs(Y)
  dev.off()
  
  ### descriptives
  measures <- c('FA','MD','RD','AD','HMOA')
  
  Y_PD <- Y[CM$GROUP == 1,]
  Y_HC<- Y[CM$GROUP == 2,]
  
  ### normality
  print("SHAPIRO PD")
  print(mshapiro.test(t(Y_PD)))
  print("SHAPIRO HC")
  print(mshapiro.test(t(Y_HC)))
  
  ### outliers
  figname <- paste(curr_path,"/Figures/Outliers_",tract,".jpg",sep="")
  jpeg(figname, width=800, height=300, bg="white")
  aq.plot(Y)
  dev.off()
  
  ### robust manova
  print("robust manova")
  CM$row <- c(1:sum(CM$GROUP==1), 1:sum(CM$GROUP==2))
  datamelt <- melt(CM[,c(1,13,4,5,12)], id = c("GROUP", "row"), measured = c("FA", "MD", "RD", "AD", "HMOA"))
  names(datamelt) <- c("Group", "Row", "Outcome", "Frequency")
  datacast <- cast(datamelt, Row ~ Group + Outcome, value = "Frequency")
  datacast$Row <- NULL
  print(mulrank(2,3,datacast))
  
  ### make boxplots
  figname <- paste(curr_path,"/Figures/Histograms_",tract,".jpg",sep="")
  jpeg(figname, width=800, height=300, bg="white")
  par(mfrow=c(1,5))
  measures <- c('FA','RD','HMOA')
  for (m in measures){
    boxplot(CM[,m]~CM$GROUP,xlab="Group",ylab=m,names=c('PD','HC'))
  }
  dev.off()

}

ps <- c(0.1888, 0.2497, 0.0703, 0.9096, 0.1777, 0.1005, 0.0205, 0.0973, 0.0141)
p.adjust(ps, method ="fdr", n = length(ps))

### descriptive stats
model <- t.test(CM$AGE~ CM$GROUP, paired = FALSE)
print(model)
print(by(CM$AGE, CM$GROUP, stat.desc, basic = FALSE))
CrossTable(CM$SEX, CM$GROUP, chisq = TRUE)