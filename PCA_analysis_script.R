#install.packages("car")
library("car")
library("ggplot2")
library("mvoutlier")
library("MASS")
library("mvnormtest")
library("pastecs")
library("reshape")
library("WRS") # install.packages("WRS", repos="http://R-Forge.R-project.org")
library("zoo") ### for mean imputation
library("corpcor")
library("GPArotation")
library("psych")
library("corrplot")
library("gmodels")

rm(list = ls())

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

curr_path <- "C:/Users/Ilona/Dropbox/Postdoc/Papers/DRTT/Stats"
tractlist <- c("L_DRT_chopped", "R_DRT_chopped", "L_SPCT_chopped", "R_SPCT_chopped", "L_MCP", "L_ICP", "R_ICP", "L_SCP", "R_SCP")
tractlist_clean <- c("L_DRT", "R_DRT", "L_SPCT", "R_SPCT", "L_MCP", "L_ICP", "R_ICP", "L_SCP", "R_SCP")

rm("all_values")

for (tract in tractlist){

  print(tract)

  sprtshtname <- paste(curr_path,"/",tract,".csv",sep="")
  CM <- read.csv(sprtshtname,sep=";")
  idx <- CM$GROUP < 3 & is.na(CM$AGE) == FALSE
  CM <- CM[idx,c(2:6,13)]

  if ( exists("all_values") == FALSE ) {
    all_values <- data.frame(CM[,c(1,2,3)]) ### create new dataframe for first occurance
    for (col in c(4,5,6)){
      columnname <- paste(strsplit(tract,"_chopped")[1],"_",names(CM[col]),sep="")
      all_values[,columnname] <- CM[,col]
    }
  } else {
    CM <- CM[,c(4,5,6)]
    for (col in c(1,2,3)){
      columnname <- paste(strsplit(tract,"_chopped")[1],"_",names(CM[col]),sep="")
      print(columnname)
      all_values[,columnname] <- CM[,col]
    }
  }
}

### now we can check frequency of tract reconstruction in SPCT and DRT
exist_L_DRT <- is.na(all_values$L_DRT_chopped_FA) == FALSE
exist_R_DRT <- is.na(all_values$R_DRT_chopped_FA) == FALSE
exist_L_SPCT <- is.na(all_values$L_SPCT_chopped_FA) == FALSE
exist_R_SPCT <- is.na(all_values$R_SPCT_chopped_FA) == FALSE
CrossTable(exist_L_DRT,exist_R_DRT, chisq = TRUE)
CrossTable(exist_L_SPCT,exist_R_SPCT, chisq = TRUE)
CrossTable(exist_L_SPCT,exist_L_DRT, chisq = TRUE)
CrossTable(exist_R_SPCT,exist_R_DRT, chisq = TRUE)

### PCA
all_fa <- stack(all_values[, seq(from = 4, to = ncol(all_values), by = 3)])
all_rd <- stack(all_values[, seq(from = 5, to = ncol(all_values), by = 3)])
all_hmoa <- stack(all_values[, seq(from = 6, to = ncol(all_values), by = 3)])
maxime_style_df <- data.frame(all_fa[,1])
maxime_style_df$all_rd <- all_rd[,1]
maxime_style_df$all_hmoa <- all_hmoa[,1]
maxime_style_df$bundles <- rep(seq(from = 1, to = length(tractlist), by = 1), each = nrow(all_values)) ### column for which tract it is
maxime_style_df$subject <- rep(seq(from = 1, to = nrow(all_values), by = 1), times = length(tractlist)) ### column for which subject it is

# METRIC PCA
for_pca_matrix = na.aggregate(maxime_style_df[,1:3]) ### only use microstructural metrics
### pre-PCA
print("CORRELATION MATRIX")
raqMatrix <- cor(for_pca_matrix) ### correlation matrix
print("BARTLETT")
print(cortest.bartlett(raqMatrix,n=nrow(all_values))) ### needs to be significant
print("IS SUITABLE")
det(raqMatrix) > .00001 ### check this test
### PCA, seems to automatically standardize values
pcModel <- principal(for_pca_matrix, nfactors=nrow(raqMatrix), rotate="none")
print(pcModel) ## Proportion Var for variance xplained, SSloadings are eigenvalues
### make screeplot
plot(pcModel$values, type="b") ### screeplot

# extract first component scores for each tract
maxime_style_df$all_scores <- pcModel$scores[,1]
comp_per_tract <- unstack(maxime_style_df, all_scores ~ bundles)
names(comp_per_tract) <- tractlist_clean

### tract correlations
t_cor <- cor(comp_per_tract)
t_cor_sorted <- mat.sort(t_cor)

### group comparison
figname <- paste(curr_path,"/Figures/Boxplot_component_score.jpg",sep="")
jpeg(figname, width=800, height=800, bg="white")
par(mfrow=c(3,3))
for (t in 1:ncol(comp_per_tract)){
  boxplot(comp_per_tract[,t]~all_values$GROUP,xlab="Group",ylab=tractlist[t],names=c('PD','HC'))
}
dev.off()


corrplot(t_cor_sorted, method="color")
figname <- paste(curr_path,"/Figures/Tract_correlation.jpg",sep="")
jpeg(figname, width=1000, height=1000, bg="white")
corrplot(t_cor, method="color")
dev.off()
### tract PCA
for_pca_matrix2 = comp_per_tract ### only use microstructural metrics
### pre-PCA
print("CORRELATION MATRIX")
raqMatrix2 <- cor(for_pca_matrix2)
print("BARTLETT")
print(cortest.bartlett(raqMatrix2,n=nrow(comp_per_tract))) ### needs to be significant
print("IS SUITABLE")
det(raqMatrix) > .00001 ### check this test
### PCA, seems to automatically standardize values
pcModel2 <- principal(for_pca_matrix2, nfactors=nrow(raqMatrix2), rotate="none")
print(pcModel2)
### make screeplot
figname <- paste(curr_path,"/Figures/Tract_PCA_screeplot.jpg",sep="")
jpeg(figname, width=1000, height=1000, bg="white")
plot(pcModel2$values, type="b",xlab="Component", ylab="Eigenvalue") ### screeplot
dev.off()

overall_val <- pcModel2$scores[,1]
boxplot(overall_val~all_values$GROUP,xlab="Group",ylab=tractlist[t],names=c('PD','HC'))
model <- t.test(overall_val ~ all_values$GROUP, paired = FALSE)
print(model)
for (c in 1:3){
  print(colnames(pcModel2$scores)[c])
  model <- t.test(pcModel2$scores[,c]~ all_values$GROUP, paired = FALSE)
  print(model)
}


## VARIMAX ROTATION
pcModel3 <- principal(for_pca_matrix2, nfactors=3, rotate="varimax")
print(pcModel3)
for (c in 1:3){
 print(colnames(pcModel3$scores)[c])
   model <- t.test(pcModel3$scores[,c]~ all_values$GROUP, paired = FALSE)
  print(model)
}
figname <- paste(curr_path,"/Figures/Boxplot_component_score_for_tract_PCA_varimax.jpg",sep="")
jpeg(figname, width=800, height=800, bg="white")
par(mfrow=c(1,3))
for (c in 1:3){
  boxplot(pcModel3$scores[,c]~all_values$GROUP,xlab="Group",ylab=paste("Component score ",c,sep=""),names=c('PD','HC'))
}
dev.off()

## OBLIQUE ROTATION
pcModel4 <- principal(for_pca_matrix2, nfactors=3, rotate="oblimin")
print(pcModel4)
for (c in 1:3){
  print(colnames(pcModel4$scores)[c])
  model <- t.test(pcModel4$scores[,c]~ all_values$GROUP, paired = FALSE)
  print(model)
}
figname <- paste(curr_path,"/Figures/Boxplot_component_score_for_tract_PCA_oblique.jpg",sep="")
jpeg(figname, width=800, height=800, bg="white")
par(mfrow=c(1,3))
for (c in 1:3){
  boxplot(pcModel4$scores[,c]~all_values$GROUP,xlab="Group",ylab=paste("Component score ",c,sep=""),names=c('PD','HC'))
}
dev.off()
