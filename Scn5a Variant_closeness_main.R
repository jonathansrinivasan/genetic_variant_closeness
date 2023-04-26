---
title: Predict LQTS Diagnosis Probability Using Structure, Function, and *In Silico*
author: "Laura Bertolami, Shoshana Kelner, and Brett Kroncke"
date: "Fall 2021"
output:
  html_document:
  theme: flatly
toc: yes
toc_depth: 4
toc_float: no
smooth_scroll: yes
code_folding: hide
highlight: zenburn
pdf_document:
  toc: yes
toc_depth: '4'
word_document:
  toc: yes
toc_depth: '4'
editor_options:
  chunk_output_type: console
---
  
```{r preamble,include=FALSE}

knitr::opts_chunk$set(echo = TRUE)
library("DBI")
library("RSQLite")
library(dplyr)
library(ggplot2)
library(ggpubr)
library(caret)
library(plotrix)
library(glmnet)
library(meta)
library(reshape2)
library(psych)
require(Hmisc)
library(tableone)
library(rms)
library(boot)
library(leaps)
library(car)
library(reticulate)
library(rootSolve)
library(pROC)
library(wCorr) 
library(MALDIquant)
library(tidyverse)      # data manipulation and visualization
library(lubridate)      # easily work with dates and times
library(fpp2)           # working with time series data
library(zoo)            # working with time series data
library(latex2exp)
library(forestplot)
library(ggplot2)
library(readxl)
library(PRROC)
library(rmarkdown)
library(lattice)

source('func_dist_seq.R')

#source('src/func_dist_seq.R') necessary? 
scn5adist<-read.csv(file = "data/Covariates/structure_files/scn5a_distances_final.csv")
scn5a.data<-read.csv(file = "scn5a_all_vars_annotated.csv")
cln.d<-scn5a.data

```

# Introduction

# Part 1: Calculate probability of LQT2 diagnosis and LQT2 Probability Density using Various Subsets of the Literature and Cohort Data

## Literature and Cohort Combined (for final predictions)

This only needs to be run if you are not loading an RData that is already clean.

```{r}

 #cln.d <- cln.d[cln.d$mut_type == "missense" & cln.d$isoform=="A",]
 cln.d[is.na(cln.d$total_carriers),"total_carriers"] <- 0
#d[is.na(d$lqt2),"lqt2"] <- 0
 # set initial weighting and penetrance
 cln.d$weight = 1-1/(0.01+cln.d$total_carriers)
 cln.d$brs1_penetrance <- cln.d$brs1/cln.d$total_carriers
 cln.d[cln.d$total_carriers < 1,"weight"] <- 0.000 # This is changed to "< 2" here to evaluate ROC-AUC of n=1 variants from the literature


```

### LQT2 empirical diagnosis probability prior

Use observed LQT2 diagnosis probability to calculate "LQTS probability density" as described in previous publication. Plot diagnosis probability density versus residue

This only needs to be run if you are not loading an RData that is already clean.

```{r}
# Mean squared error
# 
# mse <- function(sm) {
#   mean((sm$residuals)^2*(sm$weights))
# }
# # Derive alpha and beta from weighted mean and MSE (estimated variance)
# estBetaParams <- function(mu, var) {
#   alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
#   beta <- alpha * (1 / mu - 1)
#   return(params = list(alpha = alpha, beta = beta))
# }
# # Weighted mean to determine LQT2 penetrance empirical prior
# newdata = data.frame(wt=1)
# model <- lm(brs1_penetrance ~ 1, data=cln.d, weights = cln.d$weight)
# 
# summary(model)
# p<-predict(model, newdata)
# dev<- mse(model) #p*(1-p)
# # Estimated shape parameters for LQT2 empirical prior
# alpha0 = estBetaParams(p,dev)$alpha
# beta0 = estBetaParams(p,dev)$beta
# print(paste("alpha0 = ", alpha0, "  beta0 = ", beta0))
# # Bayesian LQT2 penetrance estimates from empirical priors
# # and observed affected/unaffected counts:
# cln.d$brs1_penetranceBayesian_initial <- (alpha0 + cln.d[,"brs1"])/((alpha0 + beta0 + cln.d[,"total_carriers"]))
# cln.d$brs1_penetranceBayesian<-cln.d$brs1_penetranceBayesian_initial
# 
# 
# combined.data<- cln.d
#All data is clean in Updated.Data, because the cleaned data was saved as such
# this chunk is redundant for us but necessary if new data need to run through chunk 13

```

```{r, echo=FALSE, warning=FALSE, message=FALSE}

# rm(cln.d)
# cln.d<-herg.combined.data[herg.combined.data$total_carriers>0,]
# cln.d$penetrance_lqt2<-cln.d$lqt2/cln.d$total_carriers
# 
# cln.d[, "lqt2_dist"]<-NA
# cln.d[, "lqt2_dist_weight"]<-NA
# 
# cln.d <- cln.d[cln.d$mut_type == "missense" & cln.d$isoform == "A",]
# cln.d<- cln.d[!is.na(cln.d$resnum),]
# 
# cln.d$var<- as.character(cln.d$var)
# 
# #instead of funcdist3 when looking for only observed data
# for (rec in 1:nrow(cln.d)){
#   newdata<-cln.d[is.na(match(cln.d$resnum, h2dist$V1)),]
# }
# for(rec in 1:nrow(cln.d)){
#   cln.d<-cln.d[is.na(match(cln.d$resnum, newdata$resnum)),]
# }
# 
# for(rec in 1:nrow(cln.d)){
#   cln.d[rec,
#         c("lqt2_dist", "lqt2_dist_weight", "resnum")] <- c(funcdist2(cln.d[rec, "resnum"], cln.d[rec, "var"], cln.d, cln.d$h2dist, "penetrance_lqt2", "sigmoid", 7,'V8'), cln.d[rec, "resnum"])
# }
# 
# save(cln.d, file="Seq.Obs.Data.Test.Rdata")

# #running through code for brs1
# # 
# cln.d<-cln.d[cln.d$total_carriers>0,]
# cln.d$brs1_penetrance<-cln.d$brs1/cln.d$total_carriers
# 
# cln.d[, "brs1_dist"]<-NA
# cln.d[, "brs1_dist_weight"]<-NA
# 
# #cln.d <- cln.d[cln.d$mut_type == "missense" & cln.d$isoform == "A",] ASK
# cln.d<- cln.d[!is.na(cln.d$resnum),]
# 
# cln.d$var<- as.character(cln.d$var)
# # 
# # #instead of funcdist3 when looking for only observed data
# # for (rec in 1:nrow(cln.d)){
# # newdata<-cln.d[is.na(match(cln.d$resnum, scn5adist$resnum)),]
# # }
# # for(rec in 1:nrow(cln.d)){
# #  cln.d<-cln.d[is.na(match(cln.d$resnum, newdata$resnum)),]
# # }
# # 
# for(rec in 1:nrow(cln.d)){
#   cln.d[rec,
#         c("brs1_dist", "brs1_dist_weight", "resnum")] <- c(funcdist2(cln.d[rec, "resnum"], cln.d[rec, "var"], cln.d, cln.d$scn5adist, "brs1_penetrance", "sigmoid", 7,'distance'), cln.d[rec, "resnum"])
# }
# cln.d<- cln.d[!is.na(cln.d$brs1_dist),]
# save(cln.d, file="scn5a_brs1_seq.data.Rdata")

#scn5a_brs1_all.data.rdata is everything; scn5a_brs1_cor.data.rdata is observed and scn5a_brs1_seq.data.rdata is the sequenced data 


```

### Calculate LQTS probability densities and annotate function and structural location

With the updated empirical priors applied to carrier counts, calculate "LQTS probability density" as described in previous publication. !!! NOTE: since these data are truly the "best estimates" we include all variants in the calculation such that unique scores are by residue not by variant.

# Part 2: Variance explained


```{r}


## Pearson R\^2 and Spearman Rho Against EM Posterior from Cohort


#1 

calcPval=function(xName,yName,weightName,nPerms,new.mat2){
  # Pulls out variables
  x=new.mat2[,xName] 
  y=new.mat2[,yName] 
  w=new.mat2[,weightName]
  x2=x[!is.na(x)]
  y2=y[!is.na(x)]
  w2=w[!is.na(x)]
  # Calculate the real correlation
  realCorr=weightedCorr(x2,y2,method='spearman',weights=w2)
  # Do permutations, calculate fake correlations
  permutedCorrList=c()
  for(permNum in 1:nPerms){
    permutedX=sample(x2,length(x2),replace=FALSE)
    wCorrSim=weightedCorr(permutedX,y2,method='spearman',weights=w2)
    permutedCorrList=c(permutedCorrList,wCorrSim)
  }
  permutedCorrList2=abs(permutedCorrList)
  realCorr2=abs(realCorr)
  
  # Calculate pvalue
  summ=sum(realCorr2<permutedCorrList2)
  pValue=summ/nPerms
  return(list(realCorr,pValue,length(x2)))
}


#2


calcAllPvals=function(yList,xList,nPerms,weightName,new.mat2){
  i=0
  resultTable=data.frame()
  for(yName in yList){
    for(xName in xList){
      i=i+1
      result=calcPval(xName,yName,weightName, nPerms, new.mat2)
      resultTable[i,'x']=xName
      resultTable[i,'y']=yName
      resultTable[i,'nPerms']=nPerms
      resultTable[i,'weightedCorr']=result[[1]]
      resultTable[i,'pValue']=result[[2]]
      resultTable[i,'n']=result[[3]]
      #print(resultTable[i,'pValue'])
    }
  }
  print(resultTable)
  return(resultTable)
}

```

#Forest Plot Data Manipulation (add lqt2_dist from all sets)


```{r, echo=FALSE, warning=FALSE, message=FALSE, include=FALSE}
# # # run this code once you get all of the brs1
# load("scn5a_brs1_all.data.Rdata")
# scn5a.cln.d<- cln.d
# load("scn5a_brs1_seq.data.Rdata")
# distseq<-cln.d
# load("scn5a_brs1_cor.data.Rdata")
# distcor<-cln.d
# 
# distseq<-select(distseq,var, resnum, brs1_dist, brs1_dist_weight)
# distcor<-select(distcor,var, resnum, brs1_dist, brs1_dist_weight)
# 
# # 
# scn5a.cln.d<-left_join(scn5a.cln.d, distseq, by=c('resnum','var'), suffix= c('','_sequence'))
# scn5a.cln.d<-left_join(scn5a.cln.d, distcor, by=c('resnum', 'var'),  suffix= c('','_correlated'))
# # ## herg.combined.data<-left_join(herg.combined.data, maxall, by=c('resnum','var'), suffix= c('','_max'))
# # ## herg.combined.data<-left_join(herg.combined.data, meanall, by=c('resnum','var'), suffix= c('','_mean'))
# # ## herg.combined.data<-left_join(herg.combined.data, maxobs, by=c('resnum', 'var'),  suffix= c('','_obs_max'))
# # ## herg.combined.data<-left_join(herg.combined.data, meanobs, by=c('resnum', 'var'),  suffix= c('','_obs_mean'))
# # ## herg.combined.data<-left_join(herg.combined.data, seqobs, by=c('resnum', 'var'),  suffix= c('','_obs_seq'))
# # 
# # #``{r}
# # #knitr::opts_chunk$set(fig.height = 9, fig.width = 7)
# # # Select covariates from isoform "A" in cohort dataset
# # #
#  save(scn5a.cln.d, file="all.brs1.data.trunc.Rdata")
```


# Part #3: Forest plot
##Forest Plot of all data 
 

```{r}
load("all.brs1.data.trunc.Rdata")
scn5a.cln.d <- scn5a.cln.d[,!names(scn5a.cln.d) %in%  "cardiacboost"]
scn5a.cln.d <- unique(scn5a.cln.d)
results.1.data<-scn5a.cln.d[!is.na(scn5a.cln.d$provean_score) & !is.na(scn5a.cln.d$REVEL) & !is.na(scn5a.cln.d$pamscore),]
#results.1.data<- results.1.data(results.1.data$total_carriers>0, )

brs1_penetrance <- na.omit(scn5a.cln.d$brs1_penetrance)

yList=c("brs1_penetrance")

#need to add the distance values to the scn5a.data one through the nasty file thats saved above 
xList=c('pph2_prob', 'provean_score', "blastpssm",
        'pamscore', "brs1_dist",'brs1_dist_sequence', 'brs1_dist_correlated', "REVEL")

weight <- c(scn5a.cln.d$weight )
weight <- weight[weight!=0]

resultTable<-calcAllPvals(yList, xList, 1000, 'weight', results.1.data[results.1.data$total_carriers>0,])
rm(results.1.data)

#4

i=0
FP.data<-data.frame()
for (x in xList){
  i=i+2
  FP.data[i,"Feature"]<-x
  t<-scn5a.cln.d[!is.na(scn5a.cln.d[,x]) & scn5a.cln.d$total_carriers>0,]
  
  #<-t[!is.na(t[,"var"]),]
  #FP.data[i,"Feature"]<-paste(x,"_cohort")
  #t<-herg.combined.data[!is.na(herg.combined.data[,x]) & herg.combined.data$total_carriers>0,]
  
  t<-t[!is.na(t[,"var"]),]
  foo <- boot(t, function(data,indices)
    weightedCorr(t[indices,x],t$brs1_penetrance[indices], method="Spearman", weights = t$weight[indices]), R=1000)
  
  
  FP.data[i,"Spearman"]<-foo$t0
  FP.data[i,"Spearman_low"]<-quantile(foo$t,c(0.025,0.975), na.rm = T)[1][[1]]
  FP.data[i,"Spearman_high"]<-quantile(foo$t,c(0.025,0.975), na.rm = T)[2][[1]]
  FP.data[i,"n"]<-length(t[,x]) 
}

#5

forestplot(FP.data$Feature,FP.data$Spearman,FP.data$Spearman_low,FP.data$Spearman_high, cex = 0.5, boxsize = 0.8)

rm(FP.data)

tmp<-scn5a.cln.d[!is.na(scn5a.cln.d$provean_score) & !is.na(scn5a.cln.d$REVEL) & !is.na(scn5a.cln.d$pamscore),]


# # Weighted R2 between observed LQT2 penetrance and post-test probability
##Idk if we need this, will look for p_mean_w to add on later -L
# foo <- boot(tmp, function(data,indices)
#   weightedCorr(tmp$p_mean_w[indices],tmp$brs1_penetrance[indices], method="pearson", weights = tmp$weight[indices])^2, R=1000)
# print("EM estimated BRS1 diagnosis probability versus observed cohort BRS1 diagnosis probability")
# foo$t0
# quantile(foo$t,c(0.025,0.975), na.rm = T)

brs1_dist <- na.omit(scn5a.cln.d$brs1_dist)

model<- lm(brs1_penetrance~brs1_dist, data = tmp, weights = weight)
quantile(foo$t,c(0.025,0.975), na.rm = TRUE) #added na.rm =TRUE
model <- lm(brs1_penetrance~brs1_dist, data = tmp, weights = weight)
mod<-data.frame(model$fitted.values,model$model$brs1_penetrance,model$weights)
foo <- boot(mod, function(data,indices)
  weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
print("BRS1 diagnosis probability density versus observed cohort BRS1 diagnosis probability")
foo$t0
quantile(foo$t,c(0.025,0.975))
model <- lm(brs1_penetrance~REVEL, data = tmp, weights = weight) 
mod<-data.frame(model$fitted.values,model$model$brs1_penetrance,model$weights)
foo <- boot(mod, function(data,indices)
  weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
print("REVEL score versus observed cohort BRS1 diagnosis probability")
foo$t0
quantile(foo$t,c(0.025,0.975))
rm(tmp)

# #6 up to here with renaming 
# 
# # Evaluate only variants with Heterozygous peak tail current measured.
# ##might comment this section out until i find p_mean_w -L IT CRASHES RSTUDIO MAKE SURE EVERYTHINGS SAVED
# results.2.data<-scn5a.cln.d[!is.na(scn5a.cln.d$ht_tailPeak),]
# 
# resultTable<-calcAllPvals(yList, xList, 1000, 'weight', results.2.data[results.2.data$total_carriers>0,])
# 
# foo <- boot(results.2.data, function(data,indices)
#   weightedCorr(results.2.data$p_mean_w[indices],results.2.data$brs1_penetrance[indices], method="pearson", weights = results.2.data$weight[indices])^2, R=1000)
# print("EM estimated LQT2 diagnosis probability versus observed cohort LQT2 diagnosis probability")
# foo$t0
# 
# quantile(foo$t,c(0.025,0.975), na.rm = T)
# 
# model <- lm(brs1_penetrance~ht_tailPeak, data = results.2.data, weights = weight)
# mod<-data.frame(model$fitted.values,model$model$pbrs1_penetrance,model$weights)
# foo <- boot(mod, function(data,indices)
#   weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
# print("Heterozygously measured peak tail current versus observed cohort LQT2 diagnosis probability")
# foo$t0
# quantile(foo$t,c(0.025,0.975))
# model<-lm(brs1_penetrance~brs1_dist, data = results.2.data, weights = weight)
# mod<-data.frame(model$fitted.values,model$model$brs1_penetrance,model$weights)
# foo <- boot(mod, function(data,indices)
#   weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
# print("LQT2 probability density versus observed cohort LQT2 diagnosis probability")
# foo$t0
# quantile(foo$t,c(0.025,0.975))
# model <- lm(brs1_penetrance~REVEL, data = results.2.data, weights = weight)
# mod<-data.frame(model$fitted.values,model$model$brs1_penetrance,model$weights)
# foo <- boot(mod, function(data,indices)
#   weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
# print("REVEL score versus observed cohort LQT2 diagnosis probability")
# foo$t0
# quantile(foo$t,c(0.025,0.975))
# rm(results.2.data)
```




## Variance explained from literature dataset

```{r}

tmp2<-scn5a.cln.d[!is.na(scn5a.cln.d$provean_score) & !is.na(scn5a.cln.d$brs1_penetrance) & !is.na(scn5a.cln.d$brs1_dist) & !is.na(scn5a.cln.d$REVEL),]
foo <- boot(tmp2, function(data,indices)
  weightedCorr(tmp2$p_mean_w[indices],tmp2$brs1_penetrance[indices], method="pearson", weights = tmp2$weight[indices])^2, R=1000)
print("EM estimated LQT2 diagnosis probability versus observed literature BRS1 diagnosis probability")
foo$t0
quantile(foo$t,c(0.025,0.975))
model <- lm(brs1_penetrance~brs1_dist, data = tmp2, weights = weight)
mod<-data.frame(model$fitted.values,model$model$brs1_penetrance,model$weights)
foo <- boot(mod, function(data,indices)
  weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
print("BRS1 probability density versus observed literature BRS1 diagnosis probability")
foo$t0
quantile(foo$t,c(0.025,0.975))
model <- lm(brs1_penetrance~REVEL, data = tmp2, weights = weight)
mod<-data.frame(model$fitted.values,model$model$brs1_penetrance,model$weights)
foo <- boot(mod, function(data,indices)
  weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
print("REVEL versus observed literature BRS1 diagnosis probability")
foo$t0
quantile(foo$t,c(0.025,0.975))
# # Evaluate only variants with Heterozygous peak tail current measured.
# #ht_tailPeak is not a variable in scn5a data so far and p_mean_w is here as well  -L
# tmp3<-scn5a.cln.d[!is.na(scn5a.cln.d$ht_tailPeak) & !is.na(scn5a.cln.d$provean_score) & !is.na(scn5a.cln.d$brs1_penetrance) & !is.nan(scn5a.cln.d$brs1_penetrance) & !is.na(scn5a.cln.d$brs1_dist),]
# foo <- boot(tmp3, function(data,indices)
#   weightedCorr(tmp3$p_mean_w[indices],tmp3$brs1_penetrance[indices], method="pearson", weights = tmp3$weight[indices])^2, R=1000)
# print("EM estimated BRS1 diagnosis probability versus observed literature BRS1 diagnosis probability")
# foo$t0
# quantile(foo$t,c(0.025,0.975))
# brs1_dist <- 
#   model <- lm(brs1_penetrance~brs1_dist, data = tmp3, weights = weight)
# mod<-data.frame(model$fitted.values,model$model$brs1_penetrance,model$weights)
# foo <- boot(mod, function(data,indices)
#   weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
# print("BRS1 probability density versus observed literature BRS1 diagnosis probability")
# foo$t0
# quantile(foo$t,c(0.025,0.975))
# model <- lm(brs1_penetrance~ht_tailPeak, data = tmp3, weights = weight)
# mod<-data.frame(model$fitted.values,model$model$brs1_penetrance,model$weights)
# foo <- boot(mod, function(data,indices)
#   weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
# print("Heterozygous peak tail current versus observed literature BRS1 diagnosis probability")
# foo$t0
# quantile(foo$t,c(0.025,0.975))
# model <- lm(brs1_penetrance~REVEL, data = tmp3, weights = weight)
# mod<-data.frame(model$fitted.values,model$model$brs1_penetrance,model$weights)
# foo <- boot(mod, function(data,indices)
#   weightedCorr(mod$model.fitted.values[indices],mod$model.model.brs1_penetrance[indices], method="pearson", weights = mod$model.weights[indices])^2, R=1000)
# print("REVEL versus observed literature BRS1 diagnosis probability")
# foo$t0
# quantile(foo$t,c(0.025,0.975))
```

## Forest Plot with only lqt2_dist data 
```{r preamble,}
# # code for another forest plot, ask K if he wants a different one/ use for just brs1/ just lqt3 when relevant 
#load("all.lqt2.data.trunc.Rdata")
# 
# penetrance_lqt2 <- na.omit(herg.combined.data$penetrance_lqt2)
# 
# yList=c("penetrance_lqt2")
# 
# xList=c('lqt2_dist_sequence','lqt2_dist',  'lqt2_dist_max', 'lqt2_dist_mean' ,'lqt2_dist_obs_seq', 'lqt2_dist_observed','lqt2_dist_obs_max', 'lqt2_dist_obs_mean', 'revel_score')
# 
# weight <- c(herg.combined.data$weight)
# weight <- weight[weight!=0]
# 
# # resultTable<-calcAllPvals(yList, xList, 1000, 'weight', results.1.data[results.1.data$total_carriers>0,])
# # rm(results.1.data)
# 
# #4
# 
# i=0
# FP.data<-data.frame()
# for (x in xList){
#   i=i+2
#   FP.data[i,"Feature"]<-x
#   t<-herg.combined.data[!is.na(herg.combined.data[,x]) & herg.combined.data$total_carriers>0,]
#   
#   #<-t[!is.na(t[,"var"]),]
#   #FP.data[i,"Feature"]<-paste(x,"_cohort")
#   #t<-herg.combined.data[!is.na(herg.combined.data[,x]) & herg.combined.data$total_carriers>0,]
#   
#   t<-t[!is.na(t[,"var"]),]
#   foo <- boot(t, function(data,indices)
#     weightedCorr(t[indices,x],t$penetrance_lqt2[indices], method="Spearman", weights = t$weight[indices]), R=1000)
#   
#   
#   FP.data[i,"Spearman"]<-foo$t0
#   FP.data[i,"Spearman_low"]<-quantile(foo$t,c(0.025,0.975), na.rm = T)[1][[1]]
#   FP.data[i,"Spearman_high"]<-quantile(foo$t,c(0.025,0.975), na.rm = T)[2][[1]]
#   FP.data[i,"n"]<-length(t[,x])
# }
# 
# #5
# 
# forestplot(FP.data$Feature,FP.data$Spearman,FP.data$Spearman_low,FP.data$Spearman_high, cex = 0.5)
# 
# #rm(FP.data)

```
# Part 4: ROC plots, AUC plots and PR curves.

## Observed cohort LQT2 diagnosis probability with multiple observations variants.
### ROC Plot

```{r, echo=FALSE, warning=FALSE, message=FALSE}
load("all.brs1.data.trunc.Rdata")

scn5a.cln.d$brs1_patho <- NA
scn5a.cln.d$brs1_patho[scn5a.cln.d$brs1_penetrance>=0.5] <- 1  #0.5 penetrance is threshold for determining pathogenesis 
scn5a.cln.d$brs1_patho[scn5a.cln.d$brs1_penetrance<0.5] <- 0

fglm<-scn5a.cln.d[!is.na(scn5a.cln.d$blastpssm),] #& herg.combined.data$total_carriers == 1,] #store data omitting blast_pssm NA values

colrs <- c("red", "green", "orange", "blue", "black", "grey")
colrs2 <- c("magenta", "red", "green", "orange", "blue", "black", "gray")
#produce 6 generalized linear models based on different in silico predictor models
modfunc1 <- glm(brs1_patho~pph2_prob, data = fglm, family = 'binomial') 
# rcs(revel_score,3)*
modfunc2 <- glm(brs1_patho~REVEL, data = fglm, family = 'binomial')
modfunc3 <- glm(brs1_patho~blastpssm, data = fglm, family = 'binomial')
modfunc4 <- glm(brs1_patho~SIFT_score, data = fglm, family = 'binomial')
modfunc5 <- glm(brs1_patho~brs1_dist, data = fglm, family = 'binomial')
#modfunc6 <- glm(brs1_patho~p_mean_w, data = fglm, family = 'binomial')
modfunc6 <- glm(brs1_patho~brs1_dist_sequence+brs1_dist_correlated, data = fglm, family = 'binomial')

funs <- list(modfunc2, modfunc3, modfunc4, modfunc5, modfunc6) 

i=0
par(pty="s")

ROC.tmp.data<-roc(fglm$brs1_patho[row.names(fglm) %in% as.numeric(names(predict(modfunc1)))], predict(modfunc1), ci=T, quiet=T)

plot.roc(ROC.tmp.data, col = "magenta", legacy.axes = T, xlab="False Positive Rate", ylab = "True Positive Rate", main = "ROC of Multiple Observations of Variants")

auc<- ROC.tmp.data$auc

for(li in funs){ #iterate through list 
  i=i+1
  ROC.tmp.data<-roc(fglm$brs1_patho[row.names(fglm) %in% as.numeric(names(predict(li)))], predict(li), ci=T, quiet=T)
  plot.roc(ROC.tmp.data,add = T, col = colrs[i], legacy.axes = T,xlab="False Positive Rate", ylab = "True Positive Rate")
  auc <- c(auc, ROC.tmp.data$auc)
  
}
names<- c("pph2_prob:", "REVEL:", "blastpssm:", "SIFT_score:", "brs1_dist:", "combined_dist:")

auc <- round(auc, digits = 2)
leg.data <-data.frame(names,auc)
leg.data$together <- paste(leg.data$names, leg.data$auc)

legend("bottomright", title = "Legend (AUC)", legend = leg.data$together, col = colrs, lty= 1, text.font = 2, cex = 0.7, box.lwd = 0,)

```

###Precision-Recall Curve

```{r, echo=FALSE, warning=FALSE, message=FALSE}
# herg.combined.data$lqt2_patho <- NA
# herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2>=0.5] <- 1  #0.5 penetrance is threshold for determining pathogenesis 
# herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2<0.5] <- 0
# class0 = herg.combined.data$lqt2_patho = 1
# class1 = herg.combined.data$lqt2_patho = 0
# pr<-pr.curve(class0, class1, curve=T)
# plot(pr)

#legend and auc printing 

df <- data.frame(fglm$pph2_prob, fglm$REVEL, fglm$blastpssm, fglm$SIFT_score, fglm$brs1_dist,fglm$brs1_dist_sequence+fglm$brs1_dist_correlated)

pr <- pr.curve(na.omit(df[,1]), weights.class0 = fglm$brs1_patho, curve = TRUE)
plot(pr, col = "magenta", lty = 1, main = "PR CURVE for Multiple Observations", auc.main = FALSE)
auc.list <- pr$auc.integral
n =length(df)
i = 0
for(var in 2:n){
  i = i+1
  par(new=TRUE)
  pr <- pr.curve(na.omit(df[,var]), weights.class0 = fglm$brs1_patho, curve = TRUE)
  plot(pr, col = colrs[i], lty = 1, main = "", auc.main = FALSE)
  auc.list<- c(auc.list, pr$auc.integral)
  
}

par(new = FALSE)

auc.list <- round(auc.list, digits = 2)
leg.data <-data.frame(names,auc.list)
leg.data$together <- paste(leg.data$names, leg.data$auc.list)

legend("bottomright", title = "Legend (AUC)", legend = leg.data$together, col = colrs2, lty= 1, text.font = 2, cex = 0.7, box.lwd = 0,)
```

#### AUC's from ROC's observed literature LQT2 diagnosis probability at multiple cutoffs

```{r, echo=FALSE, warning=FALSE, message=FALSE}

# Wrapper for convenience of making GLMs and outputting AUCs
glm.mod=function(fglm,independent){
  in_string <- paste(independent, collapse=" + ")
  regression_formula <- as.formula(paste("brs1_patho", in_string, sep=" ~ "))
  mod <- glm(regression_formula, data = fglm, family = 'binomial')
  #print(paste(regression_formula, tmp$auc))
  ROC.tmp.data<-roc(fglm$brs1_patho[row.names(fglm) %in% as.numeric(names(predict(mod)))], predict(mod), ci=T)
  return(ROC.tmp.data$auc)
}

brs1_patho <- NA
#scn5a.cln.d$p_mean_prior<- NA

scn5a.cln.d$brs1_patho[scn5a.cln.d$brs1_penetrance>=0.5] <- 1
scn5a.cln.d$brs1_patho[scn5a.cln.d$brs1_penetrance<0.5] <- 0
#scn5a.cln.d$p_mean_prior<-scn5a.cln.d$alpha/(scn5a.cln.d$alpha+scn5a.cln.d$beta)
cutoffs <- seq(0.1,0.8,0.05)
#colrs <- c("red", "green", "orange", "blue", "black", "gray") redundant
fglm<-scn5a.cln.d[!is.na(scn5a.cln.d$brs1_patho) & !is.na(scn5a.cln.d$blastpssm), ]# & combined.data$total_carriers == 1
dist.only <- 0
blast1 <- 0
polyphen1 <- 0
SIFT_score <- 0
EM <- 0
revel<-0
len.new<-0
i=0
for (co in cutoffs){
  fglm$brs1_patho[fglm$brs1_penetrance>=co] <- 1
  fglm$brs1_patho[fglm$brs1_penetrance<co] <- 0
  print(paste(length(fglm$brs1_patho), " ", sum(fglm$brs1_patho)))
  if (!sum(fglm$brs1_patho)<6){
    i=i+1
    blast1[i]<-glm.mod(fglm,"blastpssm")
    polyphen1[i]<-glm.mod(fglm,"pph2_prob")
    revel[i]<-glm.mod(fglm,"REVEL")
    SIFT_score[i]<-glm.mod(fglm,"SIFT_score")
    dist.only[i]<-glm.mod(fglm,"brs1_dist")
    #EM[i]<-glm.mod(fglm,"p_mean_prior")
    #test[i]<- glm.mod(fglm, "lqt2_dist_sequence"*"lqt2_dist_observed")
  }
}
par(cex=1, bty='l', lwd=2)

plot(cutoffs[1:i],polyphen1,col="magenta",type = "l",ylim = c(0.5,1), ylab = "AUC", xlab = "Cutoffs", main = "AUC of Multiple Observations of Variants") 

lines(cutoffs[1:i],revel,col=colrs[1]) # red
lines(cutoffs[1:i],blast1,col=colrs[2]) # green
lines(cutoffs[1:i],SIFT_score,col=colrs[3]) # orange
lines(cutoffs[1:i],dist.only,col=colrs[4]) # blue
#lines(cutoffs[1:i],EM,col=colrs[5]) # black
```

## Observed cohort BRS1 diagnosis probability with single observation variants.
## ROC Plot

```{r, echo=FALSE, warning=FALSE, message=FALSE}


#Second ROC plot 

n1<- scn5a.cln.d[scn5a.cln.d$total_carriers==1,]
scn5a.cln.d<- scn5a.cln.d[scn5a.cln.d$total_carriers>1,]
f<-data.frame(resnum=NA,brs1_dist=NA, brs1_dist_weight=NA)
i<-0
for(rec in seq(1,1149,1)){
  i <- i+1
  f[rec, c('brs1_dist', 'brs1_dist_weight')] <- funcdist(rec, "blaa", scn5a.cln.d, scn5a.cln.d$h2dist, "brs1_penetrance", "sigmoid", 7)
  f[i, "resnum"] <- rec
  # f[i, "lqt2_dist"] <- dist
  # f[i, "lqt2_dist_weight"] <- dist_weight
}

n1 <- merge(x = n1, y = f, by = "resnum", all = T)
n1$brs1_dist.x[!is.na(n1$brs1_dist.y)] <-
  n1$brs1_dist.y[!is.na(n1$brs1_dist.y)]
n1$brs1_dist_weight.x[!is.na(n1$brs1_dist_weight.y)] <-
  n1$brs1_dist_weight.y[!is.na(n1$brs1_dist_weight.y)]

#changed from herg.combined.data to n1 for now

n1$brs1_patho <- NA
n1$brs1_patho[n1$brs1_penetrance>=0.5] <- 1
n1$brs1_patho[n1$brs1_penetrance<0.5] <- 0
fglm<-n1[!is.na(n1$blastpssm) & n1$total_carriers == 1, ]


fglm <- fglm[,!names(fglm) %in%  "cardiacboost"]
fglm <- unique(fglm)
#colrs <- c("red", "green", "orange", "blue", "black", "gray") redundant
modfunc1 <- glm(brs1_patho~pph2_prob, data = fglm, family = 'binomial')
modfunc2 <- glm(brs1_patho~REVEL, data = fglm, family = 'binomial')
modfunc3 <- glm(brs1_patho~blastpssm, data = fglm, family = 'binomial')
modfunc4 <- glm(brs1_patho~SIFT_score, data = fglm, family = 'binomial')
modfunc5 <- glm(brs1_patho~brs1_dist.x, data = fglm, family = 'binomial')
#modfunc6 <- glm(brs1_patho~p_mean_w, data = fglm, family = 'binomial')
funs <- list(modfunc2, modfunc3, modfunc4, modfunc5)

i=0
par(pty="s")
ROC.tmp.data<-roc(fglm$brs1_patho[row.names(fglm) %in% as.numeric(names(predict(modfunc1)))], predict(modfunc1), ci=T, quiet=T)

plot.roc(ROC.tmp.data, col = "magenta",legacy.axes = T, xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC of Single Observation of Variants")

auc<- ROC.tmp.data$auc

for(li in funs){
  i=i+1
  ROC.tmp.data<-roc(fglm$brs1_patho[row.names(fglm) %in% as.numeric(names(predict(li)))], predict(li), ci=T, quiet=T)
  plot.roc(ROC.tmp.data,add = T, col = colrs[i],legacy.axes = T,xlab="False Positive Rate", ylab = "True Positive Rate")
  print(ROC.tmp.data$auc)
  auc <- c(auc, ROC.tmp.data$auc)
}

names<- c("pph2_prob:", "REVEL:", "blastpssm:","SIFT_score:", "brs1_dist:")
#added for legend
auc <- round(auc, digits = 2)
leg.data <-data.frame(names,auc)
leg.data$together <- paste(leg.data$names, leg.data$auc)
legend("bottomright", title = "Legend (AUC)", legend = leg.data$together, col = colrs, lty= 1, text.font = 2, cex = 0.7, box.lwd = 0,)

```

###PR Curve

```{r, echo=FALSE, warning=FALSE, message=FALSE, include=TRUE}
# herg.combined.data$lqt2_patho <- NA
# herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2>=0.5] <- 1  #0.5 penetrance is threshold for determining pathogenesis 
# herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2<0.5] <- 0
# class0 = herg.combined.data$lqt2_patho = 1
# class1 = herg.combined.data$lqt2_patho = 0
# pr<-pr.curve(class0, class1, curve=T)
# plot(pr)

#legend and auc printing 

df <- data.frame(fglm$pph2_prob, fglm$revel_score, fglm$blast_pssm, fglm$RMSF, fglm$lqt2_dist.x, fglm$p_mean_w)

pr <- pr.curve(na.omit(df[,1]), weights.class0 = fglm$lqt2_patho, curve = TRUE)
plot(pr, col = "magenta", lty = 1, main = "PR CURVE for Single Observations", auc.main = FALSE)
auc.list <- pr$auc.integral
n =length(df)
i = 0
for(var in 2:n){
  i = i+1
  par(new=TRUE)
  pr <- pr.curve(na.omit(df[,var]), weights.class0 = fglm$lqt2_patho, curve = TRUE)
  plot(pr, col = colrs[i], lty = 1, main = "", auc.main = FALSE)
  auc.list<- c(auc.list, pr$auc.integral)
  
}

par(new = FALSE)

auc.list <- round(auc.list, digits = 2)
leg.data <-data.frame(names,auc.list)
leg.data$together <- paste(leg.data$names, leg.data$auc.list)

legend("bottomright", title = "Legend (AUC)", legend = leg.data$together, col = colrs2, lty= 1, text.font = 2, cex = 0.7, box.lwd = 0,)
```

#### AUC for N = 1 variants from the literature.

```{r, echo=FALSE, warning=FALSE, message=FALSE}
#load("data/Covariates/KCNH2_clinical_data.RData")

# Wrapper for convenience of making GLMs and outputting AUCs
glm.mod=function(fglm,independent){
  in_string <- paste(independent, collapse=" + ")
  regression_formula <- as.formula(paste("lqt2_patho", in_string, sep=" ~ "))
  mod <- glm(regression_formula, data = fglm, family = 'binomial')
  #print(paste(regression_formula, tmp$auc))
  ROC.tmp.data<-roc(fglm$lqt2_patho[row.names(fglm) %in% as.numeric(names(predict(mod)))], predict(mod), ci=T)
  return(ROC.tmp.data$auc)
}

#changed from herg.combined data to n1 for now
n1$lqt2_patho <- NA
n1$lqt2_patho[n1$penetrance_lqt2>=0.5] <- 1
n1$lqt2_patho[n1$penetrance_lqt2<0.5] <- 0
n1$p_mean_prior<-n1$alpha/(n1$alpha+n1$beta)
cutoffs <- seq(0.1,0.8,0.05)
colrs <- c("red", "green", "orange", "blue", "black", "gray")
fglm<-n1[!is.na(n1$lqt2_patho) & !is.na(n1$blast_pssm) & n1$total_carriers == 1,]
dist.only <- 0
blast1 <- 0
polyphen1 <- 0
RMSF <- 0
EM <- 0
revel<-0

len.new<-0
i=0
for (co in cutoffs){
  fglm$lqt2_patho[fglm$penetrance_lqt2>=co] <- 1
  fglm$lqt2_patho[fglm$penetrance_lqt2<co] <- 0
  print(paste(length(fglm$lqt2_patho), " ", sum(fglm$lqt2_patho)))
  if (!sum(fglm$lqt2_patho)<6){
    i=i+1
    blast1[i]<-glm.mod(fglm,"blast_pssm")
    polyphen1[i]<-glm.mod(fglm,"pph2_prob")
    revel[i]<-glm.mod(fglm,"revel_score")
    RMSF[i]<-glm.mod(fglm,"RMSF")
    dist.only[i]<-glm.mod(fglm,"lqt2_dist.x")
    EM[i]<-glm.mod(fglm,"p_mean_prior")
  }
}
par(cex=1, bty='l', lwd=2)

plot(cutoffs[1:i],polyphen1,col="magenta",type = "l",ylim = c(0.5,1), ylab = "AUC", xlab = "Cutoffs", main = "AUC of Single Observation of Variants")

lines(cutoffs[1:i],revel,col=colrs[1]) # red
lines(cutoffs[1:i],blast1,col=colrs[2]) # green
lines(cutoffs[1:i],RMSF,col=colrs[3]) # orange
lines(cutoffs[1:i],dist.only,col=colrs[4]) # blue
lines(cutoffs[1:i],EM,col=colrs[5]) # black

```


## All Lqt2_dist Predictions
###ROC Plot
```{r, echo=FALSE, warning=FALSE, message=FALSE}
load("all.lqt2.data.trunc.Rdata")

herg.combined.data$lqt2_patho <- NA
herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2>=0.5] <- 1  #0.5 penetrance is threshold for determining pathogenesis 
herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2<0.5] <- 0

fglm<-herg.combined.data[!is.na(herg.combined.data$blast_pssm),] #& herg.combined.data$total_carriers == 1,] #store data omitting blast_pssm NA values

colrs <- c("red", "green", "orange", "blue", "black", "purple", "brown", "yellow2","skyblue", 'cyan' )
colrs2 <- c("magenta", "red", "green", "orange", "blue", "black", "purple", "brown", "yellow2","skyblue", 'cyan')
#produce 6 generalized linear models based on different in silico predictor models
modfunc1 <- glm(lqt2_patho~lqt2_dist, data = fglm, family = 'binomial') 
modfunc2 <- glm(lqt2_patho~lqt2_dist_sequence, data = fglm, family = 'binomial')
modfunc3 <- glm(lqt2_patho~lqt2_dist_observed, data = fglm, family = 'binomial')
modfunc4 <- glm(lqt2_patho~lqt2_dist_max, data = fglm, family = 'binomial')
modfunc6 <- glm(lqt2_patho~lqt2_dist_mean, data = fglm, family = 'binomial')
modfunc7 <- glm(lqt2_patho~lqt2_dist_obs_max, data = fglm, family = 'binomial')
modfunc9 <- glm(lqt2_patho~lqt2_dist_obs_mean, data = fglm, family = 'binomial')
modfunc10 <- glm(lqt2_patho~lqt2_dist_obs_seq, data = fglm, family = 'binomial')
modfunc11 <- glm(lqt2_patho~revel_score, data = fglm, family = 'binomial')



funs <- list(modfunc2, modfunc3, modfunc4, modfunc6, modfunc7, modfunc9, modfunc10, modfunc11 ) 

i=0
par(pty="s")

ROC.tmp.data<-roc(fglm$lqt2_patho[row.names(fglm) %in% as.numeric(names(predict(modfunc1)))], predict(modfunc1), ci=T, quiet=T)

plot.roc(ROC.tmp.data, col = "magenta", legacy.axes = T, xlab="False Positive Rate", ylab = "True Positive Rate", main = "ROC for lqt2_dist values for multiple observations of variants")

auc<- ROC.tmp.data$auc

for(li in funs){ #iterate through list 
  i=i+1
  ROC.tmp.data<-roc(fglm$lqt2_patho[row.names(fglm) %in% as.numeric(names(predict(li)))], predict(li), ci=T, quiet=T)
  plot.roc(ROC.tmp.data,add = T, col = colrs[i], legacy.axes = T,xlab="False Positive Rate", ylab = "True Positive Rate")#,print.auc = T, print.auc.y=50-(i*10)) 
  auc <- c(auc, ROC.tmp.data$auc)
  
}
names<- c("lqt2_dist:", "lqt2_dist_sequence:", "lqt2_dist_observed:", "lqt2_dist_max:",  "lqt2_dist_mean:", "lqt2_dist_obs_max:",  "lqt2_dist_obs_mean:","lqt2_dist_obs_seq:", "revel_score:")

auc <- round(auc, digits = 2)
leg.data <-data.frame(names,auc)
leg.data$together <- paste(leg.data$names, leg.data$auc)

legend("bottomright", title = "Legend (AUC)", legend = leg.data$together, col = colrs2, lty= 1, text.font = 2, cex = 0.7, box.lwd = 0,)


```
###Precision-Recall Curve

```{r, echo=FALSE, warning=FALSE, message=FALSE}
# herg.combined.data$lqt2_patho <- NA
# herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2>=0.5] <- 1  #0.5 penetrance is threshold for determining pathogenesis 
# herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2<0.5] <- 0
# class0 = herg.combined.data$lqt2_patho = 1
# class1 = herg.combined.data$lqt2_patho = 0
# pr<-pr.curve(class0, class1, curve=T)
# plot(pr)

#legend and auc printing 

df <- data.frame(fglm$lqt2_dist, fglm$lqt2_dist_sequence, fglm$lqt2_dist_observed, fglm$lqt2_dist_max,  fglm$lqt2_dist_mean, fglm$lqt2_dist_obs_max, fglm$lqt2_dist_obs_mean, fglm$lqt2_dist_obs_seq, fglm$revel_score)

pr <- pr.curve(na.omit(df[,1]), weights.class0 = fglm$lqt2_patho, curve = TRUE)
plot(pr, col = "magenta", lty = 1, main = "PR CURVE for lqt2_dist values for multiple observations of variants", auc.main = FALSE)
auc.list <- pr$auc.integral
n =length(df)
i = 0
for(var in 2:n){
  i = i+1
  par(new=TRUE)
  pr <- pr.curve(na.omit(df[,var]), weights.class0 = fglm$lqt2_patho, curve = TRUE)
  plot(pr, col = colrs[i], lty = 1, main = "", auc.main = FALSE)
  auc.list<- c(auc.list, pr$auc.integral)
  
}

par(new = FALSE)

auc.list <- round(auc.list, digits = 2)
leg.data <-data.frame(names,auc.list)
leg.data$together <- paste(leg.data$names, leg.data$auc.list)

legend("bottomright", title = "Legend (AUC)", legend = leg.data$together, col = colrs2, lty= 1, text.font = 2, cex = 0.7, box.lwd = 0,)
```

####AUC For All lqt2_dist predictions

```{r, echo=FALSE, warning=FALSE, message=FALSE}

# Wrapper for convenience of making GLMs and outputting AUCs
glm.mod=function(fglm,independent){
  in_string <- paste(independent, collapse=" + ")
  regression_formula <- as.formula(paste("lqt2_patho", in_string, sep=" ~ "))
  mod <- glm(regression_formula, data = fglm, family = 'binomial')
  #print(paste(regression_formula, tmp$auc))
  ROC.tmp.data<-roc(fglm$lqt2_patho[row.names(fglm) %in% as.numeric(names(predict(mod)))], predict(mod), ci=T)
  return(ROC.tmp.data$auc)
}

lqt2_patho <- NA

herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2>=0.5] <- 1
herg.combined.data$lqt2_patho[herg.combined.data$penetrance_lqt2<0.5] <- 0
herg.combined.data$p_mean_prior<-herg.combined.data$alpha/(herg.combined.data$alpha+herg.combined.data$beta)
cutoffs <- seq(0.1,0.8,0.05)
colrs <- c("red", "green", "orange", "blue", "black", "purple", "brown", "yellow2","skyblue", "cyan")
fglm<-herg.combined.data[!is.na(herg.combined.data$lqt2_patho) & !is.na(herg.combined.data$blast_pssm), ]# & combined.data$total_carriers == 1

lqt2_dist <- 0 
lqt2_dist_sequence <- 0 
lqt2_dist_observed <- 0 
lqt2_dist_max <- 0 

lqt2_dist_mean <- 0 
lqt2_dist_obs_max <- 0 

lqt2_dist_obs_mean <- 0 
lqt2_dist_obs_seq<-0
revel_score <- 0 

i=0
for (co in cutoffs){
  fglm$lqt2_patho[fglm$penetrance_lqt2>=co] <- 1
  fglm$lqt2_patho[fglm$penetrance_lqt2<co] <- 0
  print(paste(length(fglm$lqt2_patho), " ", sum(fglm$lqt2_patho)))
  if (!sum(fglm$lqt2_patho)<6){
    i=i+1
    lqt2_dist[i]<-glm.mod(fglm,"lqt2_dist")
    lqt2_dist_sequence[i]<-glm.mod(fglm,"lqt2_dist_sequence")
    lqt2_dist_observed[i]<-glm.mod(fglm,"lqt2_dist_observed")
    lqt2_dist_max[i]<-glm.mod(fglm,"lqt2_dist_max")
    
    lqt2_dist_mean[i]<-glm.mod(fglm,"lqt2_dist_mean")
    lqt2_dist_obs_max[i]<-glm.mod(fglm,"lqt2_dist_obs_max")
    
    lqt2_dist_obs_mean[i]<-glm.mod(fglm,"lqt2_dist_obs_mean")
    lqt2_dist_obs_seq[i]<- glm.mod(fglm,"lqt2_dist_obs_seq")
    revel_score[i]<-glm.mod(fglm,"revel_score")
  }
}
par(cex=1, bty='l', lwd=2)

plot(cutoffs[1:i],lqt2_dist,col="magenta",type = "l",ylim = c(0.5,1), ylab = "AUC", xlab = "Cutoffs", main = "AUC for lqt2_dist values for multiple observations of variants") 

lines(cutoffs[1:i],lqt2_dist_sequence,col=colrs[1]) # red
lines(cutoffs[1:i],lqt2_dist_observed,col=colrs[2]) # green
lines(cutoffs[1:i],lqt2_dist_max,col=colrs[3]) # orange
lines(cutoffs[1:i],lqt2_dist_mean,col=colrs[4]) # black
lines(cutoffs[1:i],lqt2_dist_obs_max,col=colrs[5]) # purple
lines(cutoffs[1:i],lqt2_dist_obs_mean,col=colrs[6]) # yellow
lines(cutoffs[1:i],lqt2_dist_obs_seq,col=colrs[7]) # light blue
lines(cutoffs[1:i],revel_score,col=colrs[8]) # cyan

```