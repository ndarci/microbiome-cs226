# compare performance of different regression methods for LDLHDL
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(glmnet)

# for principal components regression:
#install.packages("pls")
library(pls)

setwd("~/src/microbiome-cs226/data/")

ab = read.table("baseline_genusAbundance.txt")
ab_pheno = read.table("baseline_genusAbundance_pheno.txt", header = T)
pheno = ab_pheno[1:9]

ab_pheno_all = read.table("allTimePoints_abundance_pheno_visitInfo.txt", header = T, sep = ',')
ab_pheno_all = cbind(select(ab_pheno_all, starts_with("genus")),  ab_pheno_all[,c("SampleID", "SubjectID", "CollectionDate", "LDLHDL")])

ab_pheno_all = na.omit(ab_pheno_all)
pheno_all = ab_pheno_all[c("SampleID", "SubjectID", "CollectionDate", "LDLHDL")]
ab_all = select(ab_pheno_all, starts_with("genus"))
rownames(ab_all) = ab_pheno_all$SampleID

##########################
# PENALIZED REGRESSION, ALL TIME POINTS
##########################

# goal: fit a lasso regression on full data

# first we evaluate ordinary linear regression as a benchmark
data <- cbind(ab_pheno_all[,c("LDLHDL","CollectionDate")], ab_pheno_all %>% select(starts_with("genus")))
data <- as.matrix(data)

# 60-40 training split - same for each method we will test!!!!
set.seed(222)
train_rows = sample(1:nrow(data), round(nrow(data)*.6))

x.train = data[train_rows, 2:ncol(data)]
y.train = data[train_rows, 1]

x.test = data[-train_rows, 2:ncol(data)]
y.test = data[-train_rows, 1]

# fit full and null models
full.lm <- lm(y.train ~ . , data = as.data.frame(x.train))
# summary(fit.lm)

null.lm <- lm(y.train ~ 1)
# summary(null.lm)

# run thru test set
yhat.null <- predict(null.lm, newdata = as.data.frame(x.test))
yhat.full <- predict(full.lm, newdata = as.data.frame(x.test))
mse.null <- mean((y.test - yhat.null)^2)
mse.full <- mean((y.test - yhat.full)^2)

# mse of null model
mse.null
# mse of saturated model
mse.full

# fit penalized regressions for alpha = 0, 0.1, ..., 0.9, 1
for (i in 0:10) {
  assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse", 
                                            alpha=i/10,family="gaussian"))
}

# The following code can be used to make pretty plots

# fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1)
# 
# fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0)
# 
# fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5)
# 
# par(mfrow=c(3,2))
# 
# plot(fit.lasso, xvar="lambda")
# plot(fit10, main="LASSO")
# 
# plot(fit.ridge, xvar="lambda")
# plot(fit0, main="Ridge")
# 
# plot(fit.elnet, xvar="lambda")
# plot(fit5, main="Elastic Net")


# run our 10 models thru the testing data
yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test) # ridge
yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=x.test)
yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=x.test)
yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=x.test)
yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=x.test)
yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)
yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=x.test)
yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=x.test)
yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=x.test)
yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=x.test)
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test) # lasso

# using mspe as accuracy metric
mse0 <- mean((y.test - yhat0)^2)
mse1 <- mean((y.test - yhat1)^2)
mse2 <- mean((y.test - yhat2)^2)
mse3 <- mean((y.test - yhat3)^2)
mse4 <- mean((y.test - yhat4)^2)
mse5 <- mean((y.test - yhat5)^2)
mse6 <- mean((y.test - yhat6)^2)
mse7 <- mean((y.test - yhat7)^2)
mse8 <- mean((y.test - yhat8)^2)
mse9 <- mean((y.test - yhat9)^2)
mse10 <- mean((y.test - yhat10)^2)

best_lambda = c(fit0$lambda.1se, fit1$lambda.1se, fit2$lambda.1se,
                fit3$lambda.1se, fit4$lambda.1se, fit5$lambda.1se,
                fit6$lambda.1se, fit7$lambda.1se, fit8$lambda.1se,
                fit9$lambda.1se, fit10$lambda.1se)

alpha = seq(0,1,.1)

MSPE = c(mse0,mse1,mse2,mse3,mse4,mse5,mse6,mse7,mse8,mse9,mse10)

penalized_results = cbind(alpha, best_lambda, MSPE )
penalized_results

#elastic net is best by mse on test set
penalized_results[6,]
coef(fit5)

mse5
# compare to
mse.full
mse.null
# a marked improvement.....

# could visualize reduction in coefs as a fxn of alpha

##########################
# PRINCIPAL COMPONENTS REGRESSION, ALL TIME POINTS
##########################

# note: using same train/test split as in the above!!!
m1 <- pcr(y.train ~ . - CollectionDate, data = as.data.frame(x.train), validation = "CV")
summary(m1)

# interesting shapes... can 30-something components really be most optimal?
# validationplot(m1)
# validationplot(m1, val.type="MSEP")
# validationplot(m1, val.type="R2")
# still only captures ~20% of the variation in LDLHDL....


# run thru test data
pcr_pred <- predict(m1, as.data.frame(x.test))
mse.pcr = mean((pcr_pred - y.test)^2)


##########################
# MEDIAN BRAY CURTIS AS PREDICTOR, ALL TIMES
##########################

# gotta run this first:

braycurtis <- function(si, sj) {
  # sum lesser proportions of each genus
  top = sum(abs(si - sj))
  bottom = sum(si + sj)
  return(top/bottom)
}

# function to get dysbiosis scores given an abundance table and response variable
# target must be oriented so lower values are healthier
dysbiosisAnalysis <- function(abund, pheno, target) {
  # scale abundances to proportions
  abund_scale = t(apply(abund, 1, function(row) row / sum(row)))
  # rowSums(abund_scale) # verify this worked
  
  # define "healthy" samples, from people with low cholesterol
  healthythresh = 0.1
  healthy = rownames(abund)[target < quantile(target, healthythresh)]
  unhealthy = setdiff(rownames(abund), healthy)
  
  # calculate median bray-curtis dissimilarity from healthy for each sample
  medianbc = sapply(rownames(abund), function(samp) 
    median(sapply(healthy, function(h) braycurtis(abund[samp,], abund[h,])))
  )
  
  # find dysbiotic samples
  dysthresh = 0.9
  dys = pheno[medianbc > quantile(medianbc, dysthresh), "SampleID"]
  
  # merge this info back onto the phenotype data
  pheno = cbind(pheno, medianbc)
  pheno$healthy = pheno$SampleID %in% healthy
  pheno$dysbiotic = pheno$SampleID %in% dys
  
  return(list(medianbc = medianbc, pheno = pheno))
}

dysbio_all = dysbiosisAnalysis(ab_all, pheno_all, ab_pheno_all$LDLHDL)
medianbc = dysbio_all$medianbc
pheno = dysbio_all$pheno

# still using same train/test split as in the above!!!

lm1 <- (lm(LDLHDL ~ medianbc , data = pheno[train_rows,]))
bc_pred <- predict(lm1, newdata = (pheno[-train_rows,]))
mse.bc = mean((bc_pred - y.test)^2)

# and the winner is... elastic net
c(mse.null,mse.full,mse5,mse.pcr,mse.bc)

# plot results in a barplot
msevec = c(mse.null, mse.full, mse0, mse10, mse5, mse.pcr, mse.bc)
modelnames = c("Null", "Full", "Ridge", "LASSO", "ElasticNet (alpha=0.5)", 
               "PCR Regression", "Dysbiosis score regression")
barplotdf = data.frame("Model" = modelnames, "MSE" = msevec)
barplotdf = barplotdf[order(barplotdf$MSE, decreasing = T),]
barplot_comp = ggplot(barplotdf, aes(y = Model, x = MSE)) + 
  geom_bar(stat = 'identity') + 
  scale_y_discrete(limits = barplotdf$Model)
ggsave("../fig/model_comparison_barplot.png", barplot_comp, dpi = 200)
