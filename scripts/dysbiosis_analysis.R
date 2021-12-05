# test if gut microbe dysbiosis is connected to high cholesterol
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(Hmisc)
library(corrplot)

setwd("~/src/microbiome-cs226/scripts/")

# import data
# ab = read.table("../data/baseline_genusAbundance.txt")
# ab_pheno = read.table("../data/baseline_genusAbundance_pheno.txt", header = T)
# pheno = ab_pheno[1:9]

ab_pheno_all = read.table("../data/allTimePoints_abundance_pheno_visitInfo.txt", header = T, sep = ',')
desiredpheno = c("SampleID", "SubjectID", "CollectionDate", 
                 "LDLHDL", "Age", "Sex", "BMI", "Race")
# make sure we have complete data
ab_pheno_all = cbind(select(ab_pheno_all, starts_with("genus")), 
                     ab_pheno_all[desiredpheno])
ab_pheno_all = na.omit(ab_pheno_all)
# separate abundance from phenotype
ab_all = select(ab_pheno_all, starts_with("genus"))
pheno_all = ab_pheno_all[desiredpheno]
rownames(ab_all) = ab_pheno_all$SampleID

# summary(lm(ab_pheno_all$LDLHDL ~ ., data = ab_all))
# 
# pcs_all = prcomp(ab_all)
# ggplot(data.frame(pcs_all$x)) + geom_point(aes(x = PC1, y = PC2, 
#                                                color = ab_pheno_all$LDLHDL > mean(ab_pheno_all$LDLHDL)))

# function to calculate bray-curtis dissimilarity between two samples (assumes they are counted as proportions)
braycurtis <- function(si, sj) {
  # sum lesser proportions of each genus
  top = sum(abs(si - sj))
  bottom = sum(si + sj)
  return(top/bottom)
}

# function to coerce the right rows into bray-curtis and take median result
braycurtis_wrap <- function(abund, sample, healthy, pheno) {
  samesubject = pheno[pheno$SampleID == sample, "SubjectID"]
  res = sapply(healthy[!(healthy %in% samesubject)], function(h) 
    braycurtis(abund[h,], abund[sample,]))
  return(median(res))
}

# function to get dysbiosis scores given an abundance table and response variable
# target must be oriented so lower values are healthier
dysbiosisAnalysis <- function(abund, pheno, target) {
  # define "healthy" samples, from people with low cholesterol
  healthythresh = 0.1
  healthy = rownames(abund)[target < quantile(target, healthythresh)]
  unhealthy = setdiff(rownames(abund), healthy)
  
  # calculate median bray-curtis dissimilarity from healthy for each sample
  # exclude samples from the same subject
  # medianbc = sapply(rownames(abund), function(samp) 
  #   median(sapply(healthy, function(h) braycurtis(abund[samp,], abund[h,])))
  # )
  
  
  medianbc = sapply(rownames(abund), function(sample)
      braycurtis_wrap(abund, sample, healthy, pheno)
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

dysbio_all = dysbiosisAnalysis(ab_all, pheno_all, pheno_all$LDLHDL)
medianbc = dysbio_all$medianbc
pheno = dysbio_all$pheno

# compute PCs to correct for in regression
pca_all = prcomp(ab_all, center = T, scale. = T)

# find correlations with PCs and other important factors
corfeatures = cbind(pca_all$x[,1:10], 
                    pheno[c("LDLHDL", "medianbc", "Age", "Sex", "BMI")])
corfeatures$Sex = as.numeric(as.factor(corfeatures$Sex))
cor = rcorr(as.matrix(corfeatures), type = "pearson")
# fix missing diagonal
for (i in 1:nrow(cor$P)) {
    cor$P[i, i] <- 0
}
png("../fig/pca_correlation_heatmap.png", width = 1500, height = 1500, res = 200)
print(corrplot(cor$r, type = "upper", p.mat = cor$P, sig.level = 0.05/nrow(cor$r),
         pch.cex = 1, pch.col = "gray",
         tl.cex = 0.75, tl.col = "black",
         method = "color"))
dev.off()

ggplot(corfeatures, aes(x = PC6, y = medianbc)) + geom_point() + geom_smooth(method = "lm", se = F)

# fit linear model predicting LDLHDL from dysbiosis score
# correct for PCs and covariates
summary(lm(LDLHDL ~ medianbc + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Age + Sex, data = corfeatures))

# scatter plot of this result ^
dysbio_scatter = ggplot(pheno, aes(x = medianbc, y = LDLHDL)) + 
  geom_point() + 
  xlab("Dysbiosis score") + 
  geom_smooth(method = "lm", se = F)
ggsave("../fig/scatterplot_dysbiosis_LDLHDL.png", height = 7, width = 7, dpi = 200)

# are distributions of dysbiosis score different across healthy and unhealthy?
# histogram^
nbin = 20
healthycol = "blue"
unhealthycol = "red"
ggplot() + 
  geom_histogram(data = pheno[pheno$healthy,], aes(x = medianbc), 
                 fill = healthycol, color = healthycol, alpha = 0.5, bins = nbin) +
  geom_histogram(data = pheno[!pheno$healthy,], aes(x = medianbc), 
                 fill = unhealthycol, color = unhealthycol, alpha = 0.5, bins = nbin)

# box plot^
ggplot(pheno) + geom_boxplot(aes(x = healthy, y = medianbc))

# is LDLHDL different in dysbiotic samples?
dysbio_boxplot = ggplot(pheno, aes(x = dysbiotic, y = LDLHDL)) + geom_boxplot() + 
  xlab("Dysbiotic") +
  stat_compare_means(method = "t.test", label.x = TRUE)
ggsave("../fig/boxplot_Disbiotic_LDLHDL.png", dysbio_boxplot, width = 5, height = 5, dpi = 200)
t.test(pheno[pheno$dysbiotic,"LDLHDL"], pheno[!pheno$dysbiotic, "LDLHDL"])

# pcs = prcomp(ab_scale, center = T, scale. = T)
# ggplot(data.frame(pcs$x)) + geom_point(aes(x = PC1, y = PC2, color = pheno$healthy))


# ggplot(pheno_unh, aes(x = median_BC, y = LDLHDL)) + geom_point()
# summary(lm(LDLHDL ~ median_BC, data = pheno_unh))
