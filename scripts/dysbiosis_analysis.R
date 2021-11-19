# test if gut microbe dysbiosis is connected to high cholesterol
library(ggplot2)
library(ggpubr)
library(dplyr)

# setwd("~/src/microbiome-cs226/scripts/")

# import data
ab = read.table("../data/baseline_genusAbundance.txt")
ab_pheno = read.table("../data/baseline_genusAbundance_pheno.txt", header = T)
pheno = ab_pheno[1:9]

ab_pheno_all = read.table("../data/allTimePoints_abundance_pheno_visitInfo.txt", header = T, sep = ',')
ab_pheno_all = cbind(select(ab_pheno_all, starts_with("genus")), 
                     ab_pheno_all[c("SampleID", "SubjectID", "CollectionDate", "LDLHDL")])
ab_pheno_all = na.omit(ab_pheno_all)
pheno_all = ab_pheno_all[c("SampleID", "SubjectID", "CollectionDate", "LDLHDL")]
ab_all = select(ab_pheno_all, starts_with("genus"))
rownames(ab_all) = ab_pheno_all$SampleID

# summary(lm(ab_pheno_all$LDLHDL ~ ., data = ab_all))
# 
# pcs_all = prcomp(ab_all)
# ggplot(data.frame(pcs_all$x)) + geom_point(aes(x = PC1, y = PC2, 
#                                                color = ab_pheno_all$LDLHDL > mean(ab_pheno_all$LDLHDL)))

# function to calculate bray-curtis dissimilarity between two samples (assumes they are counted as proportions)
# todo: remove samples taken from the same subject
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

# fit linear model predicting LDLHDL from dysbiosis score
summary(lm(LDLHDL ~ medianbc, data = pheno))

# scatter plot of this result ^
ggplot(pheno, aes(x = medianbc, y = LDLHDL)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F)

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
ggplot(pheno, aes(x = dysbiotic, y = LDLHDL)) + geom_boxplot() + 
  stat_compare_means(method = "t.test", label.x = TRUE)
t.test(pheno[pheno$dysbiotic,"LDLHDL"], pheno[!pheno$dysbiotic, "LDLHDL"])

# pcs = prcomp(ab_scale, center = T, scale. = T)
# ggplot(data.frame(pcs$x)) + geom_point(aes(x = PC1, y = PC2, color = pheno$healthy))


# ggplot(pheno_unh, aes(x = median_BC, y = LDLHDL)) + geom_point()
# summary(lm(LDLHDL ~ median_BC, data = pheno_unh))
