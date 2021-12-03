# Use PCA to visualize data and look for patterns
library(ggpubr)
library(dplyr)
library(ggplot2)
library(glmnet)
library(tidyverse)


# file paths are relative to user's folders

base_ab_pheno <- read.table("~/microbiome-cs226/data/baseline_genusAbundance_pheno.txt", header=TRUE)
base_abund <- read.table("~/microbiome-cs226/data/baseline_genusAbundance.txt", header=TRUE)


###
# RUN PCA ON BASELINE ABUNDANCE DATA
###

# Run PCA, scaling recommended
pca_abund <- prcomp(base_abund, scale=TRUE)

PCs_df <- data.frame(pca_abund$x)

### 
# PLOT PCS 
###

# Get data of interest for plotting
LDL <- base_ab_pheno$LDL

# Bin LDL levels into 3 categories 
# (based on typical levels https://www.mayoclinic.org/tests-procedures/cholesterol-test/about/pac-20384601):
# LDL < 100 mg/dL, LDL in the 100-129 range mg/dL, and LDL 130+ mg/dL. 
# (A bin of 130-160 is typical, but with this data leaves only 5 points in the 160+ range)
v_ldl <- c(100, 130, 160)
LDL_binned <- as.character(findInterval(LDL, v_ldl))

### LDL PLOT
ggplot(data = PCs_df) + 
  geom_point(aes(x = PC1, y = PC2, color = LDL_binned)) + 
  scale_color_brewer(palette ='Set1') +
  theme_bw() 

# Run similar approach for CHOL
# Same Mayo clinic data for bins
# Bin 1: CHOL < 200 mg/dL, Bin 2: CHOL 200-239 mg/dL, Bin 3: CHOL 240+ mg/dL
CHOL <- base_ab_pheno$CHOL
v_chol <- c(200, 240)
CHOL_binned <- as.character(findInterval(CHOL, v_chol))

### CHOL PLOT
# color points based on cholesterol levels
ggplot(data = PCs_df) + 
  geom_point(aes(x = PC1, y = PC2, color = CHOL_binned)) + 
  scale_color_brewer(palette ='Set1') +
  theme_bw() 


# A1C coloring
# https://www.cdc.gov/diabetes/managing/managing-blood-sugar/a1c.html#:~:text=A%20normal%20A1C%20level%20is,for%20developing%20type%202%20diabetes.
# percentage of your red blood cells that have sugar-coated hemoglobin
# Bin 1: A1C < 5.7%, Bin 2: A1C 5.7-6.3%, Bin 3: 6.4+%
A1C <- base_ab_pheno$A1C
v_a1c <- c(5.7, 6.4)
A1C_binned <- as.character(findInterval(A1C, v_a1c))

### A1C PLOT
ggplot(data = PCs_df) + 
  geom_point(aes(x = PC3, y = PC4, color = A1C_binned)) + 
  scale_color_brewer(palette ='Set1') +
  theme_bw() 


### 
# RUN PCA ON ALL GENUS ABUNDANCE DATA ACROSS TIMEPOINTS
### 

# Import the data across all the timepoints
ab_pheno_all = read.table("~/microbiome-cs226/data/allTimePoints_abundance_pheno_visitInfo.txt", header = T, sep=',')
ab_pheno_all = cbind(select(ab_pheno_all, starts_with("genus")),
                     ab_pheno_all[c("SampleID", "SubjectID", "CollectionDate", "LDLHDL", "LDL", "CHOL", "A1C")])
ab_pheno_all = na.omit(ab_pheno_all)
pheno_all = ab_pheno_all[c("SampleID", "SubjectID", "CollectionDate", "LDLHDL", "LDL", "CHOL", "A1C")]
ab_all = select(ab_pheno_all, starts_with("genus"))

# Run PCA
pca_allTimepoints <- prcomp(ab_all, scale=TRUE)
PCs_allTimepoints_df <- data.frame(pca_allTimepoints$x)

# Bin LDLHDL for all timepoints based on levels
# this isn't an academic source persay but is useful
# https://www.healthline.com/health/cholesterol-ratio
LDLHDL_all <- pheno_all$LDLHDL
v_ldlhdl <- c(3.5, 4, 5)
LDLHDL_all_binned <- findInterval(LDLHDL_all, v_ldlhdl)
for (i in 1:length(LDLHDL_all_binned)){
  if (LDLHDL_all_binned[i] == 0)
    LDLHDL_all_binned[i] = "LDL/HDL < 3.5"
  if (LDLHDL_all_binned[i] == 1)
    LDLHDL_all_binned[i] = "LDL/HDL 3.5 - 4"
  if (LDLHDL_all_binned[i] == 2)
    LDLHDL_all_binned[i] = "LDL/HDL 4 - 5"
  if (LDLHDL_all_binned[i] == 3)
    LDLHDL_all_binned[i] = "LDL/HDL > 5"
}

### LDLHDL PLOT
myscatterplot <- ggplot(data = PCs_allTimepoints_df) + 
  geom_point(aes(x = PC1, y = PC2, color = LDLHDL_all_binned)) + 
  scale_color_brewer(palette ='Set1') +
  theme_bw() +
  labs(color='LDL HDL Cholesterol Ratio')

## ggsave
# dpi 200
ggsave("~/Desktop/ldlhdlscatter.png", myscatterplot, dpi=200)

# Bin LDL for all timepoints (can use earlier bin vector)
LDL_all <- pheno_all$LDL
LDL_all_binned <- as.character(findInterval(LDL_all, v_ldl))
### LDL PLOT
ggplot(data = PCs_allTimepoints_df) + 
  geom_point(aes(x = PC1, y = PC2, color = LDL_all_binned)) + 
  scale_color_brewer(palette ='Set1') +
  theme_bw() 

# Bin CHOL for all timepoints (can use earlier bin vector)
CHOL_all <- pheno_all$CHOL
CHOL_all_binned <- as.character(findInterval(CHOL_all, v_chol))
### CHOL PLOT
ggplot(data = PCs_allTimepoints_df) + 
  geom_point(aes(x = PC1, y = PC6, color = CHOL_all_binned)) + 
  scale_color_brewer(palette ='Set1') +
  theme_bw() 

# Bin A1C for all timepoints (can use earlier bin vector)
### A1C PLOT
A1C_all <- pheno_all$A1C
A1C_all_binned <- as.character(findInterval(A1C_all, v_a1c))
ggplot(data = PCs_allTimepoints_df) + 
  geom_point(aes(x = PC1, y = PC2, color = A1C_all_binned)) + 
  scale_color_brewer(palette ='Set1') +
  theme_bw() 


### 
# PC CORRELATIONS
###

# See which PCs correlate well with our data
LDLHDL_PC_cor = cor(PCs_allTimepoints_df, LDLHDL_all)
A1C_PC_cor = cor(PCs_allTimepoints_df, A1C_all)
CHOL_PC_cor = cor(PCs_allTimepoints_df, CHOL_all)
LDL_PC_cor = cor(PCs_allTimepoints_df, LDL_all)

# I noticed a pattern of higher correlated PCs for the different measurement types
# List of highly correlated PCs across data
PC_cors = c("PC2", "PC6", "PC11", "PC15", "PC16", "PC34", "PC35")


###
# Try running some regression models using these PCs as X data
###
data_PCs = as.matrix(cbind(PCs_allTimepoints_df[,PC_cors], A1C_all, CHOL_all, LDLHDL_all, LDL_all))

# Modeled after Regression code in model_comparison.R
set.seed(222)
train_rows = sample(1:nrow(data_PCs), round(nrow(data_PCs)*.6))
x.train = data_PCs[train_rows, 1:7]
y.train = data_PCs[train_rows, 8:11]

x.test = data_PCs[-train_rows, 1:7]
y.test = data_PCs[-train_rows, 8:11]

# Linear regression with PCs for all response variables
PC.lm = lm(y.train ~ . , data = as.data.frame(x.train))

ypred.PC <- predict(PC.lm, newdata = as.data.frame(x.test))
mse.PC <- mean((y.test - ypred.PC)^2)

# Predict singular y vectors (responses) using Elastic Net Regression
# LDLHDL Model
y_ldlhdl.train = data_PCs[train_rows, "LDLHDL_all"]
y_ldlhdl.test = data_PCs[-train_rows, "LDLHDL_all"]

PC_ldlhdl.enet = cv.glmnet(x.train, y_ldlhdl.train, type.measure="mse", alpha=0.5,family="gaussian")
fit_ldlhdl.enet <- glmnet(x.train, y_ldlhdl.train, family="gaussian", alpha=.5)
plot(fit_ldlhdl.enet, xvar="lambda")
plot(PC_ldlhdl.enet, main="Elastic Net")

ldlhdl_ypred <- predict(PC_ldlhdl.enet, s=PC_ldlhdl.enet$lambda.1se, newx=x.test)
mse_ldlhdl <- mean((y_ldlhdl.test - ldlhdl_ypred)^2)

((y_ldlhdl.test - ldlhdl_ypred)^2)/(y_ldlhdl.test )


### VARIANCE EXPLAINED PLOT

# pca_allTimepoints is our PCA
pcsummary <- summary(pca_allTimepoints)
variance <- data.frame(PC =  names(pcsummary$importance[2,1:10]), Variance = pcsummary$importance[2,1:10])
level_order = c(variance$PC)

varianceplot <- ggplot(data = variance, aes(x=factor(PC, levels=level_order), y=Variance, group=1)) +
  geom_line() +
  geom_point() +
  xlab("PC") +
  ylab("Variance Explained")

varianceplot

ggsave("~/Desktop/variance_explained.png", varianceplot, dpi=200)

