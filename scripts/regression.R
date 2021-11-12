library(glmnet)

base_ab_pheno = read.table("../data/baseline_genusAbundance_pheno.txt", header = T)
base_abund = read.table("../data/baseline_genusAbundance.txt")

# correlate each genus abundance with LDL
corrs <- apply(base_abund, 2, function (x) cor(x, base_ab_pheno$LDL))
sort(corrs)

# fit a standard linear model 
lm1 <- lm(base_ab_pheno$LDL ~ . , data = base_abund)

# fit a lasso regression
glm1 <- glmnet(x = base_abund, y = base_ab_pheno$LDL, alpha = 1)

ypred = predict(glm1, newx = as.matrix(base_abund))