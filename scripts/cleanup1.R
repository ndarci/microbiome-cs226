library(dplyr)
library(tidyverse)
library(reshape)

ct <- read.delim("~/cs226 microbiome/clinical_tests.txt", header=T, stringsAsFactors=T)
# View(ct)
gut16s <- read.delim("~/cs226 microbiome/gut_16s_abundance.txt", header=T, stringsAsFactors=T)
# View(gut16s)
subjects <- read.delim("~/cs226 microbiome/subjects.txt", header=T, stringsAsFactors=T)
# View(subjects)
visits <- read.delim("~/cs226 microbiome/visits.txt", header=T, stringsAsFactors=T)
# View(visits)

# unique(gut16s$SampleID) %in% unique(visits$VisitID)
colnames(visits)[1] <- "SampleID"

df <- left_join(gut16s,visits,by="SampleID")
df1 <- df[,c(1, 98,99,2:97,100:104)]

df2 <- df1[!is.na(df1$CollectionDate),]

colnames(ct)[1] <- "SampleID"
# trim last few cols
ct <- ct[,c(1:52)]
# glimpse(ct)
# there are quite a few NAs / missing vals here
rowSums(is.na(ct))
colSums(is.na(ct))
# 943/969 missing values for insulin... hell yea

df3 <- left_join(df2, ct, by = "SampleID")

df4 <- df3[,-c(100:104)]


baselines <- df4 %>%
  group_by(SubjectID) %>%
  mutate(
    basetime = min(CollectionDate, na.rm = T)
  ) %>%
  arrange(SubjectID)

baseline_time <- baselines %>% select(SubjectID,basetime) %>% unique()

df6 <- merge(baseline_time, df4, by.x = c("SubjectID","basetime"), by.y = c("SubjectID","CollectionDate"), all.x = TRUE)
df7 <- cbind(df6[,c("SubjectID","basetime","SampleID","A1C", "CHOL", "CHOLHDL", "HDL", "LDL", "LDLHDL")],df6 %>% select(starts_with("genus")))
df7 <- na.omit(df7)

genus <- (df7 %>% select(starts_with("genus")))
notgenus <- (df7 %>% select(-starts_with("genus")))

# corrs
corrs <- apply(genus, 2, function (x) stats::cor(x, df7$LDL))
sort(corrs)

lm1 <- lm(df7$LDL ~ . , data = genus)

library(glmnet)

glm1 <- glmnet(x = genus, y = df7$LDL,alpha = 1)

ypred = predict(glm1, newx = as.matrix(genus))

# setwd("~/cs226 microbiome")
save(df7, file = "df7.csv")


