library(dplyr)
library(tidyverse)
library(reshape)

# paths are relative to scripts/ folder

# import data
ct <- read.delim("../data/clinical_tests.txt", header=T, stringsAsFactors=T)
# View(ct)
gut16s <- read.delim("../data/gut_16s_abundance.txt", header=T, stringsAsFactors=T)
# View(gut16s)
subjects <- read.delim("../data/subjects.txt", header=T, stringsAsFactors=T)
# View(subjects)
visits <- read.delim("../data/visits.txt", header=T, stringsAsFactors=T)
# View(visits)

# unique(gut16s$SampleID) %in% unique(visits$VisitID)

# merge abundance and visit data
colnames(visits)[1] <- "SampleID"
ab_visit <- left_join(gut16s,visits,by="SampleID")

# remove event notes and missing collection date
ab_visit_clean <- ab_visit[,c(1:104)]
ab_visit_clean <- ab_visit_clean[!is.na(ab_visit_clean$CollectionDate),]

# add clinical test data
colnames(ct)[1] <- "SampleID"
# trim last few cols
ct <- ct[,c(1:52)]
# glimpse(ct)
# there are quite a few NAs / missing vals here
rowSums(is.na(ct))
colSums(is.na(ct))
# 943/969 missing values for insulin... hell yea

# merge with abundance and visit data
alldata_clean <- left_join(ab_visit_clean, ct, by = "SampleID")
# remove event notes
alldata_clean <- alldata_clean[,-c(100:104)]

# select baseline samples for each individual
# define baseline as minimum collection date
baselines <- alldata_clean %>%
  group_by(SubjectID) %>%
  mutate(
    basetime = min(CollectionDate, na.rm = T)
  ) %>%
  arrange(SubjectID)
baseline_time <- baselines %>% select(SubjectID,basetime) %>% unique()

# base_ab_pheno has complete abundance data and liver lab data for all baseline individuals
base_ab_pheno <- merge(baseline_time, alldata_clean, by.x = c("SubjectID","basetime"), by.y = c("SubjectID","CollectionDate"), all.x = TRUE)
base_ab_pheno <- cbind(base_ab_pheno[,c("SubjectID","basetime","SampleID","A1C", "CHOL", "CHOLHDL", "HDL", "LDL", "LDLHDL")], base_ab_pheno %>% select(starts_with("genus")))
base_ab_pheno <- na.omit(base_ab_pheno)

# get a dataframe with just the abundance variables
base_abund <- (base_ab_pheno %>% select(starts_with("genus")))
rownames(base_abund) <- base_ab_pheno$SampleID

write.table(base_ab_pheno, "../data/baseline_genusAbundance_pheno.txt", quote = F, row.names = F)
write.table(base_abund, "../data/baseline_genusAbundance.txt", quote = F, row.names = T)
write.table(alldata_clean, "../data/allTimePoints_abundance_pheno_visitInfo.txt", quote = F, row.names = F)


