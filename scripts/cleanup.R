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


# some of the observations in visits had sample ids that didnt align with gut 16s sample id
# like this fuckin guy
df1[grep("ZOZOW1T", df1$SampleID), c(1,2,3)]
# only about half have perfect matches, otherwise NAs
# NA values mean that match was not found for sample ID
# which means missing collection date (time var) and subjID

# itll be pretty easy to link the subject ID with the sample IDs
# but idk if anything can be done about the missing collection dates...
# other than delete them.

# for now im dropping the samples with no collection time, ie no match in SampleID btwn gut16s and visits
sum(is.na(df1$CollectionDate))
# not too many....
# bye for now
df2 <- df1[!is.na(df1$CollectionDate),]

# entroducing.... clinical tests data
colnames(ct)
colnames(ct)[1] <- "SampleID"
# trim last few cols
ct <- ct[,c(1:52)]
# glimpse(ct)
# there are quite a few NAs / missing vals here
rowSums(is.na(ct))
colSums(is.na(ct))
# 943/969 missing values for insulin... hell yea

# our final product, for now:
df3 <- left_join(df2,ct, by = "SampleID")

# View(df3)
# how many subjects end up in df3?
length(unique(df3$SubjectID))
rowSums(is.na(df3))
colSums(is.na(df3))

colnames(df3)

# focus only on hmp subjects...?
visits$SubStudy %>% summary()
df1$SubStudy %>% summary()
df2$SubStudy %>% summary()
df3$SubStudy %>% summary()


# the gut16 data already provides a solid amount of predictors... we need to think of an appropriate response.