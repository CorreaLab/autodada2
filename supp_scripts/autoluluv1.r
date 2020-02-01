#!/bin/R
#title: "LULU"
#author: "Lauren Howe-Kerr & Alex Veglia"
#date: "January 2020"
#output: pdf_document

#Need these two inputs for LULU:

#*OTU fasta file input name:* ASV_LULU_labs.fasta

#*DADA output file name:* seqtab.nochim_LULU.csv

#Necessary pre-steps after producing ASV table and associated fasta file:
#Produce a match list using BLASTn


message ("Checking for packages, installing if not done already")

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

library(devtools)

if (!requireNamespace("lulu", quietly = TRUE))
    install_github("tobiasgf/lulu")

library(lulu)

if (!requireNamespace("warppipe", quietly = TRUE))
    install_github("seankross/warppipe")

library(warppipe)

if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr")

if (!requireNamespace("RColorBrewer", quietly = TRUE))
    install.packages("RColorBrewer")

library(RColorBrewer)

if (!requireNamespace("gridExtra", quietly = TRUE))
    BiocManager::install("gridExtra")

library(gridExtra)

if (!requireNamespace("tidyr", quietly = TRUE))
     BiocManager::install("tidyr")

library(tidyr)

if (!requireNamespace("tibble", quietly = TRUE))
    BiocManager::install("tibble")

library(tibble)

if (!requireNamespace("ShortRead", quietly = TRUE))
    BiocManager::install("ShortRead")

library(ShortRead)

if (!requireNamespace("ggplot2", quietly = TRUE))
    BiocManager::install("ggplot2")

library(ggplot2)

if (!requireNamespace("phyloseq", quietly = TRUE))
    BiocManager::install("phyloseq")

library(phyloseq)

library(dplyr)
library(dada2)
message ("Packages checked/installed")

#manually add "SampleID" to header of column with sample IDs

alldat<-read.csv("seqtab.nochim_forLULU.csv") 
head(alldat)
alldat <- alldat %>% rename(SampleID=X)

match_list<-read.table("match_list.txt")
head(match_list)

#Reformat ASV table to desired LULU format
rownames(alldat)<-alldat$SampleID
ASVs<-data.frame(t(alldat[,-1])) #to remove the sequence ID column
head(ASVs)

#Now, run the LULU curation
curated_result <- lulu(ASVs, match_list) #minimum_relative_cooccurence=0.95 (this is default, can modify for more conservative values)
summary(curated_result)

curated_result$minimum_match
curated_result$minimum_relative_cooccurence

#make the results into a dataframe
alldat<-cbind(data.frame(curated_result$curated_table))
head(alldat)

alldat <- rownames_to_column(alldat, var = "Sq.id")

#make OTU fasta file into a dataframe (for merging with dataframe above)
ASV_ids <- seq_tbl("ASV_lulu_labs.fasta")
head(ASV_ids)

#merging the ASV dataframe with the big main dataframe to get the sequences attached
alldat2 <- merge(alldat, ASV_ids, by.x = 'Sq.id', by.y='Description')
#now make the Sequence the first column and delete the Sq.id column
new_df <- alldat2 %>%
  select(Sequence, everything()) %>%
  select(-2)

#transform it for Assign Taxonomy function
rownames(new_df)<-new_df$Sequence
alldat_final<-cbind(data.frame(t(new_df[,-1])))


message ("Make it into a matrix so its formatted correctly for assigntaxonomy")
LULUasv <- as.matrix(alldat_final)
head(LULUasv)


message ("Generating new fasta with LULU asvs")
path='LULUasv.fasta'
uniquesToFasta(LULUasv, path, ids = NULL, mode = "w", width = 20000)

##################################MAKING COPIES OF OUTPUT FILES########################################################################
#making a new LULUotu with "ids" so that I have a copy with the sequences as headers (LULUotu_2) and a copy with short sq names as headers (LULUotu_2ids)
LULUasv_ids <- LULUasv

#then, rename output table and write it out
ids <- paste0("sq", seq(1, length(colnames(LULUasv_ids))))
colnames(LULUasv_ids)<-ids

head(LULUasv)
head(LULUasv_ids)

write.csv(LULUasv_ids,file="LULUasv_ids_out.csv",quote=F)
write.csv(LULUasv,file="LULUasv_out.csv",quote=F)


save.image("lulu.RData")
q("yes")
