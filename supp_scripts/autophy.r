#!/bin/R

message ("Beginning phyloseq run")



library(devtools)
library(lulu)
library(warppipe)
library(RColorBrewer)
library(gridExtra)
library(tidyr)
library(tibble)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(dada2)


message ("Making taxonomy file for phyloseq")
tax <-  as.matrix(read.table("phyloseqfile.txt"))
write.csv(tax, "phyloseq_taxonomyfile.csv") 


#import dataframe holding sample information (make sure Sample Names are the same as those attached to the OTU/seqtab table)
samdf<-read.csv("variables_phyloseq.csv")
head(samdf)
rownames(samdf) <- samdf$SampleID

head(LULUasv)
head(tax)

# Construct phyloseq object (straightforward from dada2 outputs)
LULUps <- phyloseq(otu_table(LULUasv, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))

LULUps

saveRDS(LULUps, file="LULUps.rds")

rel <- transform_sample_counts(LULUps, function(OTU) OTU/sum(OTU))

pdf("rel_abundance_nonfilt.pdf")
	plot_bar(rel, x="SampleID", fill="Class") + ggtitle("Relative Symbiont Abundances") + scale_fill_brewer(palette = "Set3")
dev.off()

ps.rel <- transform_sample_counts(rel, function(OTU) OTU/sum(OTU))
ps_melted <- psmelt(ps.rel)
ps <- ps_melted %>%
  filter(Abundance >= 0.01)

pdf("rel_abundance_filt.pdf")
  ggplot(data = ps, aes(x=SampleID, y=Abundance, fill = Class, order = Phylum)) +
  geom_bar(aes(), stat="identity", position="fill") #+
  #facet_wrap(~Compartment, scales="free_x")           ####can be edited to specify what chunk of data to show
dev.off()



save.image("phy.RData")
q("yes")
