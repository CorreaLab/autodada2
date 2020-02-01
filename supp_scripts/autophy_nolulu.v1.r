#!/bin/R

message ("Beginning phyloseq run")

load(file= "amptmp.RData")

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

head(seqtab.nochim)
head(tax)

# Construct phyloseq object (straightforward from dada2 outputs)
nonLULUps <- phyloseq(otu_table(nonLULUasv, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))

nonLULUps
saveRDS(nonLULUps, file="nonLULUps.rds")


rel <- transform_sample_counts(nonLULUps, function(OTU) OTU/sum(OTU))
pdf("rel_abundance_nonfilt.plot")
	plot_bar(rel, x="SampleID", fill="Class") + ggtitle("Relative Symbiont Abundances") + scale_fill_brewer(palette = "Set3")
	dev.off()


ps.rel <- transform_sample_counts(rel, function(OTU) OTU/sum(OTU))
ps_melted <- psmelt(ps.rel)
ps <- ps_melted %>%
  filter(Abundance >= 0.01)

pdf("rel_abundance_filt.plot")
  ggplot(data = ps, aes(x=SampleID, y=Abundance, fill = Order, order = Phylum)) +
  geom_bar(aes(), stat="identity", position="fill") #+
  #facet_wrap(~Compartment, scales="free_x")
dev.off()


save.image("amptmp.RData")
q("yes")
