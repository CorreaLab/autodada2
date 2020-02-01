#!/bin/R

message ("Checking if packages available, if not, installing")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("dada2", quietly = TRUE))
    BiocManager::install("dada2")

if (!requireNamespace("ShortRead", quietly = TRUE))
    BiocManager::install("ShortRead")

if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")

if (!requireNamespace("phyloseq", quietly = TRUE))
    install.packages("phyloseq")

message ("loading packages")

library(dada2); packageVersion("dada2"); citation("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")

wd <- getwd()

message ("Working directory =")
message (wd)

path <- getwd()
fns <- list.files(path)
message ("files present in working dir")
fns

message ("setting variables corresponding to read libraries present")
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files

message ("Forward reads")
fnFs
message ("Reverse reads")
fnRs

message ("Extracting samples names....")
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1) #the last number will select the field for renaming

sample.names 

message ("Setting paths to forward and reverse read files")
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

message ("Plotting quality of forward and reverse reads")
pdf("r1_qual_plot.pdf")
	plotQualityProfile(fnFs[c(1,2,3,4)])
	plotQualityProfile(fnFs[c(10,16,20,21)])
	dev.off()
pdf("r2_qual_plot.pdf")
	plotQualityProfile(fnRs[c(1,2,3,4)])
	plotQualityProfile(fnRs[c(10,11,12,13)])
	dev.off()

message ("Creating filtered read directoy and trimming using dada ish")
filt_path <- file.path(path, "trimmed_seqs")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(210,160), #leaves ~30bp overlap
              maxN=0, #DADA does not allow Ns
              maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
              truncQ=2,
              rm.phix=TRUE, #remove reads matching phiX genome
              matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)

message ("Now looking to define error rate")
setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

message ("Plotting error results")

pdf("fwd_error_plot.pdf")
	plotErrors(errF, nominalQ=TRUE) #some issues with C2G and G2C variants being underestimated, but not terrible
	dev.off()
pdf("rev_error_plot.pdf")
	plotErrors(errR, nominalQ=TRUE) #again, worse with G2C and C2G; a little T2G, but rest err on the side of being conservative (above red line)
	dev.off()

message ("Dereplicating now doe")
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

message ("Infering sequencing variants...")
setDadaOpt(BAND_SIZE=32)   ##can be edited
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

message ("Merging los reads")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

head(mergers[[2]])
summary((mergers[[1]]))

message ("Generating Sequence Table")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

pdf("variants.plot")
	plot(table(nchar(getSequences(seqtab)))) 
	dev.off()
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(229,271)] 
table(nchar(getSequences(seqtab2)))
dim(seqtab2)

message ("Removing chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
head(seqtab.nochim)
dim(seqtab.nochim)
write.csv(seqtab.nochim,file="seqtab.nochim_nonlabel.csv",row.names=TRUE,quote=FALSE)

message ("Reads kept after chimera removal =")
sum(seqtab.nochim)/sum(seqtab2)

message ("Tracking read statistics")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="ReadFilterStats_.csv",row.names=TRUE,quote=FALSE)

message ("Generating ASV file in workdir")
path <- "ASV_og_seqs.fasta"
uniquesToFasta(seqtab.nochim, path, ids = NULL, mode = "w", width = 20000)

write.csv(seqtab.nochim,file="asv_counts_table.csv",quote=F)

message ("making a copy of seqtab to save outputs needed for LULU")
seqtab.nochim.id <- seqtab.nochim
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim.id))))
colnames(seqtab.nochim.id)<-ids

write.csv(seqtab.nochim.id,file="seqtab.nochim_forLULU.csv",quote=F)

save.image("amptmp.RData")
q("yes")
