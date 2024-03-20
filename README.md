# install DECIPHER
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(DECIPHER)

install.packages("dada2")
install.packages("ggplot2")

### load packages ###
library(dada2)
library(DECIPHER) 
library(phangorn) 
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(data.table)
library(sapply)

# need to install MSA to continue here #
if (!requireNamespace("Bioconductor", quietly=TRUE))
  install.packages("Bioconductor")
BiocManager::install("msa")
install.packages("apply")

### specify the path to bring in the FASTQ files ###
path <- "/Users/charlotte.kenneally/desktop/GF_16s_1/"
list.files(path)


### sort files and extract sample names 
# sort files to ensure forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq"))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq"))

# extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`,1) 

# specify the full path to the fnFs and rnRs
### unnecessary, only run if wanted ###
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

### Plot average quality scores for the forward reads 
plotQualityProfile(fnFs[4:5]) # trimming at 230
# now reverse
plotQualityProfile(fnRs[4:5]) # the 4:5 is the 4th and 5th file, picked this to not chose the colon and cecum, we don't need to check every file if they are all part of the same run they likely have the same cutoffs 
# deciding to trim at 160 bp 

### create a new file path to store filtered and trimmed reads
filt_path <- file.path(path, "filtered") #place filtered files in filtered/subdirectory

#rename filtered files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

# Quality filtering and trimming 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,160),
                     trimLeft = 5, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) #on Mac set multithread=TRUE, on Windows set multithread=FALSE
head(out)

### estimate the error rate for DADA2 algorithm ###
#this command creates an error model that will be used by the DADA2 algorithm
#every batch of sequencing will have a different error rate
# algorithm starts with the assumption that the error rates are the maximum possible
# alternates error rate estimation and samples composition inference until they converge at a consistent solution 
errF <- learnErrors(filtFs, multithread=TRUE) # =TRUE even on Windows
errR <- learnErrors(filtRs, multithread=TRUE) 
# plot error rates for all possible base transitions, black line = observed error rates, red line = expected error rate under nominal definition of Q score, should show trend of freq. of errors decreasing as Q score increases (black should follow red line decently closely)
plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 

### dereplicate FASTQ files to speed up computation (simplifies it so you aren't just working through every single sample but only going through unique samples)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# name the derep-class objects by the sample names 
names(derepFs) <- sample.names
names(derepRs) <- sample.names

### Sample inference ###
# Apply core sequence-variant inference algorithm to reverse reads
dadaRs <- dada(derepRs, err=errR, multithread = TRUE)
# do forward 
dadaFs <- dada(derepFs, err=errR, multithread = TRUE)

### Merge denoised paired end reads ###
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

### Organize ribosomal sequence variants (RSVs) into a sequence table (same as OTU)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# view the length of all total RSVs
table(nchar(getSequences(seqtab)))

### Perform de novo chimera sequence detection and removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

### Calculate the proportion of non-chimeric RSVs
sum(seqtab.nochim)/sum(seqtab) # gives a percent of how much of my sample is not chimeras 

# calculating the proportion of non-chimeric RSV's to total RSV's
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))

#makes table to figure out how many sequences made it through the pipeline: filtered, denoised, merged, etc. 
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered5", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) 

### Assign Taxonomy using SILVA database 
# this is performed in two steps: this first one assigns phylum to genus
taxa <- assignTaxonomy(seqtab.nochim, "/Users/charlotte.kenneally/desktop/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
unname(head(taxa)) 

### Assign species (when possible)
system.time({taxa.plus <- addSpecies(taxa, "/Users/charlotte.kenneally/desktop/silva_species_assignment_v138.1.fa.gz", verbose=TRUE)
colnames(taxa.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unname(taxa.plus)})
head(taxa.plus)

### evaluate mock community 
# figure out what your mock sample is named by running the line below. Search for the exact name of the mock sample or samples used. 
rownames(seqtab.nochim)
###Change file name each time for each control community, move reference library into folder w samples, then run below sequence:
unqs.mock <- seqtab.nochim["ZYMO-A",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) #Drop ASV's absent in Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

# we don't need the ssRNAs that are in the ref sequence, so specify to only compare with the genomes folder
mock.ref <- getSequences(file.path(path, "ZymoBIOMICS.STD.refseq.v2_copy/Genomes"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those", sum(match.ref), "were exact matches to the expected reference sequences.\n")

unqs.mock <- seqtab.nochim["ZYMO-B",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) #Drop ASV's absent in Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path(path, "ZymoBIOMICS.STD.refseq.v2_copy/Genomes")) 
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those", sum(match.ref), "were exact matches to the expected reference sequences.\n")

unqs.mock <- seqtab.nochim["ZYMO-C",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) #Drop ASV's absent in Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path(path, "ZymoBIOMICS.STD.refseq.v2_copy/Genomes"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those", sum(match.ref), "were exact matches to the expected reference sequences.\n")

unqs.mock <- seqtab.nochim["ZYMO-D",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) #Drop ASV's absent in Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path(path, "ZymoBIOMICS.STD.refseq.v2_copy/Genomes"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those", sum(match.ref), "were exact matches to the expected reference sequences.\n")
# If 10/10 are identified, then the residual error rate for the DADA2 pipeline for this sample is 0%. (you want exact matches)
# Go back and edit your parameters to try and improve this

### some of these are bad, I am skipping over this for now to simply get to abundance, will come back to this later ###



## Align sequences -- extract sequences from DADA2 output ##
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences 
# Run sequence alignment (MSA) using DECIPHER 
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

### Construct phylogenetic tree ### note these next few lines to make the tree take about 10 min to run 
# change sequence alignment output into a phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
# create distance matrix
dm <- dist.ml(phang.align)
# perform neighbor joining 
treeNJ <-NJ(dm) # note, tip order != sequence order
# Internal maximum likelihood 
fit = pml(treeNJ, data = phang.align)
# Negative edges length changed to 0! 
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

#making sample data
samples.out <- rownames(seqtab.nochim)

#stringr extract to separate and group parts of the data set
# i am setting values to what I want groups to be named as, and telling it where to pull those samples from (pulling from the 'samples.out' value)

library(tidyverse)

library(stringr)

str_extract(samples.out, "1A-1")

### have to restart it bc error message str_extract is corrupt ### 
.rs.restartR() 

# grouping by round first 
round1 <- str_extract(samples.out, "1A|1B|1C|1D")
round2 <- str_extract(samples.out, "2A|2B|2C|2D")
round3 <- str_extract(samples.out, "3A|3B|3C|3D")
FMT_single_donor <- str_extract(samples.out, "Cr-NIH|MR1KO|SCITH3|C57BL")
FMT_pooled <- str_extract(samples.out, "Pooled")
Control <- str_extract(samples.out, "PBS|DNA|water|WATER") 

df <- data.frame(round1, round2, round3, FMT_single_donor, FMT_pooled, Control)
rownames(df) <- samples.out

### make phyloseq object ###
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
               tax_table(taxa.plus), phy_tree(fitGTR$tree), sample_data(df))
# view the phyloseq object in the console
ps


# keep only the samples I want to analyze -- take out the mock samples and DNA isolation
ps <- subset_samples(ps, Aerobic =="Yes", Anaerobic =="Yes", Live =="Yes", HeatKilled =="Yes", Wash == "Yes", 
                     PMA == "Yes", Control == "Yes")


# look at top 20 tax within the mock community 
### reason for this is if we want to remove taxa based on study prevalence, lets check the mock and the negative controls first 
ctrls <- c("ZYMO-A", "ZYMO-B", "ZYMO-C", "ZYMO-D", "WATER-A", "WATER-B", "WATER-C", "WATER-D")
ps.ctrls <- prune_samples(samples = ctrls, x = ps)
top20 <- names(sort(taxa_sums(ps.ctrls), decreasing=TRUE))[1:20]
ps.ctrls.top <- prune_taxa(taxa=top20, x=ps.ctrls)
sample_sums(ps.ctrls)

myranks = c("Genus", "Species")
mylabels = apply(tax_table(ps.ctrls.top)[, myranks], 1, paste, sep="", collapse="_")
tax_table(ps.ctrls.top) <- cbind(tax_table(ps.ctrls.top), genus_species=mylabels)

# plot that so it is easier to understand and visualize
plot_bar(transform_sample_counts(ps.ctrls.top, function(x) x/sum(x)), fill = "genus_species")+ylab("Proportion")
plot_bar(ps.ctrls.top, fill="genus_species") 

### abundances graphs
library(ggplot2)


library(reshape2)
df_long <- melt(df, id.vars = "Sample", variable.name = "Phyla")
ggplot(df_long, aes(x = Sample, y = value, fill = Phyla)) + 
  geom_bar(stat = "identity")
ggplot(data = df, aes(x = round1, y = value, fill = Phyla)) + 
  geom_bar(stat = "identity")

### try 2
prune.dat <- prune_taxa(taxa_sums(df) > 2, dat) #remove less than 2%
prune_taxa(taxa, ps)

dat.aglo = tax_glom(prune.dat, taxrank = "Genus")
dat.trans = transform_sample_counts(dat.aglo, function(x) x/sum(x))
prune.dat.two = prune_taxa(taxa_sums(dat.trans) > 0.01, dat.trans)
dat.dataframe = psmelt(prune.dat.two)
dat.agr = aggregate(Abundance~Site+Location+Genus, data=dat.dataframe, FUN=mean)

ggplot(dat.agr, aes(x=Site, y=Abundance, fill=Genus)) + geom_bar(stat="identity") + facet_grid(~Location, scale="free")

### start here after making object ###
#Assessing quality and read counts
ps3 %>% 
  psmelt() %>% 
  data.table() %>% 
  .[(Abundance > 0)] %>% 
  .[, .(Prevalence = uniqueN(Sample), 
        TotalCounts = sum(Abundance)), 
    by = .(Phylum, OTU)] %>% 
  ggplot(aes(Prevalence, TotalCounts, color = Phylum)) +
  geom_point(size = 2) +
  scale_y_log10() +
  ggtitle("Total Counts v. Prevalence", "Each point is an ASV.") 

#plotting cumulative distribution of read counts vs prevalence
dataPlotCumPrev <-
  ps %>% 
  psmelt() %>% 
  data.table() %>% 
  .[(Abundance > 0)] %>% 
  .[, .(Prevalence = uniqueN(Sample), 
        TotalCounts = sum(Abundance)), 
    by = .(Phylum, OTU)] %>% 
  setorder(Prevalence, TotalCounts) %>% 
  .[, CumulativeCounts := cumsum(TotalCounts)]
# plot
install("ggrepel")
library(ggrepel)
dataPlotCumPrev %>% 
  ggplot(aes(Prevalence, CumulativeCounts), color = "gray") +
  geom_point(size = 1) +
  geom_point(
    data = dataPlotCumPrev[(Prevalence >= 20L)],
    mapping = aes(color = Phylum),
    size = 2) +
  geom_text_repel(
    size = 2,
    mapping = aes(label = OTU), 
    data = dataPlotCumPrev[(Prevalence >= 26)]) +
  scale_y_log10() +
  ggtitle("Cumulative Read Count v. Prevalence", "Each point is an OTU.")


# cumulative read count vs prevalence, cumulative number of OTUs
setorder(dataPlotCumPrev, Prevalence, TotalCounts)
dataPlotCumPrev[, CumulativeOtus := .I]
# plot
dataPlotCumPrev %>% 
  ggplot(aes(Prevalence, CumulativeOtus), color = "gray") +
  geom_point(size = 1) +
  geom_point(
    data = dataPlotCumPrev[(Prevalence >= 20L)],
    mapping = aes(color = Phylum),
    size = 2) +
  geom_text_repel(
    size = 2,
    mapping = aes(label = OTU), 
    data = dataPlotCumPrev[(Prevalence >= 23L)]) +
  ylab("Cumulative Number of OTUs") +
  ggtitle("Cumulative Read Count v. Prevalence", "Each point is an OTU.")
#converting to data table and then getting total reads
sdt <- data.table(as(sample_data(ps), "data.frame"),
                  TotalReads = sample_sums(ps), keep.rownames = TRUE)
#graphing total reads
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram(bins = 100) + ggtitle("Sequencing Depth of all samples")
pSeqDepth 
pSeqDepth + facet_wrap(round1~round2) 

view((otu_table(ps)))
# remove reads below 600 - this gets rid of dna isolation, PBS, and water negs
ps600 <- prune_samples(sample_sums(ps)>=600, ps)
sample_names(ps600)
view(otu_table(ps600))
#pruning out singletons
ps600 <- prune_taxa(taxa_sums(ps) > 1, ps)

#saving new phyloseq object with pruned data
saveRDS(ps600, "ps600.rds")
ps_derep <- readRDS("ps600.rds")


### bring in sample data frame that was made in excel 
### in excel -- save file as .csv and here direct the path to that file 
library(readr)
samdat <- read.csv("~/Desktop/GF_16s_1/GF16s_samdat.csv") %>% as.data.frame()
View(samdat)

View (otu_table(ps600)) 


OTU600 <- as.data.table(otu_table(ps600))
#OTU600 <- OTU600[1:285,]
samdat600 <- sample_data(ps600)#[1:285,]
tax600 <- as.data.frame(tax_table(ps600))

### save as an excel file 
write.csv(OTU600, "/Users/charlotte.kenneally/desktop/GF_16s_1/OTU600_2.csv")
write.csv(samdat600, "/Users/charlotte.kenneally/desktop/GF_16s_1/samdat600_2.csv")
write.csv(tax600, "/Users/charlotte.kenneally/desktop/GF_16s_1/tax-600 _2.csv")

