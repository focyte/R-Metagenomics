# Load the required libraries
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(dplyr)
library(ggpubr)
library(patchwork)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##              ██████╗  █████╗ ██████╗  █████╗ ██████╗ 
##              ██╔══██╗██╔══██╗██╔══██╗██╔══██╗╚════██╗
##              ██║  ██║███████║██║  ██║███████║ █████╔╝
##              ██║  ██║██╔══██║██║  ██║██╔══██║██╔═══╝ 
##              ██████╔╝██║  ██║██████╔╝██║  ██║███████╗
##              ╚═════╝ ╚═╝  ╚═╝╚═════╝ ╚═╝  ╚═╝╚══════╝

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


path <- "./"
list.files(path)

# Create a list of filenames where files have the extension ".fastq.gz"
file.names <- sort(list.files(path, pattern=".fastq.gz", full.names = TRUE))
# Extract sample names
sample.names <- sub("\\.fastq\\.gz$", "", basename(file.names))

# Check the quality of reads for the first two samples and save the plot
pdf(file="read_quality.pdf")
plotQualityProfile(file.names[1:2])
dev.off()

# Create filtered reads subdirectory
filtered.reads <- file.path(path, "filtered", paste0(sample.names, "_filtered.fastq.gz"))
names(filtered.reads) <- sample.names


# Check the QualityProfile plot and decide where the reads need truncating
out <- filterAndTrim(file.names, filtered.reads, 
                     rev = NULL,
                     filt.rev = NULL,
                     minLen=400,
                     truncLen=400,
                     maxN=0, 
                     maxEE=2, 
                     truncQ=2, 
                     rm.phix=TRUE,
                     compress=TRUE, 
                     multithread=FALSE)
head(out)

# Learn the error rates and save a plot of the results
error.rate <- learnErrors(filtered.reads, multithread=TRUE)
pdf(file="error_rate.pdf")
plotErrors(error.rate, nominalQ=TRUE)
dev.off()

# Apply the core sample inference algorithm to the filtered and trimmed sequence data.
# Extra steps added to account for 454 pyrosequencing technology used
dadaReads <- dada(filtered.reads, 
                  err=error.rate, 
                  multithread=TRUE, 
                  HOMOPOLYMER_GAP_PENALTY=-1, 
                  BAND_SIZE=32) # Accounts for high number of indels in 454 sequencing
dadaReads[[1]]

# Construct an amplicon sequence variant table (ASV) table 
# ASV table is a higher-resolution version of an OTU table 
seqtab <- makeSequenceTable(dadaReads)
dim(seqtab)

# Inspect the distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras from the sequence table
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaReads, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
write.csv(track,"./DADA2_processing_results.csv", row.names = TRUE)
head(track)

# Assign taxonomy using train data from: https://zenodo.org/records/1172783
taxa <- assignTaxonomy(seqtab.nochim, "./silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# Inspect the taxonomic assignments
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa,"./taxa.csv", row.names = TRUE)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##    ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ███████╗███████╗ ██████╗ 
##    ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔════╝██╔════╝██╔═══██╗
##    ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║███████╗█████╗  ██║   ██║
##    ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║╚════██║██╔══╝  ██║▄▄ ██║
##    ██║     ██║  ██║   ██║   ███████╗╚██████╔╝███████║███████╗╚██████╔╝
##    ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚══════╝╚══════╝ ╚══▀▀═╝

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Construct a simple sample data.frame
samples.out <- rownames(seqtab.nochim)

# Extract subject and Sampletype
subject <- sapply(strsplit(samples.out, "_"), `[`, 2) 
sampletype <- sapply(strsplit(samples.out, "_"), `[`, 1)

# Create the data frame
samdf <- data.frame(Subject=subject, Sampletype=sampletype)
samdf$Sampletype[samdf$Sampletype == "P"] <- "Plaque"
samdf$Sampletype[samdf$Sampletype == "S"] <- "Saliva"

# Set rownames to match the original sample names
rownames(samdf) <- samples.out

# Construct a phyloseq object directly from the dada2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Identify and remove taxa with no NA values in their taxonomic classification
# valid_taxa <- taxa_names(ps)[!apply(is.na(tax_table(ps)), 1, any)]
# ps <- prune_taxa(valid_taxa, ps)
# ps

# Estimate richness
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

# Add sample data to alpha_div for grouping information
alpha_div$Sampletype <- sample_data(ps)$Sampletype

# Shannon Diversity Plot
shannon_plot <- ggplot(alpha_div, aes(x = Sampletype, y = Shannon, color = Sampletype)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  stat_compare_means(method = "wilcox.test", label = "p.format") + 
  labs(title = "Shannon Diversity by Sample Type",
       x = "Sample Type",
       y = "Shannon Diversity") +
  theme_bw()

# Simpson Diversity Plot
simpson_plot <- ggplot(alpha_div, aes(x = Sampletype, y = Simpson, color = Sampletype)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  stat_compare_means(method = "wilcox.test", label = "p.format") + 
  labs(title = "Simpson Diversity by Sample Type",
       x = "Sample Type",
       y = "Simpson Diversity") +
  theme_bw()

# Combine plots and save to PDF
combined_plot <- shannon_plot + simpson_plot + plot_layout(ncol = 2)
pdf(file = "Richness_Shannon_Simpson_with_stats.pdf", width = 12, height = 6)
print(combined_plot)
dev.off()

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

pdf(file="Ordination.pdf")
plot_ordination(ps.prop, ord.nmds.bray, color="Sampletype", title="Bray NMDS")
dev.off()

# Normalize data to percentages for each sample
ps.percent <- transform_sample_counts(ps, function(x) 100 * x / sum(x))
write.csv(ps.percent@otu_table,"./percent_otu_table.csv", row.names = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##                        Plotting of the data

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Select the top 20 sequences by the sum of taxa and normalize the counts


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

pdf(file="Subject_Phylum_Top20.pdf")
plot_bar(ps.top20, x="Sample", fill="Phylum") + facet_wrap(~Sampletype, scales="free_x")
dev.off()

pdf(file="Subject_Family_Top20.pdf")
plot_bar(ps.top20, x="Sample", fill="Family") + facet_wrap(~Sampletype, scales="free_x")
dev.off()

pdf(file="Subject_Genus_Top20.pdf")
plot_bar(ps.top20, x="Sample", fill="Genus") + facet_wrap(~Sampletype, scales="free_x")
dev.off()


# Merge samples by Sampletype and set the factor levels
mergedPs <- merge_samples(ps, "Sampletype")
sample_data(mergedPs)$Sampletype <- factor(
  c("Plaque", "Saliva"),
  levels = c("Plaque", "Saliva")
)
write.csv(mergedPs@otu_table,"./mergedOTU.csv", row.names = TRUE)

# Transform OTU counts to percentages
mergedPs_percent <- transform_sample_counts(mergedPs, function(x) (x / sum(x)) * 100)
sample_data(mergedPs_percent)$Sampletype <- factor(
  c("Plaque", "Saliva"),
  levels = c("Plaque", "Saliva")
)
write.csv(as.data.frame(otu_table(mergedPs_percent)), "./mergedOTU_percentages.csv", row.names = TRUE)


pdf(file = "Sampletype_Phylum_Compact.pdf")
plot_bar(mergedPs_percent, x = "Sampletype", fill = "Phylum") +
  facet_wrap(~Sampletype, scales = "free_x") +
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) + 
  guides(fill = guide_legend(ncol = 1))    
dev.off()


pdf(file = "Sampletype_Genus_Compact.pdf")
plot_bar(mergedPs_percent, x = "Sampletype", fill = "Genus") +
  facet_wrap(~Sampletype, scales = "free_x") +
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  guides(fill = guide_legend(ncol = 2))
dev.off()


pdf(file = "Sampletype_Family_Compact.pdf")
plot_bar(mergedPs_percent, x = "Sampletype", fill = "Family") +
  facet_wrap(~Sampletype, scales = "free_x") +
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  guides(fill = guide_legend(ncol = 2))
dev.off()
