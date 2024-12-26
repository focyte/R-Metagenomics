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
                     truncLen=400,
                     maxN=0, 
                     maxEE= Inf, 
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
dadaReads <- dada(filtered.reads, err=error.rate, multithread=TRUE)
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

# Assign "When" column
samdf$When <- "Unknown" 
samdf$When[samdf$Sampletype == "P"] <- "Plaque"
samdf$When[samdf$Sampletype == "S"] <- "Saliva"

# Set rownames to match the original sample names
rownames(samdf) <- samples.out

# Construct a phyloseq object directly from the dada2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Identify and remove taxa with no NA values in their taxonomic classification
valid_taxa <- taxa_names(ps)[!apply(is.na(tax_table(ps)), 1, any)]
ps <- prune_taxa(valid_taxa, ps)
ps

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

##           Functions to apply to the phyloseq results during 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to aggregate data by Sampletype and calculate percentages
aggregate_by_Sampletype <- function(ps, tax_level) {
  ps.melt <- psmelt(ps)
  ps.grouped <- ps.melt %>%
    group_by(Sampletype, !!sym(tax_level)) %>%
    summarise(Abundance = sum(Abundance)) %>%
    mutate(Percentage = 100 * Abundance / sum(Abundance))
  return(ps.grouped)
}

# Function to prune to the top 20 taxa
prune_top_taxa <- function(ps, n = 20) {
  top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:n]
  ps.top <- prune_taxa(top_taxa, ps)
  return(ps.top)
}

# Function to save plots in multiple formats
save_plot <- function(plot, filename) {
  ggsave(paste0(filename, ".pdf"), plot, width = 8, height = 6)
  ggsave(paste0(filename, ".png"), plot, width = 8, height = 6)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##                        Plotting of the data

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Generate and save barplots for all taxa and top 20 taxa
for (tax_level in c("Phylum", "Family", "Genus")) {
  
  # ----- All Taxa -----
  # By Subject
  p_subject_all <- plot_bar(ps.percent, x = "Subject", fill = tax_level) +
    facet_wrap(~When, scales = "free_x") +
    labs(title = paste("Relative Abundance by Subject (All Taxa) -", tax_level),
         y = "Percentage (%)",
         x = "Subject")
  save_plot(p_subject_all, paste0("Subject_All_", tax_level))
  
  # By Sampletype
  ps.Sampletype_all <- aggregate_by_Sampletype(ps.percent, tax_level)
  p_Sampletype_all <- ggplot(ps.Sampletype_all, aes(x = Sampletype, y = Percentage, fill = !!sym(tax_level))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Relative Abundance by Sampletype (All Taxa) -", tax_level),
         y = "Percentage (%)",
         x = "Sampletype") +
    theme_minimal()
  save_plot(p_Sampletype_all, paste0("Sampletype_All_", tax_level))
  
  # ----- Top 20 Taxa -----
  ps.top20 <- prune_top_taxa(ps.percent, 20)
  
  # By Subject
  p_subject_top20 <- plot_bar(ps.top20, x = "Subject", fill = tax_level) +
    facet_wrap(~When, scales = "free_x") +
    labs(title = paste("Relative Abundance by Subject (Top 20 Taxa) -", tax_level),
         y = "Percentage (%)",
         x = "Subject")
  save_plot(p_subject_top20, paste0("Subject_Top20_", tax_level))
  
  # By Sampletype
  ps.Sampletype_top20 <- aggregate_by_Sampletype(ps.top20, tax_level)
  p_Sampletype_top20 <- ggplot(ps.Sampletype_top20, aes(x = Sampletype, y = Percentage, fill = !!sym(tax_level))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Relative Abundance by Sampletype (Top 20 Taxa) -", tax_level),
         y = "Percentage (%)",
         x = "Sampletype") +
    theme_minimal()
  save_plot(p_Sampletype_top20, paste0("Sampletype_Top20_", tax_level))
}

# Print confirmation
message("Barplots for all taxa and top 20 taxa saved as both PDF and PNG files.")
