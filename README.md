# R-Metagenomics

Metagenomics is a technique which allows researchers to infer the different microorganisms present within a microbiome. 
This is commonly achieved using universal primers to generate amplicons from bulk DNA samples representing hypervariable regions from the 16S ribosomal RNA (rRNA) gene of microorganisms. 
Sequencing is performed on these amplicons and then bioinformatics tools such as QIIME and DADA2 are used to match the sequence of these amplicons to known taxa of microorganisms. 
Phylogenetics tools such as Phyloseq can then be used to compare the taxonomic classification of the microbiome between samples and groups of samples.

More powerful techniques such as shotgun Metagenomics allow identification of novel organisms, as well as providing more information about their genomics. Scripts for this will be added soon. 

---

## Table of Contents
1. [Overview](#overview)
2. [Pipeline Stages](#pipeline-stages)
    - [Data Sources](#data-sources)
    - [DADA2](#dada2)
    - [Phyloseq](#phyloseq)
3. [Results](#results)

---

## Overview

The example data used here is taken from a research paper by Kistler J.O. et al. 2015 (ref1). 
I created an analysis pipeline by adapting the DADA2 Pipeline Tutorial (1.16) from Benjamin Callahan: https://benjjneb.github.io/dada2/tutorial.html (Figure 1A). 
This pipeline was then applied to the sequencing data from Kistler J.O. et al. 2015 using HIV positive Plaque and Saliva samples in an attempt to replicate their findings.

References
1- Kistler JO, Arirachakaran P, Poovorawan Y, Dahlén G, Wade WG. The oral microbiome in human immunodeficiency virus (HIV)-positive individuals. 
J Med Microbiol. 2015 Sep;64(9):1094-1101. doi: 10.1099/jmm.0.000128. Epub 2015 Jul 14. PMID: 26297584.

---

## Pipeline Stages

### Data Sources

Samples representing 16S rRNA gene sequencing (V1-V3 hypervariable regions) Universal primers 27FYM and 519R of HIV positive and negative patient saliva and plaque samples were obtained from the European Nucleotide Archive (ENA) Project: PRJNA285477.

The 16S training dataset consisted of DADA2-formatted training fasta files containing taxonomically assigned sequences from the Silva Project's version 132 release DOI 10.5281/zenodo.1172782.

## DADA2
The following steps were taken to prepare the data:

1. Filter and trim reads
- The plotQualityProfile function of DADA2 was used to inspect the quality of the sequencing reads (Figure 1B).
- This analysis indicates that the quality of reads drops off at around 400 bases.
- Unlike Illumina sequencing, Roche 454 pyrosequencing generates reads of quite variable lengths (Figure 1B - red line).
- Reads were filtered for those with minLen 400 and then trimmed to 400 bases (truncLen = 400) to remove poor quality regions and filtered.
  
2. Learn errors and apply inference algorithm
- The learnErrors method of DADA2 is used to model the nucleotide error rate in the input samples. 
- The estimate error rates for each possible nucleotide transition can be visualised using the plotErrors function (Figure 1C).
- To correct for amplicon sequencing errors, the core sample inference algorithm is applied to the data using the error model from the previous step.
- Code has been modified to account for the Roche 454 pyrosequencing chemistry.

3. Build ASV table and remove chimeras
- A sequence variant table (ASV) table is a higher-resolution version of an OTU table.
- An ASV is a matrix file with a row for each sample and a column for each sequence variant.
- Chimeras are artifacts created by the incorrect joining of two or more biological sequences; these are removed from the table by the removeBimeraDenovo function. 

4. Assign taxonomy
- DADA2 uses a  naive Bayesian classifier method to add taxonomic assignment to our samples based on the Silva Project version 132 training dataset.
  
5. Output Files
- read_quality.pdf :    Heat map of the frequency of each quality score at each base position (grey shading). 
                        Mean quality score at each position (green line) with quartiles (orange line).
                        Scaled proportion of reads that extend to at least that position (red line).
- error_rate.pdf :      Plots of estimated versus observed error rates for each nucleotide transition. 
- DADA2_processing_results.csv :    Table of the number of reads that made it through each step in the pipeline.
- taxa.csv :    Table of the taxonomic assignments for each sequence.

## Phyloseq

1. Build phylogenic relationship object 
- Use the DADA2 output to build a Phyloseq object.
  
2. Analyse the diversity of the samples
- Working on the sample type level, calculate the alpha diversity (richness) of the samples.
- Here I have used two different equations to calculate the Shannon Index (Figure 1D) and the Simpson Index (Figure 1E).
- A Wilcoxon Test is used to determine if the mean diversity is significantly different between sample types.
- To confirm the separation of the samples within sample types an ordination method called Bray NMDS was used to plot the data. 

3. Taxonomy
- The Phyloseq object was used to plot the relative abundance of sequences representing the different families of microbes in each sample.
- The top 20 sequences was selected for each sample, this focusses on the most highly represented families. 
- Also generated 
  
4. Output Files
- Richness_Shannon_Simpson_with_stats.pdf :  Plot displaying Shannon and Simpson diversity plots comparing sample type groups 
                                            Including p values for the Wilcoxon statistical tests.
- Ordination.pdf : Bray NMDS plot comparing samples
- Multiple Phylum, Family and Genus graphs.

## Results

The results reflect the conclusions of Kistler J.O. et al. 2015 (ref1) that there is a different microbiome in HIV positive individuals if sampled for plaque and saliva.
Simpson and Shannon diversity indices were higher in plaque versus saliva samples (Figure 1D and E). Kistler J. O. et al. also observed this difference in Simpson index.
Significant differences in abundance of specific species of microbes in plaque and saliva were observed by Kister J. O. et al. which show several Streptococcus species enriched in Saliva and Corynebacterium species enriched in Plaque. This matches analysis here (Figure 1G)


<img src="https://github.com/focyte/R-Metagenomics/blob/main/R-Metagenomics_Figure.png" alt="Results" width="500"/>
