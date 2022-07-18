# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-07-18
# Updated... 2022-07-18

# Description:
# Run an example GWAS on maize test data, isolate peaks, find
# closest genes within gff file
# Manual: https://maize-genetics.github.io/rTASSEL/
# ---------------------------------------------------------------

# Setting memory
numThreads <- 3
p_value <- 0.01
options(java.parameters = c("-Xmx10g"))

# Load package
library(rTASSEL)
library(dplyr)
library(ggplot2)

# Start logging
rTASSEL::startLogger(fullPath = "/Users/mbb262-admin/Desktop", fileName = "rtassel_logger.txt")


# ------------------------------------------------
# Gather phenotypes, genotypes
# ------------------------------------------------

# Load in genotype table
vcf <-  rTASSEL::readGenotypeTableFromPath(path = "~/Desktop/mdp_genotype.vcf", keepDepth = FALSE)

# Load in phenotype file and format
phenotype <- data.table::fread("~/Desktop/TutorialData/mdp_phenotype.txt", nThread = numThreads)
phenotype <- phenotype[-c(1:3),] # delete formatting columns 
colnames(phenotype) <- c("Taxa", "location", "EarHT", "dpoll", "EarDia", "Q1", "Q2", "Q3") # change column names

# turn columns in phenotype file back into numeric if necessary
phenotype$EarHT <- as.numeric(as.character(phenotype$EarHT))
phenotype$dpoll <- as.numeric(as.character(phenotype$dpoll))
phenotype$EarDia <- as.numeric(as.character(phenotype$EarDia))
phenotype$Q1 <- as.numeric(as.character(phenotype$Q1))
phenotype$Q2 <- as.numeric(as.character(phenotype$Q2))
phenotype$Q3 <- as.numeric(as.character(phenotype$Q3))

# Tasselize merged phenotype + PCs 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = phenotype,
  taxaID = "Taxa",
  attributeTypes = c("factor", "data", "data", "data", "covariate", "covariate", "covariate"))

# Join genotypes with (phenotype + PCs)
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDF,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.1,
  siteMaxAlleleFreq = 1.0)
print(tasGenoPhenoFilt)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = dpoll ~ Q1 + Q2 + Q3,
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = p_value,
  maxThreads = numThreads,
  outputFile = "~/Desktop/test_finding_genes.txt")


# ------------------------------------------------------------------------------------
# Work with analyzed results
# ------------------------------------------------------------------------------------

# Load in association results
gwas_result <- data.table::fread("~/Desktop/test_finding_genes.txt")

# Subset to just one trait
gwas_result <- gwas_result %>% filter(Trait == "EarHT")

# Format x axis to plot by chromosome
data_cum <- gwas_result %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)
gwas_result <- gwas_result %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Do a multiple test correction for p-values, FDR method
gwas_result$p_adj <- p.adjust(gwas_result$p, method = "fdr")

# Plot out
ggplot(gwas_result, aes(x = bp_cum, y = -log10(p_adj), color = as.factor(Chr))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "Ear Height Tassel test data") + 
  theme_minimal()


# ------------------------------------------------------------------------------------
# Find top peaks and closest genes
# ------------------------------------------------------------------------------------

# subset significant sites by chromosomes, assumes one top SNP per chromosome
peaks <- gwas_result %>% 
  filter(Chr %in% c(1, 7, 9, 10)) %>% # filter to chroms with significant peaks, if all chromosomes, delete this line
  group_by(Chr) %>% # group results by chromosome
  slice_min(p_adj) # select single SNP with the smallest adjusted p-value

# IF there are multiple peaks on the same chromosome do something like this:
# Example works if, for example, there are two peaks on chrom 1. Run this code snippet twice
# changing the filtering values by position of the peak
peaks2 <- gwas_result %>% 
  filter(Chr == 1 & Pos >= before_start_position_of_peak & Pos <= after_end_position_of_peak) %>% 
  group_by(Chr) %>% # group results by chromosome
  slice_min(p_adj) # select single SNP with the smallest adjusted p-value


## NOTE: At this step you can take these chromosome and position combinations and look through your 
# genome's genome browser (if you have one) and look around for genes in the vicinity. If your species
# does not have a genome browser, or you want a repeatable code way of doing this, use what's below


# Load in a gff file using ape
gff <- ape::read.gff("~/git_projects/mann_gwas/data/Zea_mays.B73_RefGen_v4.49.gff3.gz", 
                     na.strings = c(".", "?"), 
                     GFF3 = TRUE)

# Remove extra columns, only select chromosomes and no scaffolds
gff <- gff %>% 
  filter(type == "gene", source == "gramene", seqid %in% c(1,2,3,4,5,6,7,8,9,10)) %>% 
  select("seqid", "start", "end", "strand", "attributes")
gff$chrom <- as.numeric(as.character(gff$chrom)) # turn chromosomes into numeric

# Parse out gene names from metadata column
gff$attributes <- gsub("ID=gene:", "", gff$attributes)
gff$attributes <- gsub(";.*", "", gff$attributes)

# filter gff files to regions overlapping with GWAS hits
# The distance to nearest gene is able to be modified
# Function takes into account positive and negative strands of genes
findGene <- function(gwas_peaks, distance, gff_file){
  genes <- c()
  for(i in 1:nrow(gwas_peaks)){
    # Subset gwas results to just necessary things
    gwas_peaks_sub <- gwas_peaks %>% select("Trait", "Marker", "Chr", "Pos", "p")
    
    # For genes on negative strand
    gff_neg <- gff_file %>% filter(strand == "-")
    peak_neg_genes <- gff_neg %>% 
      filter(seqid == gwas_peaks$Chr[i] & (start >= gwas_peaks$Pos[i] + distance) & (end <= gwas_peaks$Pos[i] - distance))
    
    # If genes nearby are returned, add peak information
    if(nrow(peak_neg_genes) > 0){
      peak_neg_genes <- cbind(gwas_peaks_sub[i,], peak_neg_genes)
    }
    
    # For genes on positive strand
    gff_pos <- gff_file %>% filter(strand == "+")
    peak_pos_genes <- gff_pos %>% 
      filter(seqid == gwas_peaks$Chr[i] & (start >= gwas_peaks$Pos[i] - distance) & (end <= gwas_peaks$Pos[i] + distance))

    # If genes nearby are returned, add peak information
    if(nrow(peak_pos_genes) > 0){
      peak_pos_genes <- cbind(gwas_peaks_sub[i,], peak_pos_genes)
    } 
    
    # combine results
    pos_neg_genes <- rbind(peak_neg_genes, peak_pos_genes)
    genes <- rbind(genes, pos_neg_genes)
  }
  return(genes)
}

# Use the function, look for genes within 50 Kb of the significant SNP
findGene(gwas_peaks = peaks, distance = 50000, gff_file = gff)



