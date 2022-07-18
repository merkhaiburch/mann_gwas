# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-07-18
# Updated... 2022-07-18

# Description:
# Perform GWAS on bacterial titer measured by RT-PCR of CLas 
# within the stomachs of psyllids attacking oranges
#
# Run many different model architectures: 
# 1) titer ~ SNP + management
# 2) titer ~ SNP + PC1 + PC2
# 3) titer ~ SNP + management + PC1 + PC2
# 4) titer BLUE ~ SNP + PC1 + PC2
# 5) titer BLUE ~ SNP + Kinship + PC1 + PC2
# ---------------------------------------------------------------


# Setting memory
options(java.parameters = c("-Xmx120g"))
numThreads <- 20
p_value <- 0.01

# Load package
library(rTASSEL)
library(dplyr)
library(ggplot2)

# Start logging
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)


# ------------------------------------------------
# Gather phenotypes, genotypes
# ------------------------------------------------

# Load in genotype table
vcfpath <- "/workdir/mbb262/Mann_GWA/"
vcf <-  rTASSEL::readGenotypeTableFromPath(path = paste(vcfpath,"DC_BEAGLE_OUTPUT.vcf.gz", sep = ""), keepDepth = FALSE)

# Load in phenotype file
phenotype <- data.table::fread("/workdir/mbb262/Mann_GWA/Pheno-Absolute-log10.csv", nThread = numThreads, drop = c(1))
colnames(phenotype) <- c("Taxa", "pheno")

# Load in pcs
pc_file <- data.table::fread("/workdir/mbb262/Mann_GWA/pc_file.txt")
colnames(pc_file)[1] <- "Taxa"
pc_file <- pc_file[,1:3] # Only keep the first 3 PCs

# Management file
management <- data.table::fread("/workdir/mbb262/Mann_GWA/env.txt", drop = c(1))
colnames(management) <- c("Taxa", "field")


# ----------------------------------------------------------------------------------------------
# Run Fast Association
# Model 1: y ~ SNP + management
# ----------------------------------------------------------------------------------------------

# Intersect phenotype and PCs
pheno_management <- merge(x = phenotype, y = management)
pheno_management$field <- as.numeric(as.factor(pheno_management$field))

# Tasselize merged phenotype + PCs 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_management,
  taxaID = "Taxa",
  attributeTypes = c("data", "covariate"))

# Join genotypes with (phenotype + PCs)
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDF,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.1,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 13)
print(tasGenoPhenoFilt)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = pheno ~ .,
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = p_value,
  maxThreads = numThreads,
  outputFile = "/workdir/mbb262/fast_assoc_results_model_1.txt")

# Load in association resuilts
gwas_result_fa_model_1 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_1.txt")

# Format x axis to plot by chromosome
data_cum <- gwas_result_fa_model_1 %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)
gwas_result <- gwas_result_fa_model_1 %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Plot out
model1 <- ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "Model 1: y ~ SNP + management") + 
  theme_minimal()


# ----------------------------------------------------------------------------------------------
# Run Fast Association
# Model 2: y ~ SNP + PC1 + PC2 
# ----------------------------------------------------------------------------------------------

# Intersect phenotype and PCs
pheno_pcs <- merge(x = phenotype, y = pc_file)

# Tasselize merged phenotype + PCs 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_pcs,
  taxaID = "Taxa",
  attributeTypes = c("data", "covariate","covariate"))

# Join genotypes with (phenotype + PCs)
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDF,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.1,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 13)
print(tasGenoPhenoFilt)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = pheno ~ .,
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = p_value,
  maxThreads = numThreads,
  outputFile = "/workdir/mbb262/fast_assoc_results.txt")

# Load in association results
gwas_result_fa_model_2 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_2.txt")

# Format x axis to plot by chromosome
data_cum <- gwas_result_fa_model_2 %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)
gwas_result <- gwas_result_fa_model_2 %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Plot out
model2 <- ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "Model 2: y ~ SNP + PC1 + PC2") + 
  theme_minimal()


# ----------------------------------------------------------------------------------------------
# Run Fast Association
# Model 3: y ~ SNP + PC1 + PC2 + management
# ----------------------------------------------------------------------------------------------

# Intersect phenotype and PCs
pheno_pcs <- merge(x = phenotype, y = pc_file)
pheno_pcs_management <- merge(x = pheno_pcs, y = management)
pheno_pcs_management$field <- as.numeric(as.factor(pheno_pcs_management$field))

# Tasselize merged phenotype + PCs + management 
tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_pcs_management,
  taxaID = "Taxa",
  attributeTypes = c("data", "covariate","covariate", "covariate"))

# Join genotypes with (phenotype + PCs)
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDF,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.1,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 13)
print(tasGenoPhenoFilt)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = pheno ~ .,
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = p_value,
  maxThreads = numThreads,
  outputFile = "/workdir/mbb262/fast_assoc_results_model_3.txt")

# Load in association results
gwas_result_fa_model_3 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_3.txt")

# Format x axis to plot by chromosome
data_cum <- gwas_result_fa_model_3 %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)

gwas_result <- gwas_result_fa_model_3 %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Plot out
model3 <- ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "Model 3: y ~ SNP + PC1 + PC2 + management") + 
  theme_minimal()


# ----------------------------------------------------------------------------------------------
# Run Fast Association
# Model 4: BLUE ~ SNP + PC1 + PC2
# ----------------------------------------------------------------------------------------------

# Intersect phenotypes and PCs
pheno_manage <- merge(x = phenotype, y = management)

# Specify data types
pheno_manage <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = pheno_manage,
  taxaID = "Taxa",
  attributeTypes = c("data", "factor"))

# Calculate BLUEs
tasBLUE <- rTASSEL::assocModelFitter(
  tasObj = pheno_manage,
  formula = . ~ .,
  fitMarkers = FALSE,
  kinship = NULL,
  fastAssociation = FALSE
)

# Merge blues with pcs
tasPhenoDF <- merge(x = tasBLUE$BLUE, y = pc_file)

tasPhenoDF <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = tasPhenoDF,
  taxaID = "Taxa",
  attributeTypes = c("data", "covariate", "covariate"))

# Join genotypes with (phenotypes + g PCs)
tasPhenoDF <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDF,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.1,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 13)
print(tasGenoPhenoFilt)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = pheno ~ .,
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = TRUE,
  maxP = p_value,
  maxThreads = numThreads,
  outputFile = "/workdir/mbb262/fast_assoc_results_model_4.txt")

# Load in association resuilts
gwas_result_fa_model_4 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_4.txt")

# Format x axis to plot by chromosome
data_cum <- gwas_result_fa_model_4 %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)
gwas_result <- gwas_result_fa_model_4 %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Plot out
model4 <- ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "Model 4: BLUE ~ SNP + PC1 + PC2") + 
  theme_minimal()


# ----------------------------------------------------------------------------------------------
# Run Mixed linear model
# Model 5: titer BLUE ~ SNP + Kinship + PC1 + PC2
# ----------------------------------------------------------------------------------------------

# Merge blues with pcs (calculated in Model 4 code chunk)
tasPhenoDF_blue <- merge(x = tasBLUE$BLUE, y = pc_file)

# Tasselize merged BLUE'd phenotype + PCs 
tasPhenoDF_blue <- rTASSEL::readPhenotypeFromDataFrame(
  phenotypeDF = tasPhenoDF_blue,
  taxaID = "Taxa",
  attributeTypes = c("data", "factor"))

# Join genotypes with (BLUE'd phenotype + PCs)
tasPhenoDF_blue <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = vcf,
  phenoPathDFOrObj = tasPhenoDF_blue)

# Do a light MAF filter to remove invariant sites
tasGenoPhenoFilt <- rTASSEL::filterGenotypeTableSites(
  tasObj = tasPhenoDF_blue,
  siteRangeFilterType = "none",
  siteMinAlleleFreq = 0.01,
  siteMaxAlleleFreq = 1.0,
  siteMinCount = 13)
print(tasGenoPhenoFilt)

# Calculate kinship matrix
tasKin <- kinshipMatrix(tasObj = tasGenoPhenoFilt)

# Run fast association, write files to disk.
rTASSEL::assocModelFitter(
  tasObj = tasGenoPhenoFilt,
  formula = pheno ~ .,
  fitMarkers = TRUE,
  kinship = tasKin,
  fastAssociation = TRUE,
  maxP = p_value,
  maxThreads = numThreads,
  outputFile = "/workdir/mbb262/fast_assoc_results_model_5.txt")

# Load in association results
gwas_result_fa_model_5 <- data.table::fread("/workdir/mbb262/fast_assoc_results_model_5.txt")

# Format x axis to plot by chromosome
data_cum <- gwas_result_fa_model_5 %>% 
  group_by(Chr) %>% 
  summarise(max_bp = max(Pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(Chr, bp_add)
gwas_result <- gwas_result_fa_model_5 %>% 
  inner_join(data_cum, by = "Chr") %>% 
  mutate(bp_cum = Pos + bp_add)

# Plot everything out
model5 <- ggplot(gwas_result, aes(x = bp_cum, y = -log10(p), color = as.factor(Chr), size = -log10(p))) +
  geom_point(size = 1) +
  geom_hline(yintercept = 5, color = "grey40", linetype = "dashed") +
  labs(x = NULL, y = "-log10(p)", title = "Model 5: BLUE ~ SNP + PC1 + PC2 + K") + 
  theme_minimal()

# Plot all things together
model1/model2/model3/model4/model5+ 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
ggpubr::ggarrange(model1, model2, model3, model4, model5, nrow = 5, ncol = 1,  common.legend = TRUE)


