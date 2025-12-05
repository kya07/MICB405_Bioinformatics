#!/usr/bin/env RScript 
# MICB 405 DESeq2 Analysis, VolcanoseR, HeatMaps Visualization RScript
# October 30, 2025

# Portions of this code were developed with reference to
# MICB 405 course teaching materials and AI-assisted tools.

##### Loading required packages ##### 
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(writexl))
suppressPackageStartupMessages(library(ggplot2))

##### Loading data #####

donor1_LPS <- read.table("alignIntronMax_1000000/D1_LPS/ReadsPerGene.out.tab",
                             sep = "\t",
                             col.names = c("gene_id", "total",
                                           "antisense", "regular"),
                             skip = 4)

donor1_unstim <- read.table("alignIntronMax_1000000/D1_unstim//ReadsPerGene.out.tab",
                              sep = "\t",
                              col.names = c("gene_id", "total",
                                            "antisense", "regular"),
                              skip = 4)

donor2_LPS <- read.table("alignIntronMax_1000000/D2_LPS/ReadsPerGene.out.tab",
                             sep = "\t",
                             col.names = c("gene_id", "total",
                                           "antisense", "regular"),
                             skip = 4)

donor2_unstim <- read.table("alignIntronMax_1000000/D2_unstim//ReadsPerGene.out.tab",
                              sep = "\t",
                              col.names = c("gene_id", "total",
                                            "antisense", "regular"),
                              skip = 4)

donor3_LPS <- read.table("alignIntronMax_1000000/D3_LPS/ReadsPerGene.out.tab",
                         sep = "\t",
                         col.names = c("gene_id", "total",
                                       "antisense", "regular"),
                         skip = 4)

donor3_unstim <- read.table("alignIntronMax_1000000/D3_unstim//ReadsPerGene.out.tab",
                            sep = "\t",
                            col.names = c("gene_id", "total",
                                          "antisense", "regular"),
                            skip = 4)

donor4_LPS <- read.table("alignIntronMax_1000000/D4_LPS/ReadsPerGene.out.tab",
                         sep = "\t",
                         col.names = c("gene_id", "total",
                                       "antisense", "regular"),
                         skip = 4)

donor4_unstim <- read.table("alignIntronMax_1000000/D4_unstim//ReadsPerGene.out.tab",
                            sep = "\t",
                            col.names = c("gene_id", "total",
                                          "antisense", "regular"),
                            skip = 4)

##### DESeq2 #####
# DESeq2 is a widely used tool to assess truly differentially expressed genes.
# Let's prepare the read counts for use in DESeq2.

# Build a new data frame with gene names and read counts.

# Stranded data is used as NEBNext Ultra II Directional RNA Library Prep Kit (derives directionality) used in study.
# pick the appropriate column (regular column chosen as have less 0 read counts reads)

macro_readcounts <- data.frame(row.names = donor1_unstim$gene_id,
                               donor1_unstim = donor1_unstim$regular,
                               donor1_LPS = donor1_LPS$regular,
                               donor2_unstim = donor2_unstim$regular,
                               donor2_LPS = donor2_LPS$regular,
                               donor3_unstim = donor3_unstim$regular,
                               donor3_LPS = donor3_LPS$regular,
                               donor4_unstim = donor4_unstim$regular,
                               donor4_LPS = donor4_LPS$regular)

# DESeq2 also requires the read counts to be in matrices, not data frames.
# Matrices also contain data values in a table format, but the entire table's
# values must all be of the same data type.

# We can convert with as.matrix().
macro_matrix <- as.matrix(macro_readcounts)

# DESeq2 also requires a table that provides information about the columns.
# These must be in the order in which they appear on the count matrix. Below enable checking of names.
colnames(macro_readcounts) 

# Confusingly, column metadata is presented as a data frame. The names of each
# sample should be the row name, and be identical and in order of appearance in
# the read counts matrix.
columns_data <- data.frame(conditions = c("control_unstimulated",
                                         "LPS_stimulated",
                                         "control_unstimulated",
                                         "LPS_stimulated",
                                         "control_unstimulated",
                                         "LPS_stimulated",
                                         "control_unstimulated",
                                         "LPS_stimulated"),
                           row.names = c("donor1_unstim",
                                         "donor1_LPS",
                                         "donor2_unstim",
                                         "donor2_LPS",
                                         "donor3_unstim",
                                         "donor3_LPS",
                                         "donor4_unstim",
                                         "donor4_LPS"))

# Create our DESeq2 object.
dds_matrix <- DESeqDataSetFromMatrix(countData = macro_matrix, # Matrix 
                                     colData = columns_data, # Metadata
                                     design = ~conditions)

# Set control condition using the relevel function.
dds_matrix$conditions <- relevel(dds_matrix$conditions, ref = "control_unstimulated")

# Run DESeq2
dds <- DESeq(dds_matrix)

# Saving DESeq2 object into an RDS file for easier access next time
saveRDS(dds, file = "dds.rds")
# With this dds object, we can perform a wide variety of interesting analyses.


####### Analyses ####### 

###### Sanity Check 1: Make a PCA plot with the data######## 
# Perform log transformation on our count data
rld <- rlog(dds)

# Generate a PCA plot with DESeq2's plotPCA function
PCA <- plotPCA(rld, intgroup = "conditions") 
PCA

# ggsave("PCA_plot.png", plot = PCA, width = 8, height = 4)


########  Sanity Check 2: Make a distance matrix with the data ####### 

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate a heatmap using the pheatmap package
heatmap <- pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)
# ggsave("heatmap.png", plot = heatmap, width = 10, height = 6)


###### Further analyses of differentially expressed genes ###### 

# names of the results that DESeq2 calculated
resultsNames(dds) #conditions_LPS_stimulated_vs_control_unstimulated

# Now we will extract the results for our comparison between the LPS-stimulated
# and the Unstimulated conditions.
res <- results(dds, name = "conditions_LPS_stimulated_vs_control_unstimulated") |> as.data.frame()
# we save it as a data frame so we can apply tidyverse dplyr functions better!
head(res)

# Getting rid of NA values
res_no_NA <- res |> drop_na()
head(res_no_NA)
dim(res_no_NA) # Look at how many rows you filtered out!
nrow(res_no_NA) # How many DE genes in total (after removing NAs)?

# Filtering for adjusted p-values <=0.05
res_filtered <- res_no_NA |> filter(padj <= 0.05) #update: nov.27 changed to <0.05
head(res_filtered)
dim(res_filtered) # Look at how many rows you filtered out!
  # Y rows were filtered out
nrow(res_filtered) # How many are significant by adjusted p-value (padj <= 0.05)?

# Pull genes with more than 2x higher/lower expression
res_filtered_final <- res_filtered |>
  filter(log2FoldChange <= -1 | log2FoldChange >= 1) |> 
  rownames_to_column("gene_id") # Convert the rownames into a column so they can be saved in your CSV file
# The '|' stands for OR here!
nrow(res_filtered_final) # How many have significant fold change (padj <= 0.05 AND |log2FC| >= 1)?

##### Summary Stats ####
# How many DE genes in total (after removing NAs)?
total_DE <- nrow(res_no_NA)

# How many are significant by adjusted p-value (padj <= 0.05)?
sig_padj_DE <- nrow(res_filtered)

# How many have significant fold change (padj <= 0.05 AND |log2FC| >= 1)?
sig_FC_DE <- nrow(res_filtered_final)

# Among the significant FC genes, how many are upregulated vs downregulated?
upregulated_DE <- res_filtered_final |>
  filter(log2FoldChange > 0) |>
  nrow()

downregulated_DE <- res_filtered_final |>
  filter(log2FoldChange < 0) |>
  nrow()

# Print a small summary
cat("Total DE genes (no NA):", total_DE, "\n")
cat("Significant by padj (<= 0.05):", sig_padj_DE, "\n")
cat("Significant FC (padj <= 0.05 & |log2FC| >= 1):", sig_FC_DE, "\n")
cat("  Upregulated:", upregulated_DE, "\n")
cat("  Downregulated:", downregulated_DE, "\n")


#### Heatmaps visualization #### 
# z-scores are the amount by which a gene deviates (Standard deviations) 
# in a specific sample from the gene’s average across all samples. 
# i.e., how many standard deviations is my value away from the mean of a row (gene).

###### Heatmap 1: Top 20 genes with highest row variation ###### 
# 1) Prepare data with vst() (similar to rlog()) 
vsd <- assay(vst(dds))

# 2) Calculate z-score values from normalized count data
Z <- t(scale(t(vsd)))

# 3) Extract the row indices (akak the location in the matrix) of the rows with the highest variation
topVarGenes <- rowVars(vsd) |> 
  order(decreasing = TRUE) |> 
  head(20) 

# 4) Subset the Z-score matrix for only the row we're interested in (stored in topVarGenes)
Z_topGenes <- Z[topVarGenes,] 

# Plot using pheatmap()
pheatmap(Z_topGenes,
         main = "Top 20 Genes")


###### HeatMaps 2-3: Top 20 Up regulated & Downregulated Genes###### 
### ---- Extract DESeq2 results ----
res <- results(dds)

### ---- Identify Top 20 Upregulated Genes ----
res_df <- as.data.frame(res) %>% filter(!is.na(log2FoldChange))

top20_up <- res_df %>%
  arrange(desc(log2FoldChange)) %>%  
  head(20) %>%
  rownames()

### ---- Subset Z-score Matrix ----
Z_upGenes <- Z[top20_up,]

### ---- Plot Heatmap ----
pheatmap(
  Z_upGenes,
  main = "Top 20 Upregulated Genes (Z-score Heatmap)",
  cluster_rows = FALSE,
  cluster_cols = TRUE
)

#Downregulated
top20_down <- res_df %>%
  arrange(log2FoldChange) %>%  # most negative log2FC first
  head(20) %>%
  rownames()

### ---- Subset Z-score Matrix ----
Z_downGenes <- Z[top20_down,]

### ---- Plot Heatmap ----
pheatmap(
  Z_downGenes,
  main = "Top 20 Downregulated Genes (Z-score Heatmap)",
  cluster_rows = FALSE,
  cluster_cols = TRUE
)


##### VolcaNoseR visualization ######
# Prepare DESeq2 results table for VolcaNoseR visualization.
# Convert gene IDs from row names into a column and add a -log10(padj) column
# so that VolcaNoseR can plot log2FoldChange (x-axis) vs -log10(padj) (y-axis).

# Instead of having gene ID in the "rowname" column, convert them as the first index column.
res_filtered_final_1 <- res_filtered_final |> 
  rownames_to_column()
# Add a new column to contain the negative log-10 of each padj value.
result_minus_log10 <- res_filtered_final_1 |> 
  mutate(minus_log10_padj = -log10(res_filtered_final_1$padj))

# Write CSV
write.csv(res_filtered_final, file = "results_filtered.csv")
write.csv(result_minus_log10, file = "DESeq2_results_filtered_volcano_input.csv")
