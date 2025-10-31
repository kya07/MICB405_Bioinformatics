# Lecture 6.2 R script
# Running DESeq2
# Completed version
# MICB 405 101 (2025W1)
# October 23, 2025


# Loading required packages -----------------------------------------------

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)


# Loading data ------------------------------------------------------------

chlamy_1h_rep1 <- read.table("chlamy_1h_rep1_counts.tsv",
                             sep = "\t",
                             col.names = c("gene_id", "total",
                                           "antisense", "regular"),
                             skip = 4)

chlamy_12h_rep1 <- read.table("chlamy_12h_rep1_counts.tsv",
                              sep = "\t",
                              col.names = c("gene_id", "total",
                                            "antisense", "regular"),
                              skip = 4)

chlamy_1h_rep2 <- read.table("chlamy_1h_rep2_counts.tsv",
                             sep = "\t",
                             col.names = c("gene_id", "total",
                                           "antisense", "regular"),
                             skip = 4)

chlamy_12h_rep2 <- read.table("chlamy_12h_rep2_counts.tsv",
                              sep = "\t",
                              col.names = c("gene_id", "total",
                                            "antisense", "regular"),
                              skip = 4)

# DESeq2 is a widely used tool to assess truly differentially expressed genes.
# Let's prepare the read counts for use in DESeq2.

# Build a new data frame with gene names and read counts.

# If you are using all stranded data, pick the appropriate column.
# If using non-stranded data, use total column.

chlamy_readcounts <- data.frame(row.names = chlamy_1h_rep1$gene_id,
                                chlamy_1h_rep1 = chlamy_1h_rep1$regular,
                                chlamy_1h_rep2 = chlamy_1h_rep2$regular,
                                chlamy_12h_rep1 = chlamy_12h_rep1$regular,
                                chlamy_12h_rep2 = chlamy_12h_rep2$regular)

# DESeq2 also requires the read counts to be in matrices, not data frames.
# Matrices also contain data values in a table format, but the entire table's
# values must all be of the same data type.

# We can convert with as.matrix().
chlamy_matrix <- as.matrix(chlamy_readcounts)

# DESeq2 also requires a table that provides information about the columns.
# These must be in the order in which they appear on the count matrix.
colnames(chlamy_readcounts)

# Confusingly, column metadata is presented as a data frame. The names of each
# sample should be the row name, and be identical and in order of appearance in
# the read counts matrix.
columns_data <- data.frame(timepoint = c("1_hour",
                                         "1_hour",
                                         "12_hour",
                                         "12_hour"),
                           row.names = c("chlamy_1h_rep1",
                                         "chlamy_1h_rep2",
                                         "chlamy_12h_rep1",
                                         "chlamy_12h_rep2"))

# Create our DESeq2 object.
dds_matrix <- DESeqDataSetFromMatrix(countData = chlamy_matrix, # Matrix 
                                     colData = columns_data, # Metadata
                                     design = ~timepoint)

# Set control condition using the relevel function.
dds_matrix$timepoint <- relevel(dds_matrix$timepoint, ref = "1_hour")

# Run DESeq2
dds <- DESeq(dds_matrix)

# With this dds object, we can perform a wide variety of interesting analyses.


# Analyses covered in Tutorial 7 ------------------------------------------

# Sanity Check 1: Make a PCA plot with the data ---------------------------

# Perform log transformation on our count data
rld <- rlog(dds)

# Generate a PCA plot with DESeq2's plotPCA function
plotPCA(rld, intgroup = "condition") 


# Sanity Check 2: Make a distance matrix with the data --------------------

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate a heatmap using the pheatmap package
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)


# Further analyses of differentially expressed genes ----------------------

# names of the results that DESeq2 calculated
resultsNames(dds)

# Now we will extract the results for our comparison between the 12h timepoint
# and the 1h timepoint.
res <- results(dds, name = "timepoint_12_hour_vs_1_hour") |> as.data.frame()
# we save it as a data frame so we can apply tidyverse dplyr functions better!
head(res)

# Getting rid of NA values
res_no_NA <- res |> drop_na()
head(res_no_NA)
dim(res_no_NA) # Look at how many rows you filtered out!

# Filtering for adjusted p-values <=0.05
res_filtered <- res_no_NA |> filter(padj <= 0.05)
head(res_filtered)
dim(res_filtered) # Look at how many rows you filtered out!

# Pull genes with more than 2x higher/lower expression
res_filtered_final <- res_filtered |>
  filter(log2FoldChange <= -1 | log2FoldChange >= 1)
# The '|' stands for OR here!
head(res_filtered_final)
dim(res_filtered_final)

# Top 10 genes positively differentially expressed.
top10_genes <- res_filtered_final |>
  arrange(desc(log2FoldChange)) |>
  # We use the desc() function to organize the column in descending order.
  head(n = 10)
top10_genes

# How about the 10 genes negatively differentially expressed?
bot10_genes <- res_filtered_final |>
  arrange(log2FoldChange) |>
  # Since we don't use desc(), the column is organized in ascending order.
  head(n = 10)
bot10_genes

# Write CSV
write.csv(res_filtered_final, file = "results_filtered.csv")