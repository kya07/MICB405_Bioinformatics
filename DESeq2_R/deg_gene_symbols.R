# Sophia's code for extracting gene symbols
library(dplyr)

results <- read.csv("DESeq2_R/results_filtered.csv")
sig_res <- results %>%
  filter(padj < 0.05)

up_regulated <- sig_res %>%
  filter(log2FoldChange > 0) %>%
  pull(X)

down_regulated <- sig_res %>%
  filter(log2FoldChange < 0) %>%
  pull(X)

write(paste(up_regulated, collapse = " "), file = "all_upregulated_genes.txt")
write(paste(down_regulated, collapse = " "), file = "all_downregulated_genes.txt")

M1_markers <- c("IFNG", "TNF", "IL1A", "IL1B", "IL6", "NOS2", "TLR2", "TLR4", "CD80", "CD86")
M2_markers <- c("IL4", "IL10", "IL13", "CSF1R", "MRC1", "PPARG", "ARG1", "CD163", "CLEC10A", "CLEC7A", "PDCD1LG2", "RETNLA")
tlr <- c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "TLR10")

# Normalize case to match your gene symbols
up_upper <- toupper(up_regulated)
down_upper <- toupper(down_regulated)

# Check for overlap
M1_up <- intersect(up_upper, M1_markers)
M1_down <- intersect(down_upper, M1_markers)
M2_up <- intersect(up_upper, M2_markers)
M2_down <- intersect(down_upper, M2_markers)
TLRu <- intersect(up_upper, tlr)
TLRd <- intersect(down_upper, tlr)

# Print results
cat("M1 markers upregulated:", paste(M1_up, collapse = ", "), "\n")
cat("M1 markers downregulated:", paste(M1_down, collapse = ", "), "\n\n")
cat("M2 markers upregulated:", paste(M2_up, collapse = ", "), "\n")
cat("M2 markers downregulated:", paste(M2_down, collapse = ", "), "\n")

sig_res <- sig_res %>%
  mutate(GENE_UPPER = toupper(X))

# Extract stats for any M1 or M2 marker present
M1_hits <- sig_res %>%
  filter(GENE_UPPER %in% M1_markers) %>%
  arrange(desc(log2FoldChange))

tlr_hits <- sig_res %>%
  filter(GENE_UPPER %in% tlr) %>%
  arrange(desc(log2FoldChange))

M2_hits <- sig_res %>%
  filter(GENE_UPPER %in% M2_markers) %>%
  arrange(desc(log2FoldChange))
