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
