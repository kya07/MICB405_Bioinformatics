#!/usr/bin/env RScript 
# MICB 405 RNA-Seq Data Visualization II-GO Terms RScript
# November 13, 2025

# Goals: 
# Generating GO assignments
# Visualizing GO term enrichment with ggplot

# Portions of this code were developed with reference to
# MICB 405 course teaching materials and AI-assisted tools.

#### Load Packages ####
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap)) 
suppressPackageStartupMessages(library(topGO))

# GO term enrichment analysis 
results <- read.csv("results_filtered.csv")
sig_res <- results |> 
  filter(padj < 0.05) #consistent with Data visualization cut off

up_regulated <- sig_res |> 
  filter(log2FoldChange > 2)  |>  #consistent with Data visualization cut off
  pull(X)

down_regulated <- sig_res |> 
  filter(log2FoldChange < -2)  |> 
  pull(X)

#Save data into .txt file to be used in gProfiler
# write(paste(up_regulated, collapse = " "), file = "all_upregulated_genes.txt")
# write(paste(down_regulated, collapse = " "), file = "all_downregulated_genes.txt")


#### Visualizing gProfiler Analyses ####
##### Downregulated gene pathway enrichment ##### 
# 1) Read Data 
down_GO <- read_csv("downreg_check_github.csv", show_col_types = FALSE)

# 2) Create GeneRatio and a TopGO-like "weight01" using -log10(adjusted_p_value)
  # intersection_size / term_size (gProfiler) same as TopGO's Significant / Annotated
  # intersection_size is number of differential genes fall into this GO term
  # term_size is number genes exist in that GO term in the whole database
down_GO_filtered <- down_GO |> 
  mutate(
    GeneRatio = intersection_size / term_size,        # enrichment ratio
    weight01  = negative_log10_of_adjusted_p_value    # bigger = more significant
  )

# 3) Remove non-informative / motif terms 
down_GO_filtered <- down_GO_filtered |> 
  filter(!grepl("motif|match class|Factor:", term_name, ignore.case = TRUE))

# 4) Select the top 20 most significant terms
gprof_downreg_top20 <- down_GO_filtered  |> 
  arrange(desc(weight01)) |> #rank most significant to least significant
  slice_head(n = 20)

# 5) Set term order for plotting (least to most significant)
downreg_order_term <- gprof_downreg_top20 |> 
  arrange(weight01)|> 
  pull(term_name) # pull() extracts a column as a vector

# 6) Plot to visualize gProfiler analyses
downreg_bubble <- gprof_downreg_top20 |> 
  ggplot(aes(x= term_name, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = intersection_size)) + 
  coord_flip() +
  scale_x_discrete(limits = downreg_order_term) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() +
  labs(title = "Top 20 Enriched Pathways for Downregulated Genes (gProfiler)",
       x = "GO Term Description", 
       y = "Enrichment Ratio", 
       color = "-log(p-adj)", 
       size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), 
        panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,0.15), 
                     breaks = seq(0, 2, 0.05), # number to be printed on axis
                     expand = c(0, 0)) # Extra blank space added

#Save
# ggsave("downreg_top20_pathway_enrichment_gProfiler.png", plot = downreg_bubble, width = 12, height = 8)


##### Upregulated gene pathway enrichment ##### 
# Repeat above step 1) to 6) for unregulated genes
up_GO <- read_csv("upreg_check_github.csv", show_col_types = FALSE)

# 2) Create GeneRatio and a TopGO-like "weight01" using -log10(adjusted_p_value)
# intersection_size / term_size (gProfiler) same as TopGO's Significant / Annotated
# intersection_size is number of differential genes fall into this GO term
# term_size is number genes exist in that GO term in the whole database
up_GO_filtered <- up_GO |> 
  mutate(
    GeneRatio = intersection_size / term_size,        # enrichment ratio
    weight01  = negative_log10_of_adjusted_p_value    # bigger = more significant
  )

# 3) Remove non-informative / motif terms 
up_GO_filtered <- up_GO_filtered |>
  filter(!grepl("motif|match class|Factor:", term_name, ignore.case = TRUE))

# 4) Select the top 20 most significant terms
gprof_upreg_top20 <- up_GO_filtered  |> 
  arrange(desc(weight01)) |> #rank most significant to least significant
  slice_head(n = 20)

# 5) Set term order for plotting (least to most significant)
upreg_order_term <- gprof_upreg_top20 |> 
  arrange(weight01)|> 
  pull(term_name) # pull() extracts a column as a vector

# 6) Plot to visualize gProfiler analyses
upreg_bubble <- gprof_upreg_top20 |> 
  ggplot(aes(x= term_name, y = GeneRatio, color = weight01)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = intersection_size)) + 
  coord_flip() +
  scale_x_discrete(limits = upreg_order_term) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_light() +
  labs(title = "Enriched Pathways for Upregulated Genes (gProfiler)",
       x = "GO Term Description", 
       y = "Enrichment Ratio", 
       color = "-log(p-adj)", 
       size = "Number of Significant Genes") + 
  theme(panel.border = element_rect(color = "black"), 
        panel.grid = element_line(colour = "grey96")) +
  scale_y_continuous(limits = c(0,0.2), 
                     breaks = seq(0, 2, 0.05), # number to be printed on axis
                     expand = c(0, 0)) # Extra blank space added

#Save
# ggsave("upreg_top20_pathway_enrichment_gProfiler.png", plot = upreg_bubble, width = 12, height = 8)
