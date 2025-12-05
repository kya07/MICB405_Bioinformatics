library(gprofiler2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(DESeq2)

###
# This script contains the code to conduct gProfiler pathway analysis in R for
# group 7 MICB 405 final project. We examine pathways derived from our up 
# regulated gene set as well as our down regulated gene set and generate 
# plots to visualize the results. plot_deg_per_pathway generate the plots used
# in Figures 4, 5, and 6 of the final report. 
###

# Portions of this code were developed with reference to
# MICB 405 course teaching materials and AI-assisted tools.

# Generates a plot of the top 20 enriched pathways from a gProfiler result with
# intersectin size denoted by size of point and ontology denoted by color.
# Args:
#   gost_res: gProfiler result object returned by gost()
#   output: File path to save the plot
# Returns:
#   Nothing, saves png image of plot
plot_pathways_by_intersection <- function(gost_res, output) {
  res_filtered <- gost_res$result %>% filter(significant)
  
  top_terms <- res_filtered %>%
    group_by(source) %>%
    slice_max(order_by = intersection_size, n = 20) %>%
    ungroup()
  
  
  ggplot(top_terms, aes(x = term_name, y = -log10(p_value), 
                        size = intersection_size, color = source)) +
    geom_point() +
    theme_bw() +
    labs(color = "Ontology") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  ggsave(output)
}


# Generates a heatmap of variance-stabilized expression z score values for genes 
# in a given pathway.
# Args:
#   gost_res: g:Profiler result object returned by gost()
#   output: File path to save the heatmap
#   pathway: Name of the pathway to plot (must match term_name in gost_res)
#   dds: filtered DESeq2 result object
# Returns:
#   Nothing, saves a PNG heatmap to output
plot_deg_per_pathway <- function(gost_res, output, pathway, dds) {
  path <- gost_res$result[
    gost_res$result$term_name == pathway,
  ]
  
  genes_in_pathway <- unlist(strsplit(path$intersection, ","))
  
  vsd <- assay(vst(dds))
  Z <- vsd[rownames(vsd) %in% genes_in_pathway, ] |>
    t() |>
    scale() |>
    t()
  
  Z <- na.omit(Z)
  
  name <- paste0("Expression Levels of Genes Detected in ", pathway)
  pheatmap(
    Z,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_rownames = FALSE,
    fontsize = 9,
    color = colorRampPalette(c("blue", "white", "red"))(200),
    filename=output,
    main=name,
    width=6,
    height=8,
    fontsize_row=12,
    fontsize_col = 12
  )
}

# Load DESeq2 results and apply filtering
results <- read.csv("results_filtered.csv")
sig_res <- results %>%
  filter(padj < 0.05)

up_regulated <- sig_res %>%
  filter(log2FoldChange > 2) %>%
  pull(gene_id)

down_regulated <- sig_res %>%
  filter(log2FoldChange < -2) %>%
  pull(gene_id)

# Pathway analysis for upregulated genes
up_gost_res <- gost(
  query = up_regulated,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "KEGG"),
  correction_method = "fdr",
  evcodes = TRUE
)

# Pathway analysis for downregulated genes
down_gost_res <- gost(
  query = down_regulated,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "KEGG"),
  correction_method = "fdr",
  evcodes = TRUE
)

# Plot top upregulated and downregulated pathways
plot_pathways_by_intersection(up_gost_res, "pathway_analysis/top_upregulated_pathways_intersection_size.png")
plot_pathways_by_intersection(down_gost_res, "pathway_analysis/top_downregulated_pathways_intersection_size.png")


# Read DESeq2_R expression data
dds <- readRDS("dds.rds")
select_upregualted_pathways <- c("chemokine activity", "inflammatory response", "cytokine activity", "TNF signaling pathway", 
                                 "Rheumatoid arthritis", "T-helper 1 cell differentiation", "T-helper 17 cell differentiation",
                                 "T cell differentiation", "T cell differentiation involved in immune response", "T-helper cell differentiation",
                                 "CD4-positive, alpha-beta T cell differentiation involved in immune response", "antimicrobial humoral response",
                                 "regulation of metabolic process")

# Plot gene expression heatmap for upregulated pathways of interest
for(pathway in select_upregualted_pathways) {
  output <- paste0("pathway_analysis/heatmap_dge_up_", gsub(" ", "_", pathway), ".png")
  print(pathway)
  plot_deg_per_pathway(up_gost_res, output, pathway, dds)
}

select_downregulated_pathways <- c("mitotic cell cycle", "tubulin binding", "carbohydrate derivative binding", "lipid binding", "purine nucleotide binding")

# Plot gene expression heatmap for downregulated pathways of interest
for(pathway in select_downregulated_pathways) {
  output <- paste0("pathway_analysis/heatmap_dge_down_", gsub(" ", "_", pathway), ".png")
  plot_deg_per_pathway(down_gost_res, output, pathway, dds)
}
