library(gprofiler2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(DESeq2)

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

results <- read.csv("DESeq2_R/results_filtered.csv")
sig_res <- results %>%
  filter(padj < 0.05)

up_regulated <- sig_res %>%
  filter(log2FoldChange > 2) %>%
  pull(gene_id)

down_regulated <- sig_res %>%
  filter(log2FoldChange < -2) %>%
  pull(gene_id)

up_gost_res <- gost(
  query = up_regulated,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "KEGG"),
  correction_method = "fdr",
  evcodes = TRUE
)

down_gost_res <- gost(
  query = down_regulated,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "KEGG"),
  correction_method = "fdr",
  evcodes = TRUE
)

plot_pathways_by_intersection(up_gost_res, "pathway_analysis/top_upregulated_pathways_intersection_size.png")
plot_pathways_by_intersection(down_gost_res, "pathway_analysis/top_downregulated_pathways_intersection_size.png")


dds <- readRDS("DESeq2_R/dds.rds")
select_upregualted_pathways <- c("chemokine activity", "inflammatory response", "cytokine activity", "TNF signaling pathway", 
                                 "Rheumatoid arthritis", "T-helper 1 cell differentiation", "T-helper 17 cell differentiation",
                                 "T cell differentiation", "T cell differentiation involved in immune response", "T-helper cell differentiation",
                                 "CD4-positive, alpha-beta T cell differentiation involved in immune response", "antimicrobial humoral response",
                                 "regulation of metabolic process")

for(pathway in select_upregualted_pathways) {
  output <- paste0("pathway_analysis/heatmap_dge_up_", gsub(" ", "_", pathway), ".png")
  print(pathway)
  plot_deg_per_pathway(up_gost_res, output, pathway, dds)
}

select_downregulated_pathways <- c("mitotic cell cycle", "tubulin binding", "carbohydrate derivative binding", "lipid binding", "purine nucleotide binding")

for(pathway in select_downregulated_pathways) {
  output <- paste0("pathway_analysis/heatmap_dge_down_", gsub(" ", "_", pathway), ".png")
  plot_deg_per_pathway(down_gost_res, output, pathway, dds)
}

pathway <- "regulation of macromolecule metabolic process"
output <- paste0("pathway_analysis/heatmap_dge_up_", gsub(" ", "_", pathway), ".png")
plot_deg_per_pathway(up_gost_res, output, pathway, dds)

pathway <- "oxidative phosphorylation"
output <- paste0("pathway_analysis/heatmap_dge_down_", gsub(" ", "_", pathway), ".png")
plot_deg_per_pathway(down_gost_res, output, pathway, dds)
