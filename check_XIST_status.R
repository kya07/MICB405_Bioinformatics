library(tidyverse)
library(DESeq2)
library(ggplot2)

#load completed DESeq2 Object
dds <- readRDS("DESeq2_R/dds.rds")

# get normalized counts from the DESeq2 object
norm_counts <- counts(dds, normalized = TRUE)

# extract XIST expression
xist_expr <- norm_counts["XIST", ]

df <- data.frame(
  sample = colnames(dds),
  sex = colData(dds)$sex,
  xist = as.numeric(xist_expr)
)
