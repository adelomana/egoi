#
# -1. packages installation
#
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

setRepositories(ind=c(1:6))

BiocManager::install("igraph") # igraph seems to be failing for ReactomePA
BiocManager::install("ReactomePA")

#
# 0. load libraries
#
library(ReactomePA)

#
# 1. read data
#
gene_list = read.csv("/home/adrian/projects/hegoi/results/functional_enrichment_sleuth/hypo_A_up.tsv")
head(gene_list)