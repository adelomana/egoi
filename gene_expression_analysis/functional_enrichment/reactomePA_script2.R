#
# -1. packages installation
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# setRepositories(ind=c(1:6))
# 
# BiocManager::install("ReactomePA")
# BiocManager::install('org.Hs.eg.db')
BiocManager::install('convertid')

#
# 0. load libraries
#
library(ReactomePA)
library(convertid)

#
# 1. user-defined variables
#
gene_list_folder = '/Users/adrian/gd15/tmp/hegoi_tempo/functional_enrichment_sleuth'

#
# 2. read data
#
gene_list_file = paste(gene_list_folder, 'hypo_A_up.tsv', sep='/')
gene_list_file
gene_list = read.csv(gene_list_file)
print(gene_list)

#
# 3. enrichment
#
x <- enrichPathway(gene=gene_list, pvalueCutoff=0.05, readable=T)




