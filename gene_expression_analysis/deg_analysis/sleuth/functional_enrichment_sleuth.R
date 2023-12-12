#! BiocManager::install("enrichR")
#! BiocManager::install("enrichplot")
#! BiocManager::install('ReactomePA')
#! BiocManager::install('org.Hs.eg.db')
#! BiocManager::install('clusterProfiler')
#! BiocManager::install('OrgDb')

#
# 0. load libraries
#
library(crayon) # so messages can be in blue font
library(clusterProfiler)
library(enrichplot)

#
# 1. user-defined variables
#
gene_sets_dir = '/home/adrian/projects/hegoi/results/functional_enrichment_sleuth'
output_dir = '/home/adrian/projects/hegoi/results/functional_enrichment_sleuth/enrichment'

working_labels = list.files(gene_sets_dir, pattern='_A_|_B_|_E_')
working_labels = working_labels[c(2, 1, 4, 3)]

#
# 2. iterate gene sets
#
all_identifiers = list()
for (working_label in working_labels) {
  
  ## 2.1. retrieve genes
  working_file = paste(gene_sets_dir, working_label, sep='/')
  df = read.csv(working_file)
  identifiers = df$NCBI.gene..formerly.Entrezgene..ID
  cat(blue(paste('working with', working_label, 'with', length(identifiers), 'genes', sep=' ')), fill=TRUE)
  
  all_identifiers[[substr(working_label, 10, 24)]] = identifiers
  
  cat('\n')
}

#
# 3. run the analysis on different Ontologies
#

ck = compareCluster(all_identifiers, fun="enrichPathway", pvalueCutoff=0.05)
sm = pairwise_termsim(ck)
dotplot(ck, size='count', showCategory=10, title='enrichPathway')
emapplot(sm)
storage = paste(output_dir, 'clusterProfiler_enrichments_RP.csv', sep='/')
write.csv(ck@compareClusterResult, storage)

ck = compareCluster(all_identifiers, fun="enrichGO", pvalueCutoff=0.05, OrgDb='org.Hs.eg.db')
sm = pairwise_termsim(ck)
dotplot(ck, size='count', showCategory=10, title='enrichGO')
emapplot(sm)

ck = compareCluster(all_identifiers, fun="enrichKEGG", pvalueCutoff=0.05)
sm = pairwise_termsim(ck)
dotplot(ck, size='count', showCategory=10, title='enrichKEGG')
emapplot(sm)

ck = compareCluster(all_identifiers, fun="enrichDO", pvalueCutoff=0.05)
sm = pairwise_termsim(ck)
dotplot(ck, size='count', showCategory=5, title='enrichDO')
emapplot(sm)


