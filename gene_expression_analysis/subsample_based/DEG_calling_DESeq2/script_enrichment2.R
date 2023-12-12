BiocManager::install("enrichR")
BiocManager::install("enrichplot")

#
# 0. load libraries
#
library(crayon) # so messages can be in blue font
library(biomaRt)
library(enrichplot)

# setting up enrichR
library(enrichR)
setEnrichrSite("Enrichr") 
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)
library(knitr)
if (websiteLive) kable(head(dbs[c(1:6),-4]))


# 0. user-defined variables
gene_sets_dir = '/home/adrian/projects/hegoi/results/subsamples/DEG_filtered'
output_dir = '/home/adrian/projects/hegoi/results/subsamples/functional_enrichment'

working_labels = list.files(gene_sets_dir, pattern='_A_|_B_|_E_')
dbs <- c('GO_Biological_Process_2015', 'Reactome_2016')

# 1. iterate gene sets
for (working_label in working_labels) {
  
  ## 1.1. retrieve genes
  working_file = paste(gene_sets_dir, working_label, sep='/')
  df = read.csv(working_file)
  genes = df$Symbol
  cat(blue(paste('working with', working_label, 'with', length(genes), 'genes', sep=' ')), fill=TRUE)
  
  ## 1.2. enrichment
  if (websiteLive) {
    enriched <- enrichr(genes, dbs)
  }
  for (db in dbs) {
    hits = enriched[[db]][enriched[[db]]$Adjusted.P.value < 0.05, ]
    cat(blue(paste('\t found', dim(hits)[1], 'significant terms in', db, sep=' ')), fill=TRUE)
    
    # storage
    file_name = paste('enrichments', db, strsplit(working_label, '_')[[1]][3], strsplit(working_label, '_')[[1]][4], sep='_')
    storage = paste(output_dir, file_name, sep='/')
    write.table(hits, storage, quote=FALSE, sep='\t')
  }
  
  ## 1.3. store the ensembl IDs for manual check on AmiGO2 and enriched functions
  file_name = paste('ensembl', strsplit(working_label, '_')[[1]][3], strsplit(working_label, '_')[[1]][4], sep='_')
  storage = paste(output_dir, file_name, sep='/')
  write(df$ENSEMBL, storage)

  cat('\n')
}
