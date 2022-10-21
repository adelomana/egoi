BiocManager::install("enrichR")

#
# 0. load libraries
#
library(biomaRt)

# setting up enrichR
library(enrichR)
setEnrichrSite("Enrichr") 
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)
library(knitr)
if (websiteLive) kable(head(dbs[c(1:6),-4]))

#
# 1. read gene list
#
ensembl_file = '/home/adrian/projects/hegoi/results/functional_enrichment_sleuth/hypo_A_up.tsv'
gene_list = read.csv(ensembl_file, header=FALSE, col.names='ENSEMBL_IDs')
dim(gene_list)
View(gene_list)

#
# 2. retreive annotations to convert identifiers
#
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host='https://uswest.ensembl.org')
#! listAttributes(mart)
rosetta = biomaRt::getBM(attributes=c('ensembl_transcript_id', "ensembl_gene_id", "external_gene_name", 'entrezgene_id'), mart=mart)
rosetta$ensembl_transcript_id = NULL
dim(rosetta)
rosetta = rosetta[!duplicated(rosetta$ensembl_gene_id), ]
dim(rosetta)
View(rosetta)

#
# 3. intersect
#
working_genes = rosetta[rosetta$ensembl_gene_id %in% gene_list$ENSEMBL_IDs, ]
dim(working_genes)
View(working_genes)
length(working_genes$external_gene_name)

#
# 4. enrichment
#
dbs <- c('GO_Biological_Process_2015', 'Reactome_2016')
if (websiteLive) {
  enriched <- enrichr(working_genes$external_gene_name, dbs)
}

#
# 5. storing
#

for (db in dbs) {
  hits = enriched[[db]][enriched[[db]]$Adjusted.P.value < 0.05, ]
  print(db)
  print(dim(hits))
}
#enriched[["Reactome_2016"]][enriched[["Reactome_2016"]]$Adjusted.P.value < 0.05, ]





