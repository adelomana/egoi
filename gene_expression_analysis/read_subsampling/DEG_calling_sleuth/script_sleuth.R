library(sleuth)
library(biomaRt)

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/home/adrian/projects/hegoi/results/subsamples/kallisto.100"
metadata_file = "/home/adrian/projects/hegoi/metadata/hegoi metadata - hypotheses formatted for DESeq2 subsample.tsv"
results_dir = '/home/adrian/projects/hegoi/results/subsamples/sleuth/'

subtags = c('subsample0', 'subsample1', 'subsample2')

#
# 1. generate gene to transcript mapping
#
# mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                          dataset = "hsapiens_gene_ensembl",
#                          host = 'ensembl.org')
# t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
#                                      "external_gene_name"), mart = mart)
# t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
#                      ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# working_atributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name")
# ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
# t2g = getBM(attributes=working_atributes, mart=ensembl96)
# t2g = dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
# dim(t2g)
# View(t2g)

#
# 2. read metadata
#
metadata = read.table(metadata_file, sep='\t', header=TRUE)
View(metadata)

#
# 3. iterate hypotheses
#
for (hypothesis in 1:nrow(metadata)) {
  tag = metadata[hypothesis, "hypothesis"]
  sampleA = metadata[hypothesis, "sampleA"]
  sampleB = metadata[hypothesis, "sampleB"]
  print(c(sampleA, sampleB))
  
}

















