####################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biomaRt", force = TRUE)
#####################################################

###
### 0. load libraries
###
library(DESeq2)
library(tximport)
library(biomaRt)

###
### 1. user-defined variables
###
kallisto_dir = '/home/adrian/projects/hegoi/results/samples/kallisto/kallisto.100'
results_dir = '/home/adrian/projects/hegoi/results/tpm/'

###
### 2. annotation
###

# #! using biomart: no transcript missing but I end up with 40,320 genes
# working_atributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'gene_biotype', 'description')
# ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
# table96 = getBM(attributes=working_atributes, mart=ensembl96)
# dim(table96)
# View(table96)

mart = biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host='https://uswest.ensembl.org')
t2g = biomaRt::getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id'), mart=mart)
t2g = dplyr::rename(t2g, target_id=ensembl_transcript_id, ens_gene=ensembl_gene_id)

###
### 3. read files
###
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames

files = file.path(dirnames, 'abundance.h5')
files

labels = sapply(strsplit(files, split='/',fixed=TRUE), function(x) (x[10]))
labels

txi = tximport(files, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

###
### 4. find abundance
###
tpm = txi$abundance
colnames(tpm) = labels
dim(tpm)
View(tpm)

###
### 5. store
###
store = paste(results_dir, 'DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)