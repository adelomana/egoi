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
kallisto_dir = "/home/adrian/projects/flow/data/rna-seq/kallisto/kallisto.100"
results_dir = '/home/adrian/projects/flow/data/rna-seq/tpm/'

###
### 2. annotation
###

#! some functions I used to clear cache: 
#! biomartCacheInfo()
#! biomartCacheClear()

#! using biomart: no transcript missing but I end up with 40,320 genes
working_atributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'gene_biotype', 'description')
ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
table96 = getBM(attributes=working_atributes, mart=ensembl96)
dim(table96)
View(table96)

###
### 3. read files
###
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames

files = file.path(dirnames, 'abundance.h5')
files

labels = sapply(strsplit(files, split='/',fixed=TRUE), function(x) (x[10]))
labels

txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)

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