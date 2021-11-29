###
### 1. load libraries
###
library(rlist) # necessary for list appending
library(biomaRt)
library(tximport)
library(DESeq2)

# performance
library(BiocParallel)
library(tictoc)

###
### 2. user defined variables
###
register(MulticoreParam(20))

setwd("~/scratch/")

kallisto_dir = "/home/adrian/projects/hegoi/data/rna-seq/kallisto/kallisto.100"
metadata_file = "/home/adrian/projects/hegoi/metadata/hegoi metadata - hypotheses formatted for DESeq2.csv"
results_dir = '/home/adrian/projects/hegoi/results/DESeq2/'

###
### 3. read gene annotation
###
working_atributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'gene_biotype', 'description', 'hgnc_symbol')
ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
table96 = getBM(attributes=working_atributes, mart=ensembl96)
table96['common'] = table96$ensembl_gene_id
dim(table96)
View(table96)
#! the order should be first transcript as first column and gene as second column. The rest of columns are not used.

###
### 4. read metadata
###
metadata = read.table(metadata_file, sep=',', header=TRUE)
View(metadata)

###
### 5.1. run first hypothesis
###
tag = 'hypoA'
local_metadata = metadata[metadata$hypothesis == 'static vs laminar', ]
local_samples = local_metadata$sample
files = file.path(kallisto_dir, local_samples, "abundance.h5")
print(local_metadata)
print(files)
txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)
dds = DESeqDataSetFromTximport(txi, colData=local_metadata, design=~condition) 
dds$condition = relevel(dds$condition, ref="static")
dim(dds)

threshold = 10
keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]
print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))

dds = DESeq(dds, parallel=TRUE)

res = results(dds, lfcThreshold=1, parallel=TRUE) 
filt1 = res[which(res$pvalue < 0.05), ]
filt2 = filt1[which(filt1$padj < 0.1), ]
print(paste('DEGs found', dim(filt2)[1], sep=' '))
write.table(filt2, file=paste(results_dir, 'unformatted/unformatted_results_', tag, '.tsv', sep=''), quote=FALSE, sep='\t')

df = as.data.frame(filt2)
df['common'] = rownames(df)
annotation_table = table96[, c(3, 4, 5, 6)]
annotation_table_unique = annotation_table[!duplicated(annotation_table$common), ]
dn = merge(df, annotation_table_unique, by='common')
# check for not missing DEGs because of annotation
if (dim(df)[1] != dim(dn)[1]){
  print('ERROR: DEG missed on annotation step')
  stop()
}

up = dn[dn$log2FoldChange > 0, ]
down = dn[dn$log2FoldChange < 0, ]
sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
sorted_down = down[order(down$log2FoldChange), ]

store = paste(results_dir, tag, '_up', '.tsv', sep='')
write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
store = paste(results_dir, tag, '_down', '.tsv', sep='')
write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)

###
### 5.2. working on second hypothesis
###
tag = 'hypoB'
local_metadata = metadata[metadata$hypothesis == 'laminar vs oscillatory without Pi', ]
local_samples = local_metadata$sample
files = file.path(kallisto_dir, local_samples, "abundance.h5")
print(local_metadata)
print(files)
txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)
dds = DESeqDataSetFromTximport(txi, colData=local_metadata, design=~condition) 
dds$condition = relevel(dds$condition, ref="laminar")
dim(dds)

threshold = 10
keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]
print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))

dds = DESeq(dds, parallel=TRUE)

res = results(dds, lfcThreshold=1, parallel=TRUE) 
filt1 = res[which(res$pvalue < 0.05), ]
filt2 = filt1[which(filt1$padj < 0.1), ]
print(paste('DEGs found', dim(filt2)[1], sep=' '))
write.table(filt2, file=paste(results_dir, 'unformatted/unformatted_results_', tag, '.tsv', sep=''), quote=FALSE, sep='\t')

df = as.data.frame(filt2)
df['common'] = rownames(df)
annotation_table = table96[, c(3, 4, 5, 6)]
annotation_table_unique = annotation_table[!duplicated(annotation_table$common), ]
dn = merge(df, annotation_table_unique, by='common')
# check for not missing DEGs because of annotation
if (dim(df)[1] != dim(dn)[1]){
  print('ERROR: DEG missed on annotation step')
  stop()
}

up = dn[dn$log2FoldChange > 0, ]
down = dn[dn$log2FoldChange < 0, ]
sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
sorted_down = down[order(down$log2FoldChange), ]

store = paste(results_dir, tag, '_up', '.tsv', sep='')
write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
store = paste(results_dir, tag, '_down', '.tsv', sep='')
write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)

###
### 5.3. working on second hypothesis
###
tag = 'hypoC'
local_metadata = metadata[metadata$hypothesis == 'laminar vs oscillatory with Pi', ]
local_samples = local_metadata$sample
files = file.path(kallisto_dir, local_samples, "abundance.h5")
print(local_metadata)
print(files)
txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)
dds = DESeqDataSetFromTximport(txi, colData=local_metadata, design=~condition) 
dds$condition = relevel(dds$condition, ref="laminar")
dim(dds)

threshold = 10
keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]
print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))

dds = DESeq(dds, parallel=TRUE)

res = results(dds, lfcThreshold=1, parallel=TRUE) 
filt1 = res[which(res$pvalue < 0.05), ]
filt2 = filt1[which(filt1$padj < 0.1), ]
print(paste('DEGs found', dim(filt2)[1], sep=' '))
write.table(filt2, file=paste(results_dir, 'unformatted/unformatted_results_', tag, '.tsv', sep=''), quote=FALSE, sep='\t')

df = as.data.frame(filt2)
df['common'] = rownames(df)
annotation_table = table96[, c(3, 4, 5, 6)]
annotation_table_unique = annotation_table[!duplicated(annotation_table$common), ]
dn = merge(df, annotation_table_unique, by='common')
# check for not missing DEGs because of annotation
if (dim(df)[1] != dim(dn)[1]){
  print('ERROR: DEG missed on annotation step')
  stop()
}

up = dn[dn$log2FoldChange > 0, ]
down = dn[dn$log2FoldChange < 0, ]
sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
sorted_down = down[order(down$log2FoldChange), ]

store = paste(results_dir, tag, '_up', '.tsv', sep='')
write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
store = paste(results_dir, tag, '_down', '.tsv', sep='')
write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)

###
### 5.4. working on second hypothesis
###
tag = 'hypothesisD'
local_metadata = metadata[metadata$hypothesis == tag, ]
local_samples = local_metadata$sample
files = file.path(kallisto_dir, local_samples, "abundance.h5")
print(local_metadata)
print(files)
txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)
dds = DESeqDataSetFromTximport(txi, colData=local_metadata, design=~condition) 
dds$condition = relevel(dds$condition, ref="no")
dim(dds)

threshold = 10
keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]
print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))

dds = DESeq(dds, parallel=TRUE)

res = results(dds, lfcThreshold=1, parallel=TRUE) 
filt1 = res[which(res$pvalue < 0.05), ]
filt2 = filt1[which(filt1$padj < 0.1), ]
print(paste('DEGs found', dim(filt2)[1], sep=' '))
write.table(filt2, file=paste(results_dir, 'unformatted/unformatted_results_', tag, '.tsv', sep=''), quote=FALSE, sep='\t')

df = as.data.frame(filt2)
df['common'] = rownames(df)
annotation_table = table96[, c(3, 4, 5, 6)]
annotation_table_unique = annotation_table[!duplicated(annotation_table$common), ]
dn = merge(df, annotation_table_unique, by='common')
# check for not missing DEGs because of annotation
if (dim(df)[1] != dim(dn)[1]){
  print('ERROR: DEG missed on annotation step')
  stop()
}

up = dn[dn$log2FoldChange > 0, ]
down = dn[dn$log2FoldChange < 0, ]
sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
sorted_down = down[order(down$log2FoldChange), ]

store = paste(results_dir, tag, '_up', '.tsv', sep='')
write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
store = paste(results_dir, tag, '_down', '.tsv', sep='')
write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)


###
### 5.5. working on second hypothesis
###
tag = 'hypothesisE'
local_metadata = metadata[metadata$hypothesis == tag, ]
local_samples = local_metadata$sample
files = file.path(kallisto_dir, local_samples, "abundance.h5")
print(local_metadata)
print(files)
txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)
dds = DESeqDataSetFromTximport(txi, colData=local_metadata, design=~condition) 
dds$condition = relevel(dds$condition, ref="no")
dim(dds)

threshold = 10
keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]
print(paste('dimensions of filtered genes', dim(dds)[1], dim(dds)[2]))

dds = DESeq(dds, parallel=TRUE)

res = results(dds, lfcThreshold=1, parallel=TRUE) 
filt1 = res[which(res$pvalue < 0.05), ]
filt2 = filt1[which(filt1$padj < 0.1), ]
print(paste('DEGs found', dim(filt2)[1], sep=' '))
write.table(filt2, file=paste(results_dir, 'unformatted/unformatted_results_', tag, '.tsv', sep=''), quote=FALSE, sep='\t')

df = as.data.frame(filt2)
df['common'] = rownames(df)
annotation_table = table96[, c(3, 4, 5, 6)]
annotation_table_unique = annotation_table[!duplicated(annotation_table$common), ]
dn = merge(df, annotation_table_unique, by='common')
# check for not missing DEGs because of annotation
if (dim(df)[1] != dim(dn)[1]){
  print('ERROR: DEG missed on annotation step')
  stop()
}

up = dn[dn$log2FoldChange > 0, ]
down = dn[dn$log2FoldChange < 0, ]
sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
sorted_down = down[order(down$log2FoldChange), ]

store = paste(results_dir, tag, '_up', '.tsv', sep='')
write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
store = paste(results_dir, tag, '_down', '.tsv', sep='')
write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)














