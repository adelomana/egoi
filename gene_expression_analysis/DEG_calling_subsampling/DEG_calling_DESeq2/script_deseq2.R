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

kallisto_dir = "/home/adrian/projects/hegoi/results/subsamples/kallisto.100"
metadata_file = "/home/adrian/projects/hegoi/metadata/hegoi metadata - hypotheses formatted for DESeq2 subsample.tsv"
results_dir = '/home/adrian/projects/hegoi/results/subsamples/DESeq2/'

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
metadata = read.table(metadata_file, sep='\t', header=TRUE)
View(metadata)

###
### 5. hypothesis testing
###

subtags = c('subsample0', 'subsample1', 'subsample2')

for (hypothesis in 1:nrow(metadata)) {
  tag = metadata[hypothesis, "hypothesis"]
  sampleA = metadata[hypothesis, "sampleA"]
  sampleB = metadata[hypothesis, "sampleB"]
  print(c(sampleA, sampleB))
  
  # define the metadata table
  info = data.frame(hypothesis=rep(tag, 6), samplename=c(paste(sampleA, subtags, sep='_'), paste(sampleB, subtags, sep='_')), condition=c(rep('alpha', 3), rep('beta', 3)))
  print('info')
  print(info)
  
  # define files
  files = file.path(kallisto_dir, info$samplename, "abundance.h5")
  print('files')
  print(files)
  
  # read files
  txi = tximport(files, type="kallisto", tx2gene=table96, ignoreTxVersion=TRUE)
  dds = DESeqDataSetFromTximport(txi, colData=info, design=~condition) 
  dds$condition = relevel(dds$condition, ref="alpha")
  print(paste('original dimensions', dim(dds)[1], dim(dds)[2]))

  # preliminary filter on poorly detected 
  threshold = 10
  keep = rowMaxs(counts(dds)) >= threshold
  dds = dds[keep, ]
  print(paste('dimensions of filtered transcripts', dim(dds)[1], dim(dds)[2]))
  
  # run the model
  dds = DESeq(dds, parallel=TRUE)
  
  # store original results
  res = results(dds, lfcThreshold=1, parallel=TRUE) 
  filt1 = res[which(res$pvalue < 0.05), ]
  filt2 = filt1[which(filt1$padj < 0.1), ]
  print(paste('DEGs found', dim(filt2)[1], sep=' '))
  write.table(filt2, file=paste(results_dir, 'unformatted/unformatted_results_', tag, '.tsv', sep=''), quote=FALSE, sep='\t')

  # store formatted results
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

  print('---')
}

















