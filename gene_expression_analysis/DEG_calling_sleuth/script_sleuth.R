library(sleuth)
library(biomaRt)

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/home/adrian/projects/flow/data/rna-seq/kallisto/kallisto.100"
metadata_file = "/home/adrian/projects/flow/metadata/hegoi metadata - hypotheses formatted for sleuth.csv"

#
# 1. generate gene to transcript mapping
#
mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'ensembl.org', verbose = TRUE)
t2g = biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g = dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

#
# 2. read metadata
#
metadata <- read.table(metadata_file, sep=',', header=TRUE)
View(metadata)

#
# 3. define hypothesis object
#
s2cA = metadata[metadata$hypothesis == 'static vs laminar', c('sample', 'condition')]
s2cA$path = file.path(kallisto_dir, s2cA$sample)
View(s2cA)

#
# 4. prepare object for sleuth
#
soA = sleuth_prep(s2cA, target_mapping = t2g, aggregation_column = 'ens_gene')

#
# 5. build full and partial
#
soA = sleuth_fit(soA, ~condition, 'full')
soA = sleuth_fit(soA, ~1, 'reduced')
soA = sleuth_lrt(soA, 'reduced', 'full')
models(soA)

#
# 6. tables
#
sleuth_tableA = sleuth_results(soA, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE)
sleuth_significantA = dplyr::filter(sleuth_tableA, qval <= 0.05)
dim(sleuth_significantA)
write.csv(sleuth_significantA, 'significantA.csv')

#
# 7. go over all hypotheses
#

# 7.1. hypothesis B
s2cB = metadata[metadata$hypothesis == 'laminar vs oscillatory without Pi', c('sample', 'condition')]
s2cB$path = file.path(kallisto_dir, s2cB$sample)
View(s2cB)
soB = sleuth_prep(s2cB, target_mapping = t2g, aggregation_column = 'ens_gene')
soB = sleuth_fit(soB, ~condition, 'full')
soB = sleuth_fit(soB, ~1, 'reduced')
soB = sleuth_lrt(soB, 'reduced', 'full')
sleuth_tableB = sleuth_results(soB, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE)
sleuth_significantB = dplyr::filter(sleuth_tableB, qval <= 0.05)
dim(sleuth_significantB)
write.csv(sleuth_significantB, 'significantB.csv')

# 7.2. hypothesis C
s2cC = metadata[metadata$hypothesis == 'laminar vs oscillatory with Pi', c('sample', 'condition')]
s2cC$path = file.path(kallisto_dir, s2cC$sample)
View(s2cC)
soC = sleuth_prep(s2cC, target_mapping = t2g, aggregation_column = 'ens_gene')
soC = sleuth_fit(soC, ~condition, 'full')
soC = sleuth_fit(soC, ~1, 'reduced')
soC = sleuth_lrt(soC, 'reduced', 'full')
sleuth_tableC = sleuth_results(soC, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE)
sleuth_significantC = dplyr::filter(sleuth_tableC, qval <= 0.05)
dim(sleuth_significantC)
write.csv(sleuth_significantC, 'significantC.csv')

# 7.3. hypothesis D
s2cD = metadata[metadata$hypothesis == 'Pi or not Pi in laminar', c('sample', 'condition')]
s2cD$path = file.path(kallisto_dir, s2cD$sample)
View(s2cD)
soD = sleuth_prep(s2cD, target_mapping = t2g, aggregation_column = 'ens_gene')
soD = sleuth_fit(soD, ~condition, 'full')
soD = sleuth_fit(soD, ~1, 'reduced')
soD = sleuth_lrt(soD, 'reduced', 'full')
sleuth_tableD = sleuth_results(soD, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE)
sleuth_significantD = dplyr::filter(sleuth_tableD, qval <= 0.05)
dim(sleuth_significantD)
write.csv(sleuth_significantD, 'significantD.csv')

# 7.4. hypothesis E
s2cE = metadata[metadata$hypothesis == 'Pi or not Pi in oscillatory', c('sample', 'condition')]
s2cE$path = file.path(kallisto_dir, s2cE$sample)
View(s2cE)
soE = sleuth_prep(s2cE, target_mapping = t2g, aggregation_column = 'ens_gene')
soE = sleuth_fit(soE, ~condition, 'full')
soE = sleuth_fit(soE, ~1, 'reduced')
soE = sleuth_lrt(soE, 'reduced', 'full')
sleuth_tableE = sleuth_results(soE, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = TRUE)
sleuth_significantE = dplyr::filter(sleuth_tableE, qval <= 0.05)
dim(sleuth_significantE)
write.csv(sleuth_significantE, 'significantE.csv')






















