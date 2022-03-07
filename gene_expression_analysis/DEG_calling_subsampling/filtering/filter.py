###
### This script generates a file with filtered DEGs derived from patient-specific DEG calling.
###

"""
Filters are:
0. Filter on P < 0.05 and Q < 0.1 for DESeq2; Q < 0.05 in sleuth.
1. Filter out genes that in the comparison do not cross the 2 TPM barrier.
2. Discrete abs log2 FC > 1
3. Filter DEGs that are not consistent in n-1 patients.
"""

import pandas

def expression_reader():

    expression = pandas.read_csv(tpm_file, sep='\t', index_col=0)
    #sample_names = expression.

    return expression, sample_names

###
### MAIN
###

#
# 0. user-defined variables
#
tpm_file = '/home/adrian/projects/hegoi/results/tpm/DESeq2_TPM_values.tsv'
outputdir = '/home/adrian/projects/hegoi/results/subsamples/DEG_filtered/'

#
# 1. read expression data
#
expression, sample_names = expression_reader()

#
# 2. iterate over hypotheses and patients
#
for hypothesis in hypotheses:

    ### define filtered DEGs across patients
    DEGs_across_patients = []
    for patient in patients:

        ## read DEGs
        approach = 'DESeq2'
        DEG_up, DEG_down = DEG_reader(approach)

        ## check the filters
        filtered_DEGs = []
        for DEGs in [DEG_up, DEG_down]:
            filtered_set = DEG_filter(DEGs)
            filtered_DEGs.append(filtered_set)

        ## append filtered DEGs into patient DEG list
        DEGs_across_patients.append(filtered_DEGs)

    ### find consistent set of genes across patients
    consistent_DEGs = consistency_check(DEGs_across_patients)

    ### generate a final table
    for trend in consistent_DEGs:
        table_generator(trend)
