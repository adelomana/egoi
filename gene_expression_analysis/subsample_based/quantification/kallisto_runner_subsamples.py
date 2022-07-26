###
### usage: time python kallisto_runner_subsamples.py &> messages.txt
###

import sys, datetime, os

def kallisto_caller(label):

    printt('about to quantify {}'.format(label))

    #
    # define fastq files
    #
    found_files = os.listdir(clean_fastq_dir + label)

    all_fastq_roots = []
    for element in found_files:
        root = element.split('_R')[0]
        all_fastq_roots.append(root)
    fastq_roots = list(set(all_fastq_roots))

    modes = ['subsample0', 'subsample1', 'subsample2']
    for mode in modes:
        fastq_files = ''
        for root in fastq_roots:
            a = '{}{}/{}_R1_{}_clean.fastq'.format(clean_fastq_dir, label, root, mode)
            b = '{}{}/{}_R2_{}_clean.fastq'.format(clean_fastq_dir, label, root, mode)
            fastq_files = fastq_files + a + ' ' + b + ' '

        sample_output_dir = results_dir + label + '_' + mode
        executable = 'time kallisto quant'
        options = ' -i {} -o {} --bias -t {} -b {} {}'.format(transcriptome_index, sample_output_dir, threads, boots, strand_flag)
        command = executable + options + fastq_files

        print('')
        print(command)
        os.system(command)
        print('')

    #
    # running all samples together
    #
    fastq_files = ''
    modes = ['subsample0', 'subsample1', 'subsample2']
    for mode in modes:
        for root in fastq_roots:
            a = '{}{}/{}_R1_{}_clean.fastq'.format(clean_fastq_dir, label, root, mode)
            b = '{}{}/{}_R2_{}_clean.fastq'.format(clean_fastq_dir, label, root, mode)
            fastq_files = fastq_files + a + ' ' + b + ' '

    mode = 'all'
    sample_output_dir = results_dir + label + '_' + mode
    executable = 'time kallisto quant'
    options = ' -i {} -o {} --bias -t {} -b {} {}'.format(transcriptome_index, sample_output_dir, threads, boots, strand_flag)
    command = executable + options + fastq_files

    print('')
    print(command)
    os.system(command)
    print('')

    return None

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

#
# 0. user defined variables
#
clean_fastq_dir = '/home/adrian/projects/hegoi/data/rna-seq/subsamples/'
boots = 100
threads = 20
results_dir = '/home/adrian/projects/hegoi/results/subsamples/kallisto.{}/'.format(boots)
transcriptome_index = '/home/adrian/software/kallisto/ensembl_v96/transcriptome.idx'

strand_flag = '--rf-stranded'       # [quant] processed 45,799,936 reads, 20,315,352 reads pseudoaligned
strand_flag = '--fr-stranded'       # [quant] processed 45,799,936 reads, 20,987,918 reads pseudoaligned
strand_flag = ''                    # [quant] processed 45,799,936 reads, 41,068,344 reads pseudoaligned

#
# 1. recover labels
#
printt('recover labels...')

labels = next(os.walk(clean_fastq_dir))[1]
labels.sort()
print(labels, len(labels))

#
# 2. call kallisto quant
#
if os.path.exists(results_dir) == False:
    os.mkdir(results_dir)

for label in labels:
    kallisto_caller(label)
    #sys.exit()
