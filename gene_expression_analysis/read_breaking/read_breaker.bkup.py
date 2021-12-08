#
# This script reads the directory of clean reads and create n sets of sampled paired FASTQ files
#

import gzip, datetime, numpy, os, sys

def breaker(label):

    #
    # 1. detect the number of lines
    #
    printt('\t detecting the number of lines')

    # detecting the number of reads in R1 file
    R1_file = reads_directory + label + '/{}_R1_clean.fastq.gz'.format(label)
    with gzip.open(R1_file, 'rb') as f:
        R1_lines = f.readlines()
    R1_number_of_lines = len(R1_lines)
    printt('\t\t detected {} lines for R1'.format(R1_number_of_lines))

    R2_file = reads_directory + label + '/{}_R2_clean.fastq.gz'.format(label)
    with gzip.open(R2_file, 'rb') as f:
        R2_lines = f.readlines()
    R2_number_of_lines = len(R2_lines)
    printt('\t\t detected {} lines for R2'.format(R2_number_of_lines))

    if R1_number_of_lines != R2_number_of_lines:
        raise ValueError('not the same number of reads in R1 and R2 files')

    #sys.exit()

    #
    # 2. define the lines to read
    #
    printt('\t sampling reads')

    number_of_reads = int(R1_number_of_lines/4)
    all_indices = numpy.arange(number_of_reads)
    numpy.random.shuffle(all_indices)
    splitted_indices = numpy.array_split(all_indices, 3)

    # generate line indices instead of read indices
    for i in range(len(splitted_indices)):

        expanded = numpy.expand_dims(splitted_indices[i], axis=1)
        expanded = expanded * 4
        new = numpy.hstack((expanded, expanded+1, expanded+2, expanded+3))
        flatten = new.flatten()
        splitted_indices[i] = flatten

    #
    # 3. create new files
    #
    printt('\t creating sampled files')

    for i in range(len(splitted_indices)):

        target = R1_file.replace('_R1_clean.fastq.gz', '_R1_subsample{}_clean.fastq.gz'.format(i))
        content = numpy.array(R1_lines)[splitted_indices[i]]
        with gzip.open(target, 'wb') as f:
            for line in content:
                f.write(line)

        target = R2_file.replace('_R2_clean.fastq.gz', '_R2_subsample{}_clean.fastq.gz'.format(i))
        content = numpy.array(R2_lines)[splitted_indices[i]]
        with gzip.open(target, 'wb') as f:
            for line in content:
                f.write(line)

    return None

def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

###
### MAIN
###

#
# 0. user-defined variables
#
reads_directory = '/home/adrian/projects/hegoi/data/rna-seq/uncompressed/'

#
# 1. read working samples
#
labels = next(os.walk(reads_directory))[1]
labels.sort()
print(labels, len(labels))

#
# 2. iterate over samples
#
for label in labels:
    printt('working with label {}'.format(label))
    breaker(label)
    sys.exit()
