#
# This script reads the directory of clean reads and create n sets of sampled paired FASTQ files
#

import datetime, numpy, os, sys

def breaker(label):

    #
    # 1. detect the number of lines
    #
    printt('\t detecting the number of lines')

    # detecting the number of reads in R1 file
    R1_file = reads_directory + label + '/{}_R1_clean.fastq'.format(label)
    f = open(R1_file, 'r')
    R1_lines = f.readlines()
    f.close()
    R1_number_of_lines = len(R1_lines)
    printt('\t\t detected {} lines for R1'.format(R1_number_of_lines))

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
        printt('\t\t subsample {}'.format(i))

        target = R1_file.replace('_R1_clean.fastq', '_R1_subsample{}_clean.fastq'.format(i))
        g = open(target, 'w')
        for index in splitted_indices[i]:
            g.write(R1_lines[index])
        g.close()

    # now working with R2 from scratch
    printt('\t\t working wih R2')
    del(R1_lines)

    printt('\t\t\t reading R2')
    R2_file = reads_directory + label + '/{}_R2_clean.fastq'.format(label)
    f = open(R2_file, 'r')
    R2_lines = f.readlines()
    f.close()

    printt('\t\t\t saving R2 samples')
    for i in range(len(splitted_indices)):
        printt('\t\t\t\t subsample {}'.format(i))

        target = R2_file.replace('_R2_clean.fastq', '_R2_subsample{}_clean.fastq'.format(i))
        g = open(target, 'w')
        for index in splitted_indices[i]:
            g.write(R2_lines[index])
        g.close()

    del(R2_lines)

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
reads_directory = '/home/adrian/scratch/'

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
