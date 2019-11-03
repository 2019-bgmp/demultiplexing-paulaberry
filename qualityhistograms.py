#!/usr/bin/env python3
import gzip
import numpy as np
import math
import matplotlib.pyplot as plt
import argparse

def get_args():
    """Function to pass in filenames to demultiplex."""
    getfiles = argparse.ArgumentParser(description="A program demultiplex dual-indexed fastq file reads.")
    getfiles.add_argument("-iR1", "--indexR1", help = "read 1 index fastq file", required = True)
    getfiles.add_argument("-iR2", "--indexR2", help = "read 2 index fastq file", required = True)
    getfiles.add_argument("-sR1", "--sequenceR1", help = "read 1 sequence fastq file", required = True)
    getfiles.add_argument("-sR2", "--sequenceR2", help = "read 2 sequence fastq file", required = True)
    return getfiles.parse_args()
args = get_args()

R1index = str(args.indexR1)
R2index = str(args.indexR2)
R1sequence = str(args.sequenceR1)
R2sequence = str(args.sequenceR2)

def convert_phred(letter):
    '''A function to convert a Phred 33 score to a quality score value.'''
    qscore = ord(letter) - 33
    return qscore

def read_length(filename):
    '''A function that scans a fastq file and returns the length of the reads.'''
    counter = 0
    for line in filename:
        print(line)
        counter = counter + 1
        if counter == 4:
            line = line.strip()
            read_length = len(line)
            return int(read_length)

def qscore_means(filename, dictionary):
    '''A function that takes a fastq file, and returns a dictionary of each nucleotide position as keys with values of mean qscore for that position.'''
    line_counter = 0
    read_counter = 0
    for line in filename:
        line_counter = line_counter + 1
        if line_counter % 4 == 0: # extracts line of quality scores
            read_counter = read_counter + 1
            line = line.strip()
            nuc_position = 0
            for character in line: # creates running sums in the dictionary
                nuc_position = nuc_position + 1
                qscore = convert_phred(character)
                if nuc_position not in dictionary.keys():
                    dictionary[nuc_position] = qscore
                else:
                    dictionary[nuc_position] = dictionary[nuc_position] + qscore
    for key in dictionary:
        dictionary[key] = dictionary[key] / read_counter

def create_histogram(dictionary, title, figurename):
    '''takes a dictionary of nucleotide positions and mean quality scores and returns a histogram.'''
    x = list(dictionary.keys())
    y = list(dictionary.values())
    plt.bar(x, y)
    #plt.xticks(range(len(dictionary)), list(dictionary.keys()))
    plt.xlabel("Nucleotide Position")
    plt.ylabel("Average Quality Score")
    plt.title(str(title))
    plt.savefig(figurename + ".png")



# actual fast q files to process
#sequence1 = gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz", "rt")
#index1 = gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz", "rt")
#sequence2 = gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz", "rt")
#index2 = gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz", "rt")

index1 = gzip.open(R1index, "rt")
index2 = gzip.open(R2index, "rt")
sequence1 = gzip.open(R1sequence, "rt")
sequence2 = gzip.open(R2sequence, "rt")

# qscore lists for each file
sequence1_qscores = {}
index1_qscores = {}
sequence2_qscores = {}
index2_qscores = {}

# populate the dictionaries of mean quality scores
qscore_means(sequence1, sequence1_qscores)
qscore_means(index1, index1_qscores)
qscore_means(sequence2, sequence2_qscores)
qscore_means(index2, index2_qscores)

#print(index1_qscores)

# create histogram images for each dictionary
create_histogram(sequence1_qscores, "Sequence 1 Quality Scores Histogram", "seq1_histogram")
plt.savefig:("sequence1.png")
plt.close()

create_histogram(sequence2_qscores, "Sequence 2 Quality Scores Histogram", "seq2_histogram")
plt.savefig:("sequence2.png")
plt.close()

create_histogram(index1_qscores, "Index 1 Quality Scores Histogram", "index1_histogram")
plt.savefig:("index1.png")
plt.close()

create_histogram(index2_qscores, "Index 2 Quality Scores Histogram", "index2_histogram")
plt.savefig:("index2.png")
plt.close()

print("Sequence 1 Histogram Values")
print("Nucleotide position  Mean quality score")
s1_counter = 0
s1_sum = 0
for key in sequence1_qscores:
    print(str(key) + "  " + str(sequence1_qscores[key]))
    s1_sum = s1_sum + sequence1_qscores[key]
    s1_counter = s1_counter + 1
s1_mean = s1_sum / s1_counter
print("Average quality score for Sequence 1: " + str(s1_mean))
print()

print("Index 1 Histogram Values")
print("Nucleotide position  Mean quality score")
i1_counter = 0
i1_sum = 0
for key in index1_qscores:
    print(str(key) + "  " + str(index1_qscores[key]))
    i1_sum = i1_sum + index1_qscores[key]
    i1_counter = i1_counter + 1
i1_mean = i1_sum / i1_counter
print("Average quality score for Index 1: " + str(i1_mean))
print()

print("Sequence 2 Histogram Values")
print("Nucleotide position  Mean quality score")
s2_counter = 0
s2_sum = 0
for key in sequence2_qscores:
    print(str(key) + "  " + str(sequence2_qscores[key]))
    s2_sum = s2_sum + sequence2_qscores[key]
    s2_counter = s2_counter + 1
s2_mean = s2_sum / s2_counter
print("Average quality score for Sequence 2: " + str(s2_mean))
print()

print("Index 2 Histogram Values")
print("Nucleotide position  Mean quality score")
i2_counter = 0
i2_sum = 0
for key in index2_qscores:
    print(str(key) + "  " + str(index2_qscores[key]))
    i2_sum = i2_sum + index2_qscores[key]
    i2_counter = i2_counter + 1
i2_mean = i2_sum / i2_counter
print("Average quality score for Index 2: " + str(i2_mean))
print()

sequence1.close()
index1.close()
sequence2.close()
index2.close()
