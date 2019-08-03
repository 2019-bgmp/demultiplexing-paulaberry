#!/usr/bin/env python3
import gzip
#import numpy as np
import math
import matplotlib.pyplot as plt

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
sequence1 = gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz", "rt")
index1 = gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz", "rt")
sequence2 = gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz", "rt")
index2 = gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz", "rt")

# test fast q files to process
#sequence1 = gzip.open("/mnt/c/Users/Paula/Documents/BGMP/Bi622/demultiplexing-paulaberry/test_R1.fastq.gz", "rt")
#index1 = gzip.open("/mnt/c/Users/Paula/Documents/BGMP/Bi622/demultiplexing-paulaberry/test_R2.fastq.gz", "rt")
#sequence2 = gzip.open("/mnt/c/Users/Paula/Documents/BGMP/Bi622/demultiplexing-paulaberry/test_R4.fastq.gz", "rt")
#index2 = gzip.open("/mnt/c/Users/Paula/Documents/BGMP/Bi622/demultiplexing-paulaberry/test_R3.fastq.gz", "rt")

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

sequence1.close()
index1.close()
sequence2.close()
index2.close()
