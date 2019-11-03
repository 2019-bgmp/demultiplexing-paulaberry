#!/usr/bin/env python3
import gzip
import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import os

def get_args():
    """Function to pass in filenames to demultiplex."""
    getfiles = argparse.ArgumentParser(description="A program demultiplex dual-indexed fastq file reads.")
    getfiles.add_argument("-iR1", "--indexR1", help = "read 1 index fastq file", required = True)
    getfiles.add_argument("-iR2", "--indexR2", help = "read 2 index fastq file", required = True)
    getfiles.add_argument("-sR1", "--sequenceR1", help = "read 1 sequence fastq file", required = True)
    getfiles.add_argument("-sR2", "--sequenceR2", help = "read 2 sequence fastq file", required = True)
    getfiles.add_argument("-i", "--indexlist", help = "file containing list of indexes used in the sequencing run", required = True)
    getfiles.add_argument("-q", "--qualitycutoff", help = "quality score cutoff for indexes", required = True)
    return getfiles.parse_args()
args = get_args()

R1index = str(args.indexR1)
R2index = str(args.indexR2)
R1sequence = str(args.sequenceR1)
R2sequence = str(args.sequenceR2)
indices = str(args.indexlist)
quality_cutoff = int(args.qualitycutoff)

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

def qscore_means(index1_line4, index2_line4):
    '''A function that takes a fastq file, and returns a dictionary of each nucleotide position as keys with values of mean qscore for that position.'''
    nuc_position = 0
    index1_qscore = 0
    index2_qscore = 0
    for character in index1_line4: # creates running sum of quality score of line
        index1_qscore = convert_phred(character) + index1_qscore
    for character in index2_line4:
        index2_qscore = convert_phred(character) + index2_qscore
    index1q_mean = index1_qscore / len(index1_line4)
    index2q_mean = index2_qscore / len(index2_line4)
    return min(index1q_mean, index2q_mean)

def importindices(indexfile, dictionary):
    """reads in a text file of a list of indices used in the sequencing run, generates a dictionary of indexes and reverse complements"""
    f = open(indexfile, "r")
    complements = {"A":"T", "G":"C", "C":"G", "T":"A", "N":"N"}
    for line in indexfile:
        line = line.strip()
        char_count = 0
        for char in line:
            if char_count == 0:
                dictionary[line] = complements[char]
            else:
                dictionary[line] = str(complements[char]) + str(dictionary[line])
    f.close()

def generate_index_dict(indexes, dictionary):
    """Takes a list of indexes as keys to a dictionary and then generates values for those keys with the reverse complement to match the second read."""
    line_counter = 0
    complements = {"A":"T", "G":"C", "C":"G", "T":"A", "N":"N"}
    for line in indexes:
        line = line.strip()
        revcomp = str("")
        for char in line:
            revcomp = str(complements[char]) + str(revcomp)
        dictionary[line] = str(revcomp)

def index_match_check(dictionary, index1_line2, index2_line2):
  """Takes sequence reads from the index files and returns True if the indexes match, False if they do not.
  Example: input of "AAAAA" and "TTTTT" outputs True, input of "CCCCC" and "GGTGG" returns False"""
  if index1_line2 in dictionary:
      if dictionary[index1_line2] == index2_line2:
          return True
      else:
        return False
  else:
      return 0

def index_N_check(index_dict, index1_line2, index2_line2):
  """Takes sequence reads from the index files and returns False if neither has any N base calls, True otherwise.
  Example: input of "AAAAA" and "TTNTT" outputs False, input of "AAAAA" and "TTTTT" outputs True.
  return True or False"""
  indexsequences = str(index1_line2) + str(index2_line2)
  if indexsequences.find("N") == -1:
      return True
  else:
      return False

def dual_index(index1, index2, sequence1):
  """Takes indices and paired-end reads and adds the index sequences to the headers of both reads
  copy both directly to match the index dictionary key:value pairs
  return R1_header, R2_header"""
  new_header = str(sequence1) + "::" + str(index1) + ":" + str(index2)
  return new_header, new_header

# variables for files we are working with
index1 = gzip.open(R1index, "rt")
index2 = gzip.open(R2index, "rt")
sequence1 = gzip.open(R1sequence, "rt")
sequence2 = gzip.open(R2sequence, "rt")
index_list = open(indices, "r")

# create empty dictionary for indexes(keys) and reverse complements(values)
index_dict = {}
generate_index_dict(index_list, index_dict)
index_list.seek(0)

# create empty directories for output files by index
#path = str(os.getcwd())
for line in index_list:
    line = line.strip()
    path = str(os.getcwd()) + "/" + str(line)
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s" % path)
index_list.seek(0)

# create an empty dictionary for indexes(keys) and filenames for barcodes
# create empty dictionaries for filenames and filehandles
file_list = {}
fh_list_R1 = {}
fh_list_R2 = {}

while True:
    barcode = str(index_list.readline().strip())
    if barcode == "":
        break
    file_list[barcode] = [barcode + "/" + barcode + "_R1.fq", barcode + "/" + barcode + "_R2.fq"]
    fh_list_R1[barcode] = barcode + "/" + barcode + "_R1"
    fh_list_R2[barcode] = barcode + "/" + barcode + "_R2"
index_list.seek(0)

# create empty R1 and R2 FASTQ files for each barcode, leave open
#for key in file_list:
    #with open(file_list[key][0], "a+") as fh_list_R1[key]:
        #continue
        #fh_list_R1[key].open()
    #with open(file_list[key][1], "a+") as fh_list_R2[key]:
        #continue
        #fh_list_R2[key].open()

# open the files we just created
#for key in fh_list_R1:
    #filename = str(fh_list_R1[key])
    #filename.open()

#for key in fh_list_R2:
    #fh_list_R2[key].open()

# create empty files to hold the demultiplexed reads and stats
unknownR1 = open("unknown_reads_R1.fq","a+")
unknownR2 = open("unknown_reads_R2.fq","a+")
indexhoppedR1 = open("hopped_reads_R1.fq","a+")
indexhoppedR2 = open("hopped_reads_R2.fq","a+")


# Open 4 files one line at a time, keep count of total sequences
sequence_count = 0
unknown_count = 0
hopped_count = 0
matched_count = 0
hopped_dictionary = {}  # empty dictionary to hold number of hopped reads
matched_dictionary = {} # empty dictionary to hold number of reads per index
complements = {"A":"T", "G":"C", "C":"G", "T":"A", "N":"N"} # global dictionary to generate reverse complements

for line in index_list: # initialize counter dictionaries to 0
    line = line.strip()
    hopped_dictionary[line] = 0
    matched_dictionary[line] = 0

while True:
    index1_line1 = index1.readline().strip()
    index1_line2 = index1.readline().strip()
    index1_line3 = index1.readline().strip()
    index1_line4 = index1.readline().strip()
    index2_line1 = index2.readline().strip()
    index2_line2 = index2.readline().strip()
    index2_line3 = index2.readline().strip()
    index2_line4 = index2.readline().strip()
    sequence1_line1 = sequence1.readline().strip()
    sequence1_line2 = sequence1.readline().strip()
    sequence1_line3 = sequence1.readline().strip()
    sequence1_line4 = sequence1.readline().strip()
    sequence2_line1 = sequence2.readline().strip()
    sequence2_line2 = sequence2.readline().strip()
    sequence2_line3 = sequence2.readline().strip()
    sequence2_line4 = sequence2.readline().strip()

    if index1_line1 == "": # stop while loop at end of index1 file
        break

    sequence_count = sequence_count + 1 # increment total sequence counter
    sequence1_line1, sequence2_line1 = dual_index(index1_line2, index2_line2, sequence1_line1) # make new sequence header

    if index_N_check(index_dict, index1_line2, index2_line2) == False: # check for uncalled bases
        unknownR1.write(sequence1_line1+"\n")
        unknownR1.write(sequence1_line2+"\n")
        unknownR1.write(sequence1_line3+"\n")
        unknownR1.write(sequence1_line4+"\n")
        unknownR2.write(sequence2_line1+"\n")
        unknownR2.write(sequence2_line2+"\n")
        unknownR2.write(sequence2_line3+"\n")
        unknownR2.write(sequence2_line4+"\n")

        unknown_count = unknown_count + 1 # increment unknown counter

    elif index_match_check(index_dict, index1_line2, index2_line2) == False:
        indexhoppedR1.write(sequence1_line1+"\n")
        indexhoppedR1.write(sequence1_line2+"\n")
        indexhoppedR1.write(sequence1_line3+"\n")
        indexhoppedR1.write(sequence1_line4+"\n")
        if index1_line2 in hopped_dictionary:
            hopped_dictionary[index1_line2] = hopped_dictionary[index1_line2] + 1
        else:
            hopped_dictionary[index1_line2] = 1
        indexhoppedR2.write(sequence2_line1+"\n")
        indexhoppedR2.write(sequence2_line2+"\n")
        indexhoppedR2.write(sequence2_line3+"\n")
        indexhoppedR2.write(sequence2_line4+"\n")
        rev_comp = ""
        for char in index2_line2:
            rev_comp = complements[char] + rev_comp
            if rev_comp not in hopped_dictionary:
                hopped_dictionary[rev_comp] = 1
            else:
                hopped_dictionary[rev_comp] = hopped_dictionary[rev_comp] + 1

        hopped_count = hopped_count + 1 # increment hopped index counter

    elif qscore_means(index1_line4, index2_line4) <= quality_cutoff:
        unknownR1.write(sequence1_line1+"\n")
        unknownR1.write(sequence1_line2+"\n")
        unknownR1.write(sequence1_line3+"\n")
        unknownR1.write(sequence1_line4+"\n")
        unknownR2.write(sequence2_line1+"\n")
        unknownR2.write(sequence2_line2+"\n")
        unknownR2.write(sequence2_line3+"\n")
        unknownR2.write(sequence2_line4+"\n")

        unknown_count = unknown_count + 1 # increment unknown counter

    elif index_match_check(index_dict, index1_line2, index2_line2) == True: # for correctly matched, high quality index reads
        R1_file = open(str(file_list[index1_line2][0]), "a+")
        R1_file.write(sequence1_line1+"\n")
        R1_file.write(sequence1_line2+"\n")
        R1_file.write(sequence1_line3+"\n")
        R1_file.write(sequence1_line4+"\n")
        R2_file = open(str(file_list[index1_line2][1]), "a+")
        R2_file.write(sequence2_line1+"\n")
        R2_file.write(sequence2_line2+"\n")
        R2_file.write(sequence2_line3+"\n")
        R2_file.write(sequence2_line4+"\n")

        matched_count = matched_count + 1 # increment total matched counter
        if index1_line2 in matched_dictionary:
            matched_dictionary[index1_line2] = matched_dictionary[index1_line2] + 1 # increment index matched counter
        else:
            matched_dictionary[index1_line2] = 1
    else: # anything that doesn't have matched or hopped indices
        unknownR1.write(sequence1_line1+"\n")
        unknownR1.write(sequence1_line2+"\n")
        unknownR1.write(sequence1_line3+"\n")
        unknownR1.write(sequence1_line4+"\n")
        unknownR2.write(sequence2_line1+"\n")
        unknownR2.write(sequence2_line2+"\n")
        unknownR2.write(sequence2_line3+"\n")
        unknownR2.write(sequence2_line4+"\n")

        unknown_count = unknown_count + 1 # increment unknown counter

# run statistics output
with open("runstats.txt", "w+") as runstats:
    runstats.write("Overall Sequencing Run Stats" + "\n")
    runstats.write("Total reads: " + str(sequence_count) + "\n")
    runstats.write("Total unknown/low quality reads: " + str(unknown_count) + " (" + str((int(unknown_count) / int(sequence_count)) * 100) + "%)" + "\n")
    runstats.write("Total index-hopped reads: " + str(hopped_count) + " (" + str((hopped_count / sequence_count) * 100) + "%)" + "\n")
    runstats.write("Total matched reads: " + str(matched_count) + " (" + str((matched_count / sequence_count) * 100) + "%)" + "\n")
    runstats.write("\n")
    runstats.write("\n")
    runstats.write("Individual file index quality stats:" + "\n")
    for key in matched_dictionary:
        runstats.write(key + ": " + str(matched_dictionary[key]) + " reads (" + str((matched_dictionary[key] / sequence_count) * 100) + "%)" + "\n")
        if key in hopped_dictionary:
            runstats.write(key + ": " + str(hopped_dictionary[key]) + " indices hopped." + "\n")

# produce heatmap of hopped indices



sequence1.close()
index1.close()
sequence2.close()
index2.close()
index_list.close()
unknownR1.close()
unknownR2.close()
indexhoppedR1.close()
indexhoppedR2.close()
runstats.close()
#for key in file_list:
    #file_list[key][0].close()
    #file_list[key][1].close()
