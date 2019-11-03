Paula Berry
<br>Demultiplexing assignment


**Part 1: Quality Score Distribution per Nucleotide**

1.
a) Table of files

| filename | contents |
|------------|--------------|
| 1294_S1_L008_R1_001.fastq.gz | read 1 |
| 1294_S1_L008_R2_002.fastq.gz | index 1 |
| 1294_S1_L008_R2_003.fastq.gz | index 2 |
| 1294_S1_L008_R4_001.fastq.gz | read 2 |

b) Histograms:

Data is contained in [qualityhistograms.out](qualityhistograms.out). The average quality score for each file is used as the cut off.

Scripts used to generate data:<br>
[qualityhistograms.sh](qualityhistograms.sh) - script with slurm and bash commands to initialize correct python environment and call the python script.<br>
[qualityhistograms.py](qualityhistograms.py) - python script that performed the statistical analysis on the read quality and generated the histograms.

Sequence 1 average quality score cut off: 39.315685332524005

![Sequence 1 Histogram](/images/seq1_histogram.png "Sequence 1 Histogram")


Index 1 average quality score cut off: 35.86209954646942

![Index 1 Histogram](/images/index1_histogram.png "Index 1 Histogram")


Sequence 2 average quality score cut off: 37.427600426257854

![Sequence 2 Histogram](/images/seq2_histogram.png "Sequence 2 Histogram")


Index 2 average quality score cut off: 34.976062685394815

![Index 2 Histogram](/images/index2_histogram.png "Index 2 Histogram")


c) 7304664 indexes were sequenced with an "N" base call. I used n_index.sh to wrap the following command in slurm script.

```$ /usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | awk 'NR % 4 == 2' | grep '[N]\+' | wc -l > n_index_count.txt```

2. Algorithm written in pseudocode is available on [demultiplex_algorithm.txt](demultiplex_algorithm.txt).

Test input files:
* [testR1index.fastq](/testfiles/testR1index.fastq)
* [testR2index.fastq](/testfiles/testR2index.fastq)
* [testR1sequence.fastq](/testfiles/testR1sequence.fastq)
* [testR2sequence.fastq](/testfiles/testR2sequence.fastq)

Test output files:
* [testR1dualmatched.fastq](/testfiles/testR1dualmatched.fastq)
* [testR2dualmatched.fastq](/testfiles/testR2dualmatched.fastq)
* [testR1unknown.fastq](/testfiles/testR1unknown.fastq)
* [testR2unknown.fastq](/testfiles/testR2unknown.fastq)
* [testR1hopped.fastq](/testfiles/testR1hopped.fastq)
* [testR2hopped.fastq](/testfiles/testR2hopped.fastq)
