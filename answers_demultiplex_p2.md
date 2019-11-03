Paula Berry
<br>Demultiplexing assignment, Part 2

**Scripts**<br>
[demultiplex.sh](/scripts/demultiplex.sh) - script with slurm and bash commands to inialize correct Python3 envrionment and call the Python script.<br>
[demultiplex.py](/scripts/demultiplex.py) - Python script that performed the analysis and demultiplexing of the reads.

The python script only analyzed the quality score of the indices, and used a cut off of 35 for the lowest quality score accepted for either R1 or R2 reads. The end user can alter the quality cut off score using the ```-q``` option in their bash command. A list of indexes used for the sequencing run is required by the ```-i``` option. An option to have the script automatically figure out which indexes are correct based on this option having an integer input indicating the number of different samples is currently in development.

Totals and percentages are output in the [runstats.txt](runstats.txt) file.

**Overall Sequencing Run Stats**
| Total reads | 363246735 |
| Total unknown/low quality reads | 104584689 (28.79%) |
| Total index-hopped reads | 27286519 (7.51%) |
| Total matched reads | 231375527 (63.70%) |

**Percentage of reads from each sample, QS cut off = 35**
| index | number of reads that passed QC | percentage of reads |
|------------|--------------|--------------|
| GTAGCGTA | 5861439 | 1.61% |
| CGATCGAT | 4315853 | 1.19% |
| GATCAAGG | 4755585 | 1.31% |
| AACAGCGA | 6414541 | 1.77% |
| TAGCCATG | 7729032 | 2.13% |
| CGGTAATC | 2987857 | 0.82% |
| CTCTGGAT | 25336662 | 6.98% |
| TACCGGAT | 48414288 | 13.33% |
| CTAGCTCA | 13365218 | 3.68% |
| CACTTCAC | 2639732 | 0.73% |
| GCTACTCT | 4588455 | 1.26% |
| ACGATCAG | 6113028 | 1.68% |
| TATGGCAC | 7563947 | 2.08% |
| TGTTCCGT | 11950887 | 3.29% |
| GTCCTAAG | 6500149 | 1.79% |
| TCGACAAG | 2689280 | 0.74% |
| TCTTCGAC | 30356308 | 8.36% |
| ATCATGCG | 7248822 | 2.00% |
| ATCGTGGT | 4865615 | 1.34% |
| TCGAGAGT | 7406313 | 2.04% |
| TCGGATTC | 3042115 | 0.84% |
| GATCTTGC | 2793834 | 0.77% |
| AGAGTCCA | 8027830 | 2.21% |
| AGGATAGC | 6408737 | 1.76% |

**Number of hopped reads from each sample, QS cut off = 35**
| index | number of hopped reads (passed QC)|
|------------|--------------|
| GTAGCGTA | 485599  |
| CGATCGAT | 300436  |
| GATCAAGG | 381425  |
| AACAGCGA | 559257  |
| TAGCCATG | 563164  |
| CGGTAATC | 279583  |
| CTCTGGAT | 2026529  |
| TACCGGAT | 7025787  |
| CTAGCTCA | 862614  |
| CACTTCAC | 324200  |
| GCTACTCT | 511275  |
| ACGATCAG | 468023  |
| TATGGCAC | 1041132  |
| TGTTCCGT | 1020483  |
| GTCCTAAG | 555121  |
| TCGACAAG | 278822  |
| TCTTCGAC | 3137482  |
| ATCATGCG | 602488  |
| ATCGTGGT | 442389  |
| TCGAGAGT | 836198  |
| TCGGATTC | 268845  |
| GATCTTGC | 208289  |
| AGAGTCCA | 681708  |
| AGGATAGC | 500525  |
