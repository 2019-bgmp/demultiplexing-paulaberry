#!/usr/bin/env bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=demultiplex
#SBATCH --output=demultiplex.out
#SBATCH --error=demultiplex.err
#SBATCH --time=0-12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=pberry@uoregon.edu
#SBATCH --mail-type=ALL

conda deactivate
conda deactivate
conda deactivate
conda activate bgmp_py3


#./demultiplex.py -iR1 testR1index.fastq.gz -iR2 testR2index.fastq.gz -sR1 testR1sequence.fastq.gz -sR2 testR2sequence.fastq.gz -i indexlist.txt -q 35
/usr/bin/time -v ./demultiplex.py -iR1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
-iR2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
-sR1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
-sR2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
-i indexlist.txt -q 35
