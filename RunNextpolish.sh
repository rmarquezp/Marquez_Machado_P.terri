#!/bin/bash

## Usage: RunNextpolish.sh input.fa Illumina_R1.fastq.gz Illumina_R2.fastq.gz nRounds nThreads

#Set input and parameters
round=$4
threads=$5
readPE1=$2
readPE2=$3
input=$1

for ((i=1; i<=${round};i++)); do

	#step 1:
	bwa index ${input};
	bwa mem -t ${threads} ${input} ${readPE1} ${readPE2}|samtools view --threads 3 -F 0x4 -b -|samtools fixmate -m --threads 3  - -|samtools sort -m 24g --threads 12 -|samtools markdup --threads 5 -r - sgs.PE.sort.mkdup.bam
	samtools index -@ ${threads} sgs.PE.sort.mkdup.bam;
	samtools faidx ${input};

	#polish genome file
	python ~/software/NextPolish/lib/nextpolish1.py -g ${input} -t 1 -p ${threads} -s sgs.PE.sort.mkdup.bam > genome.polishtemp.fa;
	input=genome.polishtemp.fa;

	#step2:
	bwa index ${input};
	bwa mem -t ${threads} ${input} ${readPE1} ${readPE2}|samtools view --threads 3 -F 0x4 -b -|samtools fixmate -m --threads 3  - -|samtools sort -m 24g --threads 12 -|samtools markdup --threads 5 -r - sgs.PE.sort.mkdup.bam
	samtools index -@ ${threads} sgs.PE.sort.mkdup.bam;
	samtools faidx ${input};

	#polish genome file
	python ~/software/NextPolish/lib/nextpolish1.py -g ${input} -t 2 -p ${threads} -s sgs.PE.sort.mkdup.bam > genome.nextpolish.fa;
	input=genome.nextpolish.fa;
done;

echo "\n\n----DONE----\n"
