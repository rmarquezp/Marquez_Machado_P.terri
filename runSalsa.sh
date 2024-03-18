#!/bin/bash

##Usage: runSALSA.sh inputAssembly.fa HiC_R1.fastq.gz HIC_R2.fastq.gz threads

assembly=$1
R1=$2
R2=$3
threads=$4

bwa index -a bwtsw $assembly
samtools faidx $assembly

# Map R1 and R2 separately (Also tried with the uncleaned reads to see what happens)

bwa mem -t $threads $assembly $R1 | samtools view -@ 24 -Sb - > SALSA_HiC_UF_1.bam
bwa mem -t $threads $assembly $R2 | samtools view -@ 24 -Sb - > SALSA_HiC_UF_2.bam

#Run the Arima custom filters

#Find reads that have pieces of the 5' and 3' sides of the ligation junction and keep only the piece in the 5' end
samtools view -h SALSA_HiC_UF_1.bam | perl ~/software/mapping_pipeline/filter_five_end.pl | samtools view -Sb - > SALSA_HiC_5F_1.bam&
samtools view -h SALSA_HiC_UF_2.bam | perl ~/software/mapping_pipeline/filter_five_end.pl | samtools view -Sb - > SALSA_HiC_5F_2.bam&
wait

### Join filtered files
perl ~/software/mapping_pipeline/two_read_bam_combiner.pl SALSA_HiC_5F_1.bam SALSA_HiC_5F_2.bam ~/software/samtools-1.11/samtools 10 | samtools view -bS -t "$assembly".fai - | samtools sort -@ 120 -m 8gb -o SALSA_HiC_joint.bam -

## Add read groups

java -Xmx480G -jar ~/software/picard.jar AddOrReplaceReadGroups INPUT=SALSA_HiC_joint.bam OUTPUT=SALSA_joint.rg.bam ID=SALSA_HiC LB=Pter SM=SALSAHiC PL=ILLUMINA PU=none

## RMdup

java -Xmx480G -XX:-UseGCOverheadLimit -jar -jar ~/software/picard.jar MarkDuplicates INPUT=SALSA_joint.rg.bam OUTPUT=SALSA_joint.rg.dedup.bam METRICS_FILE=HiCDedup.txt ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
samtools index SALSA_joint.rg.dedup.bam

#Stats

perl ~/software/mapping_pipeline/get_stats.pl SALSA_joint.rg.dedup.bam > HiC_stats.txt #NOTE: This needs samtools in the path to work!!

## clean up

rm SALSA_HiC_UF_?.bam SALSA_HiC_5F_?.bam SALSA_HiC_joint.bam SALSA_joint.rg.bam

## Now SALSA (needs fine tuning for sure)

## Make bam into bed

~/software/bedtools bamtobed -i SALSA_joint.rg.dedup.bam > SALSA_joint.rg.dedup.bed
sort -k 4 SALSA_joint.rg.dedup.bed > tmp && mv tmp SALSA_joint.rg.dedup.bed

## Run pipeline

python ~/software/SALSA/run_pipeline.py -a $assembly -l "$assembly".fai -b SALSA_joint.rg.dedup.bed -e GATC -o SALSA_HiCscaff -c 500 -m yes -i 30
