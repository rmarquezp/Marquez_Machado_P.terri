## Run MaSurCa (config file is on repo)

masurca Masurca_config.txt
bash ./assemble.sh


## Convert PacBio files to Fastq. Since I'm not using Illumina dta to scaffold I took all reads longer than 100. 
 

for i in *.bam; do bamtools filter -length ">100" -in $i | bamtools convert -format fastq -out $i.over100.fastq; gzip $i.over100.fastq; done

## RNA-seq based scaffolding. Downloaded SRX8741407 from SRA for this, files are called PbicRNAseq....


#Read cleanup
java -jar -Xmx16g ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 24 -phred33 PbicRNAseq_R1_001.fastq.gz PbicRNAseq_R2_001.fastq.gz PbicRNAseq_trim_R1.fastq.gz PbicRNAseq_trim_R1U.fastq.gz PbicRNAseq_trim_R2.fastq.gz PbicRNAseq_trim_R2U.fastq.gz ILLUMINACLIP:/home/phyllobates/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:20:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:36
~/software/skewer/skewer -t 32 -m pe -l 36 -n -Q 20 -z -o PbicRNAseq PbicRNAseq_trim_R1.fastq.gz  PbicRNAseq_trim_R2.fastq.gz

## Count to make sure equal numbers of reads

for i in *pair*.gz; do echo $i; echo $(zcat $i|wc -l)/4|bc; done

## Map

~/software/bwa/bwa index pterr-scaffold-masurca.fasta 
~/software/bwa/bwa mem -t 72 ../Massurca_IlluminaOnly/pterr-scaffold-masurca.fasta PbicRNAseq-trimmed-pair1.fastq.gz PbicRNAseq-trimmed-pair2.fastq.gz > PbicRNAseq.sam

## P_RNA_Scaffolder

bash ~/software/P_RNA_scaffolder/P_RNA_scaffolder.sh -d /home/phyllobates/software/P_RNA_scaffolder/ -i PbicRNAseq.sam -j ../Massurca_IlluminaOnly/pterr-scaffold-masurca.fasta -F PbicRNAseq-trimmed-pair1.fastq.gz -R PbicRNAseq-trimmed-pair2.fastq.gz -o ./P_RNA_scaff -s yes -b yes -p 0.95 -t 48 -e 500000 -f 4

## Nextpolish (script is on repo)

bash ~/software/NextPolish/RunNextpolish.sh P_RNA_scaffold.pilon.fasta /home/phyllobates/reads/Illumina/PE/AllPE_R1.fastq.gz /home/phyllobates/reads/Illumina/PE/AllPE_R2.fastq.gz 1 120

## Scaffolding/Gapfilling Pipeline -  run four times, plus a final NextPolish at the end

input=~/assemblies/RNAseqScaff/P_RNA_scaff/NextPolish/genome.nextpolish.fa
PacBio=~/reads/PB/fasta/readCorrection/P.terribilisPacBio_bamtools.over100.fmlrc_corrected.fasta
PE1=~/reads/Illumina/PE/AllPE_R1.fastq.gz
PE2=~/reads/Illumina/PE/AllPE_R2.fastq.gz
HiC1=~/reads/HiC/Pterr_HiC-trimmed-pair1.fastq.gz
HiC2=~/reads/HiC/Pterr_HiC-trimmed-pair2.fastq.gz
round=4

for ((i=1; i<=${round};i++)); do

	mkdir Round"$i"
	cd Round"$i"
	
	echo "Running pipeline on $input"
	
	bash ~/software/runSALSA.sh $input $HiC1 $HiC2 120
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $PacBio > PacBio_singleLine.fa
	bash ~/software/LR_Gapcloser/src/LR_Gapcloser.sh -i SALSA_HiCscaff/scaffolds_FINAL.fasta -l PacBio_singleLine.fa -s p -t 120 -o LR -r 5 -v 150
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' LR/iteration-5/gapclosed.fasta > LR_assembly_SingleLine.fa
	bash ~/software/RAILS_v1.5.1/runRAILSminimapSTREAM.sh LR_assembly_SingleLine.fa PacBio_singleLine.fa 100 0.90 100 1 pacbio ~/software/samtools-1.11
	mv *_rails.scaffolds.fa Rails.final.fa
	#ragtag.py scaffold Rails.final.fa Rails.final.fa -r -t 120 -u
	bash ~/software/NextPolish/RunNextpolish.sh Rails.final.fa $PE1 $PE2 1 120
	
	input=../Round"$i"/genome.nextpolish.FINAL.fa
	cd ..
	
done

cp Round"$round"/genome.nextpolish.FINAL.fa genome.round"$round".fa
bash ~/software/NextPolish/RunNextpolish.sh genome.round"$round".fa $PE1 $PE2 3 120

## This resulted in the final assembly