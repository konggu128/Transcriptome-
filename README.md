# Transcriptome-

#soybean transcriptome 
#download soybean reference genome:
wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/fasta/glycine_max/dna/Glycine_max.V1.0.28.dna.chromosome.*.fa.gz

#soybean transcriptome:
wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/fasta/glycine_max/cdna/Glycine_max.V1.0.28.cdna.all.fa.gz

#soybean genome annotation:
wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/gtf/glycine_max/Glycine_max.V1.0.34.gtf.gz

#install sratools for fast download SRA data;
wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
tar -xzf sratoolkit.current-centos_linux64.tar.gz

#add the toolkit bin directory into the PATH:

export PATH=/path/to/sratoolkit/bin:$PATH

#want to have a loop to download multiple sra files;
#generate a txt file with all the run-names:

nano download.txt
SSR2079326
SSR2079327
SSR2079328
SSR2079329

#then, write sh to download the SRA files:
#in NCBI SRA, the paired reads were joined in one fastq, therefore, flag --split-files would be used to split reads;

nano download.sh

#in sh, write the loop:

for i in $(cat download.txt)
do
 fastq-dump --split-files -A $i
done

#####################################################################################################################################
#1_quality assessment:

qrsh -q medium*
module load fastqc
mkdir 1_fastqc

#run fastqc:
#-o:tells fastqc where we want to put output files;
#'./' means 'here' as in the directory you are running the command from. If we don't do this, fastqc will default to put the output files where the input files are, which happens to be in the raw_data folder in this case;

fastqc ../raw_data/*.fastq -o ./

#two output files were generated for each fastq file (*_fastqc.html and *_fastqc.zip);


#####################################################################################################################################
#2_trimmomatic
mkdir 2_trimmomatic
module load trimmomatic

#create sh for trimmomatic to run (take 329 for example);
nano trimmo.sh

trimmomatic \
PE \
-trimlog soyheat329 \
-phred33 \
../0_raw_data/SRR2079329_1.fastq \
../0_raw_data/SRR2079329_2.fastq \
SRR329_1.paired.trimmed.fastq \
SRR329_1.unpaired.trimmed.fastq \
SRR329_2.paired.trimmed.fastq \
SRR329_2.unpaired.trimmed.fastq \
SLIDINGWINDOW:4:15 \
MINLEN:30 \
>& 329.output

#####################################################################################################################################
#3_fastqc
#4_alignment

############################
########1.toptat############
############################
module load toptat
module load bowtie2

#index reference genome first;
#soybean has lots of chromosomes and therefore has multiple fa files to index;
#while the input file for bowtie2-build can be single fasta file or zip archive of multiple fasta files;
#therefore in genome raw folder, zip fa files first;
zip genome.zip *.fa

#build index files:
bowtie2-build genome.zip genome

#errors were shown up with "empty reference" or "reference only gaps"
#therefore, combine all the fa files into 1 fa file:
cat *.fa > genome.fa
bowtie2-build genome.fa genome

#seems working fine but when run tophat, tophat cannot find the .bt2 index file;
#error was about cannot find .bt2l file (long-index file);
#therefore, white a for loop to index these fa files one by one;
nano index.sh
for i `seq 1 20`
do 
 bowtie2-build ../0_raw_data/soygenome/gma_ref_Glycine_max_v2.0_chr${i}.fa ../0_raw_data/soygenome/gma_ref_Glycine_max_v2.0_chr${i}
done

############################
########2.STAR############
############################
#Index reference genome
mkdir alignment_STAR && cd alignment_STAR
mkdir genomeDir

##--runMode genomeGenerate: run genome indices generation job, default is to run alignment;
##--genomeDir: specify the directory for storing genome indices
##--genomeFastaFiles: one or more FASTA files with genome reference sequences
##--runThreadN: the number of threads to use.
##--sjdbGTFfile: The annotation file that STAR uses to build splice junctions database
##--sjdbOverhang: specifies the length of genomic sequence around the annotated junction. Usually it is set to Readlength - 1.
#get the readlength (result is 100 in this case, therefore sjdbOverhang set to 100-1=99):
head -2 ../0_raw_data/DRR016140_1.1percent.fastq | awk "{print length}" | tail -2
head -2 ../0_raw_data/DRR016140_2.1percent.fastq | awk "{print length}" | tail -2

STAR --runMode genomeGenerate \
    --genomeDir ./genomeDir \
    --genomeFastaFiles ../0_raw_data/soygenome/*.fa \
    --runThreadN 2 \
    --sjdbGTFfile ../0_raw_data/Glycine_max.V1.0.34.gtf \
    --sjdbOverhang 99

#align the reads:
#only 4 pairs of reads files (326-329);
#a for loop can be used to get the job done;
#in alignment_STAR directory:

mkdir alignment_output          ## create a directory to store the alignment output files

##--runMode genomeGenerate: run genome indices generation job, default is to run alignment.
##--genomeDir: specifies the directory where you put your genome indices
##--readFilesIn: your paired RNASeq reads files.
##--outFileNamePrefix: your output file name prefix.
##--outSAMtype: your output file type. Here we want the generated bam file to be sorted by coordination.
##--runThreadN: the number of threads to be used.

for i in `seq 6 9`
do
 STAR --genomeDir ./genomeDir \
      --readFilesIn ../2_trimmomatic/SRR32${i}_1.paired.trimmed.fastq ../2_trimmomatic/SRR32${i}_2.paired.trimmed.fastq \
      --outFileNamePrefix ./alignment_output/SRR32${i}_  \
      --outSAMtype BAM SortedByCoordinate     \
      --runThreadN 2
done

############################
########3.hisat2############
############################
#Index reference genome
mkdir alignment_hisat2 && cd alignment_hisat2
mkdir genomeDir          ##again, create a directory for the genome indices

##-p: specifies the number of threads to use
##../0_raw_data/*.fa: path to the reference genome
##./genomeDir/Athal_index: the base of indices files that will be generated.

hisat2-build -p 2 ../0_raw_data/soygenome/*.fa ./genomeDir/soybean_index

##align the reads:
##./genomeDir/soybean_index: path to the directory of genome indices
##-1: specifies the first paired reads file
##-2: specifies the second paired reads file
##-S: output to a SAM file. Default is stdout.
##-p: the number of threads to use.

for i in `seq 6 9`
do
 hisat2 ./genomeDir/soybean_index \
 -1 ../2_trimmomatic/SRR32${i}_1.paired.trimmed.fastq \
 -2 ../2_trimmomatic/SRR32${i}_2.paired.trimmed.fastq \
 -S ./alignment_output/SRR32${i}.sam \
 -p 2
done

#hisat2 does not generate sorted bam file automatically;
#use samtools to convert the sam files to sorted bam files;

mkdir sorted_bam

for i in `seq 6 9`
do
 samtools sort ./alignment_output/SRR32${i}.sam -o ./sorted_bam/SRR32${i}_sorted.bam
done


############################
########4.rapmap############
############################
cd ~/RNASeq
mkdir 4_rapmap && cd 4_rapmap
mkdir transcriptomeDir

#Index reference genome
##quasiindex: builds a suffix array-based index
##t: specifies the path to the reference transcriptome file
##i: specifies the location where the index should be written

rapmap quasiindex -t ../0_raw_data/Glycine_max.V1.0.cdna.all.fa -i ./transcriptomeDir

#Align the reads
##quansimap: map reads using the suffix-array based method, should match the method you used for indexing.
##-1: specifies the first set of reads from a paired library.
##-2: specifies the second set of reads from a paired library.
##-t: number of threads to use.
##-o: path to the file where the output should be written.

mkdir alignment_output
nano alignment.qsh
#$ -l mem=7G
#$ -q medium*

for i in `seq 6 8`
do
 rapmap quasimap \
 -i ./transcriptomeDir \
 -1 ../0_raw_data/SRR32${i}_1.paired.trimmed.fastq \
 -2 ../0_raw_data/SRR32${i}_2.paired.trimmed.fastq \
 -t 2 \
 -o ./alignment_output/SRR32${i}.sam    
done

#submit the job:
qsub alignment.qsh

#generate sorted bam files
mkdir sorted_bam
nano sortedbam.sh
#$ -l mem=6G
#$ -q medium*

for i in `seq 6 8`
do 
 samtools sort ./alignment_output/SRR32${i}.sam -o ./sorted_bam/SRR32${i}_sorted.bam
done

#####################################################################################################################
#so, 3-4 methods have been used to align the reads into genome/cdna; 
#BAM files comparisons using samtools flagstat tool;
#now, we can compare the resutls of these alignment results;
#for each method, in their corresponding folder:

mkdir flagstat_output
for i in `seq 6 9`
do
   samtools flagstat ./alignment_output/SRR32${i}_sortedBy.bam > ./flagstat_output/SRR32${i}_flagstat.txt
done


#check the flag field of each sam file;

awk '{print $2}' DRR016125.sam | egrep '^[0-9]' | sort -n | uniq


cat alignment_STAR/flagstat_output/SRR326_flagstat.txt
cat alignment_hisat2/flagstat_output/SRR326_flagstat.txt
cat alignment_rapmap/flagstat_output/SRR326_flagstat.txt

#one output example:
[Newton:sigma00 RNASeq_lab_I]$ cat alignment_rapmap/flagstat_output/SRR326_flagstat.txt
335196 + 0 in total (QC-passed reads + QC-failed reads)
109792 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
331260 + 0 mapped (98.83% : N/A)
225404 + 0 paired in sequencing
112702 + 0 read1
112702 + 0 read2
220072 + 0 properly paired (97.63% : N/A)
220072 + 0 with itself and mate mapped
2666 + 0 singletons (1.18% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

##QC-passed reads: platform/aligner specific flag.
##secondary: the same read may have multiple alignments. One of these alignments is considered primary. Others are secondary.
##supplementary: parts of a reads mapped to different regions (chimeric alignment). One of these aligned regions is representative, others are supplementary.
##duplicates: PCR duplicates. Read the same sequence multiple times.
##with itself and mate mapped: both paired reads are mapped.
##singletons: one of the paired reads is mapped.

##########################################################
#important stats include:         ########################
#properly paired                  ########################
#with itself and mate mapped      ########################
#secondary                        ########################
##########################################################



####################################################################################################################
#after alignment, count the reads;
#install htseq
mkdir 5_reads_count && cd 5_reads_count
conda install htseq -y
conda install pysam -y

#create a script file count_reads.sh

#!/bin/bash

## USAGE
## ./count_reads.sh hisat /path/to/the/sorted_bam_file_directory bam_file_ORDER_TYPE
## ./count_reads.sh STAR /path/to/the/sorted_bam_file_directory bam_file_ORDER_TYPE

mkdir counts_$1

for sorted_bam_path in $(find $2 -name *.bam)
do
    counts_file=$(echo $sorted_bam_path | grep -o "SRR32[6-8]*")_$1_ct
    echo "The target bam file is: "$sorted_bam_path
    echo "==================================================="
    htseq-count -f bam \
                -t gene \
                -i gene_id  \
                -r $3 \
                $sorted_bam_path \
                /lustre/projects/qcheng1RNA/0_raw_data/ Glycine_max.V1.0.34.gtf \
                | \
                grep -v '^__' > ./counts_$1/$counts_file
    echo "The count data has been written into: $counts_file"
    echo "==================================================="
done


#we did need "grep -v '^__'"(-v: select non-matching lines) in the command line, because:
...
...
...
ATMG01370	0
ATMG01380	0
ATMG01390	272
ATMG01400	0
ATMG01410	0
__no_feature	8798
__ambiguous	4637
__too_low_aQual	1423
__not_aligned	2794
__alignment_not_unique	52464

#count reads from STAR alignment:
./count_reads.sh STAR /lustre/projects/4_star/alignment_output/ pos

#count reads from hisat2 alignment:
./count_reads.sh hisat2 /lustre/projects/4_hisat2/alignment_output/ name

#count matrix
echo gene_ID $(ls SRR* | grep -o "SRR32[6-8]*" | tr "\n" ' ') | tr -s [:blank:] ',' > count_data.csv
paste $(ls SRR* | sort) | awk '{for(i=3;i<=NF;i+=2) $i=""}{print}' | tr -s [:blank:] ',' >> count_data.csv







