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


#1_quality assessment:

qrsh -q medium*
module load fastqc
mkdir 1_fastqc

#run fastqc:
#-o:tells fastqc where we want to put output files;
#'./' means 'here' as in the directory you are running the command from. If we don't do this, fastqc will default to put the output files where the input files are, which happens to be in the raw_data folder in this case;

fastqc ../raw_data/*.fastq -o ./

#two output files were generated for each fastq file (*_fastqc.html and *_fastqc.zip);


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








