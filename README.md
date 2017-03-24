### Transcriptome/Soybean transcriptome 
### 0.1 Download soybean reference genome:
```{php}
wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/fasta/glycine_max/dna/Glycine_max.V1.0.28.dna.chromosome.*.fa.gz
```
### 0.2 Soybean transcriptome:
```{php}
wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/fasta/glycine_max/cdna/Glycine_max.V1.0.28.cdna.all.fa.gz
```
### 0.3 Soybean genome annotation:
```{php}
wget ftp://ftp.ensemblgenomes.org/pub/release-34/plants/gtf/glycine_max/Glycine_max.V1.0.34.gtf.gz
```
### 0.4 Install sratools for fast download SRA data;
```{php}
wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
tar -xzf sratoolkit.current-centos_linux64.tar.gz
```
### 0.5 Add the toolkit bin directory into the PATH:
```{php}
export PATH=/path/to/sratoolkit/bin:$PATH
```
## 0.6 Raw data
Want to have a loop to download multiple sra files;
Generate a txt file with all the run-names:
```{php}
nano download.txt
SSR2079326
SSR2079327
SSR2079328
SSR2079329
```
then, write sh to download the SRA files:
in NCBI SRA, the paired reads were joined in one fastq, therefore, flag --split-files would be used to split reads;
```{php}
nano download.sh
for i in $(cat download.txt)
do
 fastq-dump --split-files -A $i
done
```

### 1.0 Quality assessment:
```{php}
qrsh -q medium*
module load fastqc
mkdir 1_fastqc
```
run fastqc:
-o:tells fastqc where we want to put output files;
'./' means 'here' as in the directory you are running the command from. If we don't do this, fastqc will default to put the output files where the input files are, which happens to be in the raw_data folder in this case;
```{php}
fastqc ../raw_data/*.fastq -o ./
```
two output files were generated for each fastq file (*_fastqc.html and *_fastqc.zip);

### 2.0 Trimmomatic
```{php}
mkdir 2_trimmomatic
module load trimmomatic
```
create sh for trimmomatic to run (take 329 for example);
```{php}
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
```
### 3.0 Fastqc
### 4.0 alignment

### 4.1 Toptat
```{php}
module load toptat
module load bowtie2
```
index reference genome first;
soybean has lots of chromosomes and therefore has multiple fa files to index;
while the input file for bowtie2-build can be single fasta file or zip archive of multiple fasta files;
therefore in genome raw folder, zip fa files first;
```{php}
zip genome.zip *.fa
```
build index files:
```{php}
bowtie2-build genome.zip genome
```
errors were shown up with "empty reference" or "reference only gaps"
therefore, combine all the fa files into 1 fa file:
```{php}
cat *.fa > genome.fa
bowtie2-build genome.fa genome
```
seems working fine but when run tophat, tophat cannot find the .bt2 index file;
error was about cannot find .bt2l file (long-index file);
therefore, white a for loop to index these fa files one by one;
```{php}
nano index.sh
for i `seq 1 20`
do 
 bowtie2-build ../0_raw_data/soygenome/gma_ref_Glycine_max_v2.0_chr${i}.fa ../0_raw_data/soygenome/gma_ref_Glycine_max_v2.0_chr${i}
done
```

### 4.2 STAR
Index reference genome
```{php}
mkdir alignment_STAR && cd alignment_STAR
mkdir genomeDir
```

get the readlength (result is 100 in this case, therefore sjdbOverhang set to 100-1=99):
```{php}
head -2 ../0_raw_data/DRR016140_1.1percent.fastq | awk "{print length}" | tail -2
head -2 ../0_raw_data/DRR016140_2.1percent.fastq | awk "{print length}" | tail -2

STAR --runMode genomeGenerate \
    --genomeDir ./genomeDir \
    --genomeFastaFiles ../0_raw_data/soygenome/*.fa \
    --runThreadN 2 \
    --sjdbGTFfile ../0_raw_data/Glycine_max.V1.0.34.gtf \
    --sjdbOverhang 99
```
* --runMode genomeGenerate: run genome indices generation job, default is to run alignment;
* --genomeDir: specify the directory for storing genome indices
* --genomeFastaFiles: one or more FASTA files with genome reference sequences
* --runThreadN: the number of threads to use.
* --sjdbGTFfile: The annotation file that STAR uses to build splice junctions database
* --sjdbOverhang: specifies the length of genomic sequence around the annotated junction. Usually it is set to Readlength - 1.

align the reads:
only 4 pairs of reads files (326-329);
a for loop can be used to get the job done;
in alignment_STAR directory:
```{php}
mkdir alignment_output          ## create a directory to store the alignment output files

for i in `seq 6 9`
do
 STAR --genomeDir ./genomeDir \
      --readFilesIn ../2_trimmomatic/SRR32${i}_1.paired.trimmed.fastq ../2_trimmomatic/SRR32${i}_2.paired.trimmed.fastq \
      --outFileNamePrefix ./alignment_output/SRR32${i}_  \
      --outSAMtype BAM SortedByCoordinate     \
      --runThreadN 2
done
```
* --runMode genomeGenerate: run genome indices generation job, default is to run alignment.
* --genomeDir: specifies the directory where you put your genome indices
* --readFilesIn: your paired RNASeq reads files.
* --outFileNamePrefix: your output file name prefix.
* --outSAMtype: your output file type. Here we want the generated bam file to be sorted by coordination.
* --runThreadN: the number of threads to be used.

#### 4.3 Hisat2
Index reference genome
```{php}
mkdir alignment_hisat2 && cd alignment_hisat2
mkdir genomeDir          ##again, create a directory for the genome indices
```
* -p: specifies the number of threads to use
* ../0_raw_data/*.fa: path to the reference genome
* ./genomeDir/Athal_index: the base of indices files that will be generated.
```{php}
hisat2-build -p 2 ../0_raw_data/soygenome/*.fa ./genomeDir/soybean_index
```
align the reads:
* ./genomeDir/soybean_index: path to the directory of genome indices
* -1: specifies the first paired reads file
* -2: specifies the second paired reads file
* -S: output to a SAM file. Default is stdout.
* -p: the number of threads to use.

```{php}
for i in `seq 6 9`
do
 hisat2 ./genomeDir/soybean_index \
 -1 ../2_trimmomatic/SRR32${i}_1.paired.trimmed.fastq \
 -2 ../2_trimmomatic/SRR32${i}_2.paired.trimmed.fastq \
 -S ./alignment_output/SRR32${i}.sam \
 -p 2
done
```
hisat2 does not generate sorted bam file automatically;
use samtools to convert the sam files to sorted bam files;
```{php}
mkdir sorted_bam

for i in `seq 6 9`
do
 samtools sort ./alignment_output/SRR32${i}.sam -o ./sorted_bam/SRR32${i}_sorted.bam
done
```

### 4.4 Rapmap
```{php}
cd ~/RNASeq
mkdir 4_rapmap && cd 4_rapmap
mkdir transcriptomeDir
```
Index reference genome
```{php}
rapmap quasiindex -t ../0_raw_data/Glycine_max.V1.0.cdna.all.fa -i ./transcriptomeDir
```
* quasiindex: builds a suffix array-based index
* t: specifies the path to the reference transcriptome file
* i: specifies the location where the index should be written
Align the reads
quansimap: map reads using the suffix-array based method, should match the method you used for indexing.
* -1: specifies the first set of reads from a paired library.
* -2: specifies the second set of reads from a paired library.
* -t: number of threads to use.
* -o: path to the file where the output should be written.
```{php}
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
```
submit the job:
```{php}
qsub alignment.qsh
```
generate sorted bam files:
```{php}
mkdir sorted_bam
nano sortedbam.sh
#$ -l mem=6G
#$ -q medium*

for i in `seq 6 8`
do 
 samtools sort ./alignment_output/SRR32${i}.sam -o ./sorted_bam/SRR32${i}_sorted.bam
done
```
### 4.5 Different alignment methods comparison;
so, 3-4 methods have been used to align the reads into genome/cdna 
BAM files comparisons using samtools flagstat tool
now, we can compare the resutls of these alignment results
for each method, in their corresponding folder:
```{php}
mkdir flagstat_output
for i in `seq 6 9`
do
   samtools flagstat ./alignment_output/SRR32${i}_sortedBy.bam > ./flagstat_output/SRR32${i}_flagstat.txt
done
```
check the flag field of each sam file
```{php}
awk '{print $2}' DRR016125.sam | egrep '^[0-9]' | sort -n | uniq

cat alignment_STAR/flagstat_output/SRR326_flagstat.txt
cat alignment_hisat2/flagstat_output/SRR326_flagstat.txt
cat alignment_rapmap/flagstat_output/SRR326_flagstat.txt
```
one output example:
```{php}
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
```

* QC-passed reads: platform/aligner specific flag.
* secondary: the same read may have multiple alignments. One of these alignments is considered primary. Others are secondary.
* supplementary: parts of a reads mapped to different regions (chimeric alignment). One of these aligned regions is representative, others are supplementary.
* duplicates: PCR duplicates. Read the same sequence multiple times.
* with itself and mate mapped: both paired reads are mapped.
* singletons: one of the paired reads is mapped.
* important stats include: properly paired; with itself and mate mapped; secondary                       

### 5.1 Reads counting
after alignment, count the reads
install htseq
```{php}
mkdir 5_reads_count && cd 5_reads_count
conda install htseq -y
conda install pysam -y
```
create a script file count_reads.sh
```{php}
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
```
we did need "grep -v '^__'"(-v: select non-matching lines) in the command line, because:
```{php}
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
```
count reads from STAR alignment:
```{php}
./count_reads.sh STAR /lustre/projects/qcheng1RNA/4_star/alignment_output/ pos
```
count reads from hisat2 alignment:
```{php}
./count_reads.sh hisat2 /lustre/projects/qcheng1RNA/4_hisat2/alignment_output/ name
```
count matrix
```{php}
echo gene_ID $(ls SRR* | grep -o "SRR32[6-8]*" | tr "\n" ' ') | tr -s [:blank:] ',' > count_data.csv
paste $(ls SRR* | sort) | awk '{for(i=3;i<=NF;i+=2) $i=""}{print}' | tr -s [:blank:] ',' >> count_data.csv
```

### 5.2 Expression analysis
copy these count_data.csv in local computer;
in R, analyze the count data with DESeq2 package;
all the following code would be run in R;
### 5.2.1 Install DESq2 in R
```{R}
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
install.packages("colorspace")
library(DESeq2)
```
### 5.2.2 Set working directory and load data into R
```{R}
setwd("C:/Transcriptome")

countData = read.csv('count_data.csv', header = TRUE, row.names = 1)
colData = read.csv("experimental_info.csv", header = TRUE, row.names = 2)[, c("phenotype", "stress")]
colnames(countData) = paste0(colData$phenotype, '_', colData$stress)
rownames(colData) = paste0(colData$phenotype, '_', colData$stress)

## construct the data that analyzing functions from DESeq2 can recognize.
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ phenotype + stress)
dds
```
delete rows that have 0 for all treatments
```{R}
dim(dds)
dds=dds[rowSums(counts(dds))>1,]
dim(dds)
```
### 5.2.3 Differential expression analysis
```{R}
dds = DESeq(dds)
res = results(dds)
res
```
Sample results:
```{R}
log2 fold change (MAP): stress saline vs ABA 
Wald test p-value: stress saline vs ABA 
DataFrame with 19453 rows and 6 columns
             baseMean log2FoldChange     lfcSE        stat       pvalue       padj
            <numeric>      <numeric> <numeric>   <numeric>    <numeric>  <numeric>
AT1G01010   0.9246572     0.04106919 1.1664635  0.03520829   0.97191365         NA
AT1G01020   0.7345512    -0.50228525 1.1997038 -0.41867437   0.67545413         NA
AT1G01030   0.4918562     0.53399011 1.2159565  0.43915231   0.66055118         NA
AT1G01040   4.1881473    -0.88294298 0.7686287 -1.14872497   0.25066941  0.5074783
AT1G01050  12.9203449     0.95553494 0.4842900  1.97306342   0.04848834  0.1769671
```
* baseMean: the average of normalized counts across all samples. This represents the intercept of your GLM.
* #log2Foldchange: (stress saline vs ABA): log2(treated/untreated). Here it is log2(saline/ABA)
* lfcSE: standard error of the log2FoldChange estimate
* stat: statistic for the hypothesis test. For Wald test, it is log2FoldChange/lfcSE.
* pvalue: the corresponding p-value from Wald test (or likelihood ratio test)
* padj: adjusted p value due to multiple comparisons.
* reasons for NA values: By default, independent filtering is performed to select a set of genes for multiple test correction which will optimize the number of adjusted p-values less than a given critical value ‘alpha’ (by default 0.1). The adjusted p-values for the genes which do not pass the filter threshold are set to ‘NA’. By default, the mean of normalized counts is used to perform this filtering, though other statistics can be provided. Several arguments from the ‘filtered_p’ function of genefilter are provided here to control or turn off the independent filtering behavior.
* By default, ‘results’ assigns a p-value of ‘NA’ to genes containing count outliers, as identified using Cook's distance. See the ‘cooksCutoff’ argument for control of this behavior. Cook's distances for each sample are accessible as a matrix "cooks" stored in the ‘assays()’ list. This measure is useful for identifying rows where the observed counts might not fit to a Negative Binomial distribution.

if we do not want NA, these two arguments can be added
```{R}
results(dds,independentFiltering=FALSE, cooksCutoff=Inf)
```
Sample results
```{R}
log2 fold change (MAP): stress saline vs ABA 
Wald test p-value: stress saline vs ABA 
DataFrame with 19453 rows and 6 columns
             baseMean log2FoldChange     lfcSE        stat       pvalue        padj
            <numeric>      <numeric> <numeric>   <numeric>    <numeric>   <numeric>
AT1G01010   0.9246572     0.03501531 1.1418071  0.03066658   0.97553545   1.0000000
AT1G01020   0.7345512    -0.45421205 1.1660316 -0.38953665   0.69687920   1.0000000
AT1G01030   0.4918562     0.45730016 1.1722862  0.39009260   0.69646808   1.0000000
AT1G01040   4.1881473    -0.86833117 0.7735613 -1.12251113   0.26164518   0.7854604
AT1G01050  12.9203449     0.94929113 0.4899696  1.93744912   0.05269047   0.3171373
```
these two results can be compared. Higher pvalue in result1 would be "filtered out", therefore, padj would be NA;

* By default, the results() function display comparison of the last level of the last variable over the first level of this variable.
We can extract results of comparisons between other levels:
```{R}
results(dds, contrast=c("stress", "mock", "saline"))
results(dds, contrast=c("phenotype", "ros1-3", "wildtype"))MA
```
MA-plot of the results
```{R}
plotMA(res)
```

![plotma](https://cloud.githubusercontent.com/assets/23083944/24308124/36cc327a-109d-11e7-873f-799267d72515.jpeg)

We can select the gene which has the smallest padj value to plot its count at different levels
```{R}
mostSigGene = rownames(dds)[which.min(res$padj)]
mostSigGene
par(mfcol=c(1,2))
plotCounts(dds, gene=mostSigGene, intgroup="stress")
plotCounts(dds, gene=mostSigGene, intgroup="phenotype")
```
export results (padj less than 0.05) into csv files
```{R}
orderedRes = res[order(res$padj), ]
sigOrderedRes = subset(orderedRes, padj < 0.05)
write.csv(as.data.frame(sigOrderedRes), "STRESS_saline_vs_ABA.csv")
```
### 5.2.4 Likelihood Ratio Test (LRT) method
The LRT method compares the full model with the reduced model to see if the removed variable (factor) has significant effect on the fitted model
```{R}
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ phenotype + stress)
dds = DESeq(dds, test="LRT", reduced = ~phenotype)
resLRT = results(dds)
resLRT
```
### 5.2.5 Wald test vs Likelihood Ratio test
How many genes are significant in LRT test?
```{R}
orderedResLRT = resLRT[order(resLRT$padj), ]
sigOrderedResLRT = subset(orderedResLRT, padj<0.05)
dim(sigOrderedResLRT)
```
How many genes are significant in Wald test?
```{R}
orderedRes = res[order(res$padj), ]
sigOrderedRes = subset(orderedRes, padj < 0.05)
dim(sigOrderedRes)
```
__Why is the number from LRT much larger than the number from Wald test?__

Let's pick up a gene that is significant in LRT but not in Wald test.

```{R}
LRTsigGenes = rownames(sigOrderedResLRT)
Res_genesSigInLRT = res[LRTsigGenes, ]
notSigWaldGene = rownames(Res_genesSigInLRT)[which.max(Res_genesSigInLRT$padj)]
notSigWaldGene
```

padj comparison for the selected gene
```{R}
resLRT[notSigWaldGene, ]
res[notSigWaldGene, ]
```

Let's plot the counts for this gene

```{R}
plotCounts(dds, gene=notSigWaldGene, intgroup="stress")
```

![at4g39090](https://cloud.githubusercontent.com/assets/23083944/24307997/a61ff69e-109c-11e7-9b51-180e22876c80.jpeg)

Let's extract the results of comparisons between __dehydration__ and "ABA".

```{R}
res2 = results(dds, contrast=c("stress", "dehydration", "ABA"))
res2[notSigWaldGene, ]
```






