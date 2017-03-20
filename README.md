# Transcriptome-

#soybean transcriptome 
#install sratools for fast download SRA data;
wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
tar -xzf sratoolkit.current-centos_linux64.tar.gz

#add the toolkit bin directory into the PATH:
export PATH=/path/to/sratoolkit/bin:$PATH

#want to have a loop to download multiple sra files;
#generate a txt file with all the run-names:
nano download.txt

#then, write sh to download the SRA files:
nano download.sh

#in sh, write the loop:
for i in $(cat download.txt)
do
 fastq-dump -A $i
done

