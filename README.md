# Transcriptome-

#soybean transcriptome 
#install sratools for fast download SRA data;
wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
tar -xzf sratoolkit.current-centos_linux64.tar.gz

#add the toolkit bin directory into the PATH:
export PATH=/path/to/sratoolkit/bin:$PATH

#then, download the SRA files:
