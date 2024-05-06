# made a directory for all of the final project materials to go in
mkdir final_project/raw_data

# copied the raw data into the new directory
cp -r /tmp/gen711_project_data/genome-assembly-fqs/. raw_data

# moved just the reads we wanted into a separate directory
mkdir our_reads
cp ./raw_data/36_S3_L001_R1_001.fastq.gz ./raw_data/36_S3_L001_R2_001.fastq.gz ./our_reads

# ran fastqc to analyze read quality, then got the html files onto our personal devices
mkdir fastqc_output
tmux new -s fastqc_finalproject 
fastqc ./our_reads/*fastq.gz -o fastqc_output
# (In local terminal)
sftp oms1034@ron.sr.unh.edu
lcd Desktop
cd fastqc_output
get *.html 

# ran trimmomatic to get rid of adapter sequences since those were most prevalent in our fastqc analysis
mkdir reads_trimmed
cp /tmp/trim_scriptV2.sh .
conda create -n trim 
tmux new -s trim
./trim_scriptV2.sh ../our_reads/36_S3_L001_R1_001.fastq.gz ../our_reads/36_S3_L001_R2_001.fastq.gz

