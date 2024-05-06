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

# ran fastqc on the trimmed reads to see how trimming improved quality, brought html files to our personal devices again
mkdir trimmed_fastqc
tmux new -s fastqc_trimmed
fastqc ./trimmed-reads/*fastq.gz -o trimmed_fastqc
# (In local terminal)
sftp oms1034@ron.sr.unh.edu
lcd Desktop
cd final_project/reads_trimmed/trimmed_fastqc
get *.html 

# ran SPAdes to assemble the trimmed reads and the unpaired reads that trimmomatic gave us
tmux attach -t trim
conda activate genomics
spades.py -1 ../reads_trimmed/trimmed-reads/36_S3_L001_R1_001.fastq.gz -2 ../reads_trimmed/trimmed-reads/36_S3_L001_R2_001.fastq.gz -s ../reads_trimmed/trimmed-reads/unpaired-36_S3_L001_R1_001.fastq.gz -s ../reads_trimmed/trimmed-reads/unpaired-36_S3_L001_R2_001.fastq.gz -o spades_assembly_default -t 24
cp contigs.fasta spades.log ..

# ran QUAST to look for contiguity in our assembled genome
quast.py ./contigs.fasta -o ../quast_results

# ran BUSCO to examine the completeness of our assembled genome using the bacteria dataset available on RON
conda activate busco 
busco -i spades_assembly/contigs.fasta -m genome -o busco_results -l bacteria
less -S 

# ran prokka to annotate our genome, the .gff file that is an output of PROKKA is used later for visualization
prokka spades_assembly/contigs.fasta --outdir prokka_output --cpus 24 --mincontiglen 200
grep -o "product=.*" ./prokka_output/PROKKA_*.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt
mv protein_abundances.txt prokka_output





