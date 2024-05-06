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

# we extracted the 16s rRNA information from the PROKKA .ffn file to use in our BLAST analysis
grep 16S prokka_output/*.ffn
extract_sequences "16S ribosomal RNA" prokka_output/PROKKA_*.ffn > organism_identification/16S_sequence.fasta

# we used BLAST on the NCBI website to look for species that match our 16s rRNA 
# https://blast.ncbi.nlm.nih.gov/Blast.cgi
# to do this, we copied and pasted the 16s rRNA sequence into BLAST and waited for the program to match with other species in the database
nano 16S_sequence.fasta

# we did read mapping using samtools to map the forward and reverse reads onto our completed reference genome
mkdir read_mapping
cp ../spades_assembly/contigs.fasta .
bwa index contigs.fasta
bwa mem -t 24 contigs.fasta ../reads_trimmed/trimmed-reads/36_S3_L001_R1_001.fastq.gz ../reads_trimmed/trimmed-reads/36_S3_L001_R2_001.fastq.gz > raw_mapped.sam
samtools view -@ 24 -Sb  raw_mapped.sam  | samtools sort -@ 24 - sorted_mapped
samtools flagstat sorted_mapped.bam
samtools index sorted_mapped.bam
bedtools genomecov -ibam sorted_mapped.bam > coverage.out
gen_input_table.py  --isbedfiles contigs.fasta coverage.out >  coverage_table.tsv

# ran BLAST again using command line because we realized we were missing a file to use in the non-target contig removal step
mkdir blast
blob_blast.sh ../spades_assembly/contigs.fasta 

# we completed non-target contig removal to prepare our data for visualization
# visualization was completed using blobtools and our ooutput files were brought to our local devices
mkdir contig_removal
blobtools create -i ../spades_assembly/contigs.fasta -b ../read_mapping/sorted_mapped.bam -t ../blast/contigs.fasta.vs.nt.cul5.1e5.megablast.out -o blob_out
blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy
grep -v '##' blob_taxonomy.blob_out.blobDB.table.txt
blobtools plot -i blob_out.blobDB.json -r genus
# (In local terminal)
cd contig_removal
lcd Desktop
get *.png

# we filtered our outputs by both length and coverage
mkdir filtered_assembly
cp contig_removal/blob_taxonomy.blob_out.blobDB.table.txt .
mv blob_taxonomy.blob_out.blobDB.table.txt ./filtered_assembly/
grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2 > 500' | wc -w
grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2 < 500' | wc -w
grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 < 5' | wc -w
grep -v '#' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 5' | wc -w

# we constructed a list of contigs to keep from filtering and put them into a file called "list_of_contigs_to_keep_len500_cov20.txt"
grep -v '##' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | less
grep -v '##' blob_taxonomy.blob_out.blobDB.table.txt | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | awk -F'\t' '{print $1}' > list_of_contigs_to_keep_len500_cov20.txt

# we filtered the genome assembly by the list of contigs that we wanted to keep from filtering
# called our final, filtered assembly "Microbacterium.fasta"
filter_contigs_by_list.py spades_assembly/contigs.fasta filtered_assembly/list_of_contigs_to_keep_len500_cov20.txt Microbacterium.fasta

# we looked at the coverage across our contigs from the blobtools table that we created
grep -f filtered_assembly/list_of_contigs_to_keep_len500_cov20.txt filtered_assembly/blob_taxonomy.blob_out.blobDB.table.txt | awk '{w = w + $2; e = e + $5 * $2;} END {print e/w}'

# we ran BLAST again on our final alignment against UniVec to make sure there was no contamination
wget "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"
blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query ../Microbacterium.fasta -subject UniVec  -outfmt 6 -out genome_vs_univec.6

# ran QUAST and BUSCO again to see how our original assembly compared to our filtered assembly
quast.py Microbacterium.fasta -o quast_again
conda activate busco 
busco -i Microbacterium.fasta -m genome -o busco_again -l bacteria

# indexed the Microbacterium.fasta file using samtools to make it suitable for IGV Web Viewer visualization
samtools faidx ../Microbacterium.fasta

