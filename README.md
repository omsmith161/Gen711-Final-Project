# Gen711-Final-Project

## Background
The Original Team's Data: 
- Samples of seaweed eating microbes were collected from Acadia National Park
https://docs.google.com/presentation/d/1OBzO8tTlOovftlic2hYZ7EWInaVoaPCbp_dOabNGVWY/edit#slide=id.g13fdf044cda_11_289
 
Steps Used to Obtain Samples:
1. Isolating bacteria
2. DNA Extraction
3. Prepare Library
4. Sequence using Illumina HiSeq 2500
- there were forward & reverse reads for each sample in FASTQ format

We used data that was found in a shared directory on ron.
-- The Goal: To have a complete genome assembly and assessment of the chosen reads from the given data.

-the chosen forward & reverse reads:

36_S3_L001_R1_001.fastq.gz
36_S3_L001_R2_001.fastq.gz

## Methods
We completed this genome assembly and analysis using a variety of packages on RON. Our pipeline is outlined in the figure below. 
[Pipeline.pdf](https://github.com/omsmith161/Gen711-Final-Project/files/15194420/Pipeline.pdf)
Figure 1. Bioinformatics pipeline that we used for our genome assembly and analysis. 

## Results 
Running trimmomatic successfully allowed us to trim the adapter sequences off of our forward and reverse reads. This was shown in our fastqc analysis. The most overrepresented sequence in the forward read was the TruSeq Adapter, also called the Illumina Universal Adapter. This sequence was not considered overrepresented in the reverse read, however, the Illumina Universal Adapter was also present in the "Adapter Content" section of the fastqc output. In both the forward and the reverse read, 0 sequences were flagged for poor quality, so trimmomatic was not needed for that purpose. After using trimmomatic, our additional fastqc analysis showed that it successfully removed adapter sequences since it was no longer flagging anything overrepresented. Additionally, the Illumina Universal Adapter became a lot less prevalent in the "Adapter Content" section of the output. 
<img width="972" alt="Screen Shot 2024-05-02 at 6 08 38 PM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/116675192/24aacbc9-b880-4e76-bd91-4289944b65ec">
<img width="936" alt="Screen Shot 2024-05-02 at 6 12 01 PM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/116675192/cc178afb-3ac0-4214-b924-a306ccb528b5">
Figure 2. Comparing the forward read Adapter Content before and after trimming adapters using trimmomatic. 

## Pre-Filtered Quast Data: 
Assembly------------------------contigs
- Number of Contigs (>= 0 bp)      80
- Number of Contigs (>= 1000 bp)   25
- Number of Contigs (>= 5000 bp)   22
- Number of Contigs (>= 10000 bp)  19
- Number of Contigs (>= 25000 bp)  17
- Number of Contigs (>= 50000 bp)  14
- Total length (>= 0 bp)           3112414
- Number of Contigs                28

## Pre-Filtered Busco Data:
***** Results: *****

	C:99.2%[S:99.2%,D:0.0%],F:0.0%,M:0.8%,n:124	   
	123	Complete BUSCOs (C)			   
	123	Complete and single-copy BUSCOs (S)	   
	0	Complete and duplicated BUSCOs (D)	   
	0	Fragmented BUSCOs (F)			   
	1	Missing BUSCOs (M)			   
	124	Total BUSCO groups searched

## Prokka
One of our outputs was a file called protein_abundances.txt, and this gave us a list of all the hypothetical proteins that could possibly come from the genes in this genome. The most important gene that we recieved from the analysis was a 16s rRNA gene. This gene was used and inserted into the BLAST database. 

## BLAST Results
<img width="695" alt="Screen Shot 2024-05-13 at 2 05 45 PM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/158241303/50c46d6b-0c94-49df-b080-c14905fd198f">

## Blobtools 

![blob_out blobDB json bestsum genus p8 span 100 blobplot bam0](https://github.com/omsmith161/Gen711-Final-Project/assets/158241303/f20add11-23fc-4731-b471-d5b15b1cc5d1)

Figure 3. Comparing the GC Content vs the Coverage of Microbacterium and the Span in (Kb) of Microbacterium vs No-Hit data using Blobtools (Pre-filtering) 

![blob_out blobDB json bestsum genus p8 span 100 blobplot read_cov bam0](https://github.com/omsmith161/Gen711-Final-Project/assets/158241303/ea10f1f8-0464-4d72-988d-4ecf10eec3b3)

Figure 4. Left Graph Comparing Mapped vs Unmapped Percentages of the Assembly. Right Graph Comparing the Percentages of No-Hit vs Microbacterium Data (Pre-filtering). Both used Blobtools. 


## Conclusions
Below is a shortened Quast data set for our fully filtered and assembled genome. It is clear that when comparing the two sets of Quast data that the filtering commands worked and cleaned up our data.   
## Final Filtered Quast Data
Assembly-------------microbacterium
- contigs (>= 0 bp)         26            
- contigs (>= 1000 bp)      25            
- contigs (>= 5000 bp)      22            
- contigs (>= 10000 bp)     19            
- contigs (>= 25000 bp)     17            
- contigs (>= 50000 bp)     14            
- Total length (>= 0 bp)      3090831              
- contigs                   26            

In comparing the two Busco results from before and after filtering, they were both 99.2% complete and were only missing 1 Busco group. 

## Final Busco Results:
***** Results: *****

        C:99.2%[S:99.2%,D:0.0%],F:0.0%,M:0.8%,n:124        
        123     Complete BUSCOs (C)                        
        123     Complete and single-copy BUSCOs (S)        
        0       Complete and duplicated BUSCOs (D)         
        0       Fragmented BUSCOs (F)                      
        1       Missing BUSCOs (M)                         
        124     Total BUSCO groups searched  


 ## Circular Genome
 - We were able to create a circular genome visual using our Prokka results with the 16s rRNA gene into a site called Prokksee.
<img width="590" alt="Screen Shot 2024-05-13 at 2 22 48 PM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/158241303/e0d30992-d541-43b6-8490-a2d47154f3cb">




