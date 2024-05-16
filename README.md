# Gen711-Final-Project

## Background
The original reads that we used for this project came from sequencing seaweed-eating microbes that were collected from Acadia National Park. 
https://docs.google.com/presentation/d/1OBzO8tTlOovftlic2hYZ7EWInaVoaPCbp_dOabNGVWY/edit#slide=id.g13fdf044cda_11_289
 
Steps Used to Obtain Reads:
1. Isolating bacteria from the samples
2. DNA Extraction
3. Prepare Library
4. Sequence using Illumina HiSeq 2500

We copied a matching forward and reverse read from the gen711_project_data directory in tmp. 

The Goal of this Project: To assemble the bacterial genome with these reads, identify the bacterial species that the genome belongs to, and assess the quality of our assembly.

-The chosen forward & reverse reads:
36_S3_L001_R1_001.fastq.gz
36_S3_L001_R2_001.fastq.gz

## Methods
All of our code can be found in this GitHub repository in the .sh file titled "code-used.sh". Comments are included to help the reader follow along with wach step of the assembly. 

Our process for this assembly is outlined in the figure below. 
[Pipeline.pdf](https://github.com/omsmith161/Gen711-Final-Project/files/15328181/Pipeline.pdf)
Figure 1. Bioinformatics pipeline that we used for our genome assembly and analysis. 

We were guided by the sample genome assembly tutorial that was provided for us.
https://github.com/Joseph7e/MDIBL-T3-WGS-Tutorial

## Results 
Running trimmomatic successfully allowed us to trim the adapter sequences off of our forward and reverse reads. This was shown in our fastqc analysis. The most overrepresented sequence in the forward read was the TruSeq Adapter, also called the Illumina Universal Adapter. This sequence was not considered overrepresented in the reverse read, however, the Illumina Universal Adapter was also present in the "Adapter Content" section of the fastqc output. In both the forward and the reverse read, 0 sequences were flagged for poor quality, so trimmomatic was not needed for that purpose. After using trimmomatic, our additional fastqc analysis showed that it successfully removed adapter sequences since it was no longer flagging anything overrepresented. Additionally, the Illumina Universal Adapter became a lot less prevalent in the "Adapter Content" section of the output. 
<img width="972" alt="Screen Shot 2024-05-02 at 6 08 38 PM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/116675192/24aacbc9-b880-4e76-bd91-4289944b65ec">
<img width="936" alt="Screen Shot 2024-05-02 at 6 12 01 PM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/116675192/cc178afb-3ac0-4214-b924-a306ccb528b5">
Figure 2. Comparing the forward read Adapter Content before and after trimming adapters using trimmomatic. 


<img width="1281" alt="Screen Shot 2024-05-08 at 11 29 52 AM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/116675192/20a75691-1890-452e-964c-fe10e8725e15">
<img width="1300" alt="Screen Shot 2024-05-08 at 11 25 33 AM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/158241303/530b1c10-ebe7-4efc-be17-825b36cf59d5">

Figure 3. Comparing the reverse reads Adapter Content before and after trimming adapters using trimmomatic.

## Pre-Filtered Quast Data: 
Following our assembly of the trimmed forward and reverse reads using SPAdes, we used Quast and Busco to assess contiguity and completeness. Our Quast analysis showed that we had 28 contigs in total, the majority of which were greater than or equal to 1000 basepairs in length. The total length of the genome was 3,112,414 basepairs. 

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
After Quast, BUSCO was used to examine the completeness of this stage of assembly. Here, you can see that our genome was 99.2% complete, with only 1 missing BUSCO. We were pleased to see that our completeness was over 99%. 

***** Results: *****

	C:99.2%[S:99.2%,D:0.0%],F:0.0%,M:0.8%,n:124	   
	123	Complete BUSCOs (C)			   
	123	Complete and single-copy BUSCOs (S)	   
	0	Complete and duplicated BUSCOs (D)	   
	0	Fragmented BUSCOs (F)			   
	1	Missing BUSCOs (M)			   
	124	Total BUSCO groups searched

## Prokka
Prokka was used to annotate the genome with common bacterial genes. One of our outputs was a file called protein_abundances.txt, and this gave us a list of all the hypothetical proteins that could possibly come from the genes in this genome. The most important gene that we annotated with PROKKA was the 16s rRNA gene. This gene was used and inserted into the BLAST database. The 16s rRNA gene is very highly conserved in bacteria, and it can be used for species identification because this region is highly variable between species.

## BLAST Results
This was our output from running BLAST on the NCBI website. After PROKKA annotation, we pulled just the 16s rRNA gene sequence and copied and pasted that into BLAST. From this analysis, we saw that our 16s rRNA gene matched with both __Microbacterium oryzae_ and __Microbacterium barkeri_. These two bacteria species are probably very closely related, so it is not surprising that our gene matched with more than one species. We chose to identify our species as _Microbacterium barkeri___ because it had a greater percent match and came up more frequently. 
<img width="695" alt="Screen Shot 2024-05-13 at 2 05 45 PM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/158241303/50c46d6b-0c94-49df-b080-c14905fd198f">

## Blobtools 

![blob_out blobDB json bestsum genus p8 span 100 blobplot bam0](https://github.com/omsmith161/Gen711-Final-Project/assets/158241303/f20add11-23fc-4731-b471-d5b15b1cc5d1)

Figure 4. Comparing the GC Content vs the Coverage of Microbacterium and the Span in (Kb) of Microbacterium vs No-Hit data using Blobtools (Pre-filtering) 

![blob_out blobDB json bestsum genus p8 span 100 blobplot read_cov bam0](https://github.com/omsmith161/Gen711-Final-Project/assets/158241303/ea10f1f8-0464-4d72-988d-4ecf10eec3b3)

Figure 5. Left Graph Comparing Mapped vs Unmapped Percentages of the Assembly. Right Graph Comparing the Percentages of No-Hit vs Microbacterium Data (Pre-filtering). Both used Blobtools. 


## Conclusions
Below is a shortened Quast data set for our fully filtered and assembled genome. It is clear that when comparing the two sets of Quast data that the filtering commands worked and cleaned up our data. Filtering removed 2 contigs, we started with 28 and ended with 26 contigs. This brought the total length of our genome down to 3,090,831 basepairs. 

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

In comparing the two Busco results from before and after filtering, they were both 99.2% complete and were only missing 1 Busco group. This means the filtering steps did not help with improving the completeness of our genome, but it did help with improving contiguity, as was seen in the Quast results above. 

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
 The final step of our analysis was to visualize the genome. We used Proksee alongside one of our PROKKA files because it contained the annotations of all of the genes in the genome (Grant et al., 2023). The circular genome for our species can be found below. 
<img width="590" alt="Screen Shot 2024-05-13 at 2 22 48 PM" src="https://github.com/omsmith161/Gen711-Final-Project/assets/158241303/e0d30992-d541-43b6-8490-a2d47154f3cb">

## Resources 
Blast Web Data:
- https://blast.ncbi.nlm.nih.gov/Blast.cgi#

Original Study’s Data: 
- https://docs.google.com/presentation/d/1OBzO8tTlOovftlic2hYZ7EWInaVoaPCbp_dOabNGVWY/edit#slide=id.g13fdf044cda_11_289

Github Tutorial:
- https://github.com/Joseph7e/MDIBL-T3-WGS-Tutorial 

Papers:
- Grant, J. R., Enns, E., Marinier, E., Mandal, A., Herman, E. K., Chen, C., Graham, M., Van 
Domselaar, G., & Stothard, P. (2023). Proksee: in-depth characterization and visualization of bacterial genomes. Nucleic Acids Research, 51(W1), 
W484–W492. https://doi.org/10.1093/nar/gkad326

- Sanchez, S., Snider, E. v., Wang, X., & Kearns, D. B. (2022). Identification of Genes 
Required for Swarming Motility in Bacillus subtilis Using Transposon Mutagenesis and High-Throughput Sequencing (TnSeq). Journal of Bacteriology, 
204(6). https://doi.org/10.1128/jb.00089-22

- P0A6Y8 · DNAK_ECOLI. (n.d.). UniProt.
https://www.uniprot.org/uniprotkb/P0A6Y8/entry





