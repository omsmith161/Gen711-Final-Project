# Gen711-Final-Project

## Background
The Original Team's Data: 

We used data that was found in a shared directory on ron. 
Steps Used to Obtain Samples:
1. Isolating bacteria
2. DNA Extraction
3. Prepare Library
4. Sequence using Illumina HiSeq 2500
- there were forward & reverse reads for each sample in FASTQ format

The Goal: To have a complete genome assembly and assessment of the chosen reads from the given data.

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

## Conclusions
