UnO-nf is a nextflow pipeline for second tier of UnO workflow, the analysis of shotgun metagenomic seqeuncing reads
Steps for UnO Tier 2:
1. FastQC- QC of raw reads
2. Trimmomatic- Trimming illumina adapters and low quality sequences 
--Include baked in parameters--
3. FastQC- QC of trimmed reads 
4. MegaHit- Co-Assembly of reads from the outbreak
--include baked in paramters---
5. Bowtie2 Index and Mapping- Reads are mapped back to the co-assembly for MetaBat2 Binning
6. MaxBin2- 
7. MetaBat2-
8. MetaQUAST-
9. DasTool Contig2Bin and Refining
10. CheckM-
11. (In progress) Bowtie2 Indexing and Mapping-Mapping reads back to bins
12. Prodigal-
Usage:
nextflow run QC_Trim_Mega_BT_ST.nf -profile conda, sge
NextFlow version--- currently running on version 22.10.6
Input format-- currently set by params (directory containing filenames with wild cards) will be changed to --input flag

Output- sets of bins, gene predictions
MIDAS2 module is run seperately on trimmed reads but will soon be integrated into the UnO-nf workflow
Input format-currently set by params 
