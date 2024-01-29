## Introduction
UnO-nf is draft pipeline for the second tier of the UnO workflow, the analysis of shotgun metagenomic seqeuncing reads. 
UnO-nf is built using Nextflow. This pipeline is under active development. Future releases of this pipeline will include containerized versions of the tools utilized in this pipeline. At present tools for this pipeline are installed via UnO.yml. 
## Running UnO-nf 
1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) ('22.10.6')
2. Download UnO-nf source code for latest release
2. Install DASTool via conda [https://github.com/cmks/DAS_Tool]
Add bioconda channel:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
Install DAS Tool using conda:
```
conda install -c bioconda das_tool
```
3. Install MIDAS2 from source [https://midas2.readthedocs.io/en/latest/installation.html]
```
git clone https://github.com/czbiohub/MIDAS2.git
cd MIDAS2

conda env create -n midas2.0 -f midas2.yml
conda activate midas2
cpanm Bio::SearchIO::hmmer --force # Temporary fix for Prokka

pip install .

```
## How to run UnO-nf 
nextflow run UnO.nf -profile conda, sge
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
