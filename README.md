<h2 align="center"> nf-HiChIP Pipeline </h2>

<p align="center">
<img align="center" src="https://github.com/SFGLab/nf-hichip/blob/main/nf_HiChIP_pipeline.png">
</p>

-------
## Introduction

We have developed an nf-HiChIP pipeline that combines the analytical approach designed for ChIP-seq data processing (mapping, filtering, peak calling, coverage tracks calculations) with HiChIP-specific analysis (MAPS pipeline, Juric, Ivan, et al.). This pipeline enables users to conduct thorough and efficient analysis of multiple HiChIP datasets simultaneously, eliminating the requirement for additional ChIP-seq experiments. This workflow is based on the reference implementation of the method designed by Zofia Tojek. The original version is available [here](https://github.com/Zojka/luigi_seq).


-------
## Working with nf-HiChIP pipeline

#### Step 1.
[Docker](https://hub.docker.com/) image available:
```
https://hub.docker.com/repository/docker/mateuszchilinski/hichip-nf-pipeline/general
```
Command to run Docker image (use -v to bind folder with data):
```
docker run -v /path_to_your_data/:/data_in_container/ -it mateuszchilinski/hichip-nf-pipeline:latest bash
```
#### Step 2. 
Required Files for Reference Folder (total 6 files) -
```
1. Reference fasta files -
    > Homo_sapiens_assembly38.fasta

2. BWA Reference Index files -
    > Homo_sapiens_assembly38.fasta.amb
    > Homo_sapiens_assembly38.fasta.ann
    > Homo_sapiens_assembly38.fasta.bwt
    > Homo_sapiens_assembly38.fasta.pac
    > Homo_sapiens_assembly38.fasta.sa
```

#### Step 3.
To run, use the command inside the container use: 

```
/opt/nextflow run main.nf --design design.csv
```
#### Step 4.
Example for design.csv file:

sample | fastq_1 |fastq_2 | replicate | chipseq
-- | ------ |------ | ------ | ------
S1 | /data/SAMPLE1_1_R1.fastq.gz | /data/SAMPLE1_1_R2.fastq.gz | 1 | None
S1 | /data/SAMPLE1_2_R1.fastq.gz | /data/SAMPLE1_2_R2.fastq.gz | 2 | None
S2 | /data/SAMPLE2_1_R1.fastq.gz | /data/SAMPLE2_1_R2.fastq.gz | 1 | /data/SAMPLE2.narrowPeak
S2 | /data/SAMPLE2_2_R1.fastq.gz | /data/SAMPLE2_2_R2.fastq.gz | 2 | /data/SAMPLE2.narrowPeak

Very important information regarding chipseq data!

If you don't have additional chipseq experiment results (in form of narrowPeaks), you need to put "None" in the last column - and the pseudo-chipseq will be called from HiChIP data. Also, remember that the pipeline is using chromosomes in "chrX" form. Peaks files need to be consistent with that as well.

If you have chipseq data, but don't have the peaks called yet, you can use ChIP-Seq processing part of the pipeline by calling:

```
/opt/nextflow run main_chipseq.nf --design design_chipseq.csv
```

In that case, the design_chipseq.csv should be in format of:

sample | fastq_1 |fastq_2 | input_1 | input_2 | replicate
-- | ------ |------ | ------ | ------ | --
S1 | /data/SAMPLE1_1_R1.fastq.gz | /data/SAMPLE1_1_R2.fastq.gz | /data/SAMPLE1_INPUT_R1.fastq.gz | /data/SAMPLE1_INPUT_R2.fastq.gz | 1
S1 | /data/SAMPLE1_2_R1.fastq.gz | /data/SAMPLE1_2_R2.fastq.gz | /data/SAMPLE1_INPUT_R1.fastq.gz | /data/SAMPLE1_INPUT_R2.fastq.gz | 2
S2 | /data/SAMPLE2_1_R1.fastq.gz | /data/SAMPLE2_1_R2.fastq.gz | /data/SAMPLE2_INPUT_R1.fastq.gz | /data/SAMPLE2_INPUT_R2.fastq.gz | 1
S2 | /data/SAMPLE2_2_R1.fastq.gz | /data/SAMPLE2_2_R2.fastq.gz | /data/SAMPLE2_INPUT_R1.fastq.gz | /data/SAMPLE2_INPUT_R2.fastq.gz | 2

If you want to call chipseq peaks first, you can use main_chipseq.nf with exactly the same parameters. The input should be exactly the same for all replicates in a given sample.

#### Step 5.
The parameters of the pipeline can be found in the following table. All of them are optional: 

Parameter | Description | Default |
-- | ------ |------ |
--ref | Reference genome for the analysis. | /workspaces/hichip-nf-pipeline/ref/Homo_sapiens_assembly38.fasta
--outdir | Folder with the final results. | results
--design | .csv file containing information about samples and replicates. | /workspaces/hichip-nf-pipeline/design_high.csv
--chrom_sizes | Sizes of chromosomes for the specific reference genome. | /workspaces/hichip-nf-pipeline/hg38.chrom.sizes
--threads | Threads to use in each task. | 4
--mem | Memory to use (in GB) for all samtools tasks (per-sample - e.g., 4 samples with 4 threads with 4GB would result in consumption of 64GB of memory). | 16
--mapq | MAPQ for MAPS. | 30
--peak_quality | Quality parameter (q-value (minimum FDR) cutoff) for MACS3. | 0.05
--genome_size | Genome size string for MACS3. | hs

#### Step 5.
For Post-processing and figure recreation, please follow the scripts in the folder [post_processing](https://github.com/SFGLab/nf-hichip/tree/main/post_processing)

-------
## Citation
If you use nf-HiChIP in your research (the idea, the algorithm, the analysis scripts, or the supplemental data), please give us a star on the GitHub repo page and cite our paper as follows:    

- Official version  --
- Preprint bioRxiv --
-------
