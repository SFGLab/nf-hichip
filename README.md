<h2 align="center"> nf-HiChIP Pipeline </h2>

<p align="center">
<img align="center" src="https://github.com/SFGLab/nf-hichip/blob/main/nf_HiChIP_pipeline.png">
</p>

-------
## Introduction

We have developed an nf-HiChIP pipeline that combines the analytical approach for ChIP-seq data processing (mapping, filtering, peak calling, coverage tracks calculations) with HiChIP-specific analysis (MAPS pipeline, Juric, Ivan, et al.). This pipeline enables users to conduct thorough and efficient analysis of multiple HiChIP datasets simultaneously, eliminating the requirement for additional ChIP-seq experiments. This workflow is based on the reference implementation of the method designed by Zofia Tojek. The original version is available [here](https://github.com/Zojka/luigi_seq).


-------
## Working with nf-HiChIP pipeline

#### Step 0. (optional)
1) You can get familiar with Nextflow options.
2) ```-resume``` flag allows you to execute the pipeline from the last successful step. 
3) For more details, see [Nextflow documentation](https://www.nextflow.io/docs/latest/cache-and-resume.html).

#### Step 1.
[Docker](https://hub.docker.com/r/mateuszchilinski/hichip-nf-pipeline) image available:
```
https://hub.docker.com/repository/docker/mateuszchilinski/hichip-nf-pipeline/general
```
Command to run Docker image (use -v to bind folder with data):
```
docker run -v /path_to_your_data/:/data_in_container/ -it mateuszchilinski/hichip-nf-pipeline:latest bash
```
#### Step 2. 
Required Files for Reference Folder (Total 6 files) -
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

**Example 1 for design.csv file** 
> If you do not have raw and processed results (narrow peaks) from ChIP-Seq experiment

sample | fastq_1 |fastq_2 | replicate | chipseq
-- | ------ |------ | ------ | ------
S1 | /data/SAMPLE1_1_R1.fastq.gz | /data/SAMPLE1_1_R2.fastq.gz | 1 | None
S1 | /data/SAMPLE1_2_R1.fastq.gz | /data/SAMPLE1_2_R2.fastq.gz | 2 | None
S2 | /data/SAMPLE2_1_R1.fastq.gz | /data/SAMPLE2_1_R2.fastq.gz | 1 | None
S2 | /data/SAMPLE2_2_R1.fastq.gz | /data/SAMPLE2_2_R2.fastq.gz | 2 | None

**Note** - 
1) "None" (note the capital letter) in the last column.
2) In this case, pseudo-ChIP-Seq data will be generated from HiChIP data.

**Example 2 for design.csv file**
> If you have raw ChIP-Seq data but the peaks have not been called yet:

sample | fastq_1 |fastq_2 | input_1 | input_2 | replicate
-- | ------ |------ | ------ | ------ | --
S1 | /data/SAMPLE1_1_R1.fastq.gz | /data/SAMPLE1_1_R2.fastq.gz | /data/SAMPLE1_INPUT_R1.fastq.gz | /data/SAMPLE1_INPUT_R2.fastq.gz | 1
S1 | /data/SAMPLE1_2_R1.fastq.gz | /data/SAMPLE1_2_R2.fastq.gz | /data/SAMPLE1_INPUT_R1.fastq.gz | /data/SAMPLE1_INPUT_R2.fastq.gz | 2
S2 | /data/SAMPLE2_1_R1.fastq.gz | /data/SAMPLE2_1_R2.fastq.gz | /data/SAMPLE2_INPUT_R1.fastq.gz | /data/SAMPLE2_INPUT_R2.fastq.gz | 1
S2 | /data/SAMPLE2_2_R1.fastq.gz | /data/SAMPLE2_2_R2.fastq.gz | /data/SAMPLE2_INPUT_R1.fastq.gz | /data/SAMPLE2_INPUT_R2.fastq.gz | 2

**Example 3 for design.csv file**
> If you have processed ChIP-Seq experiment results (in the form of narrow peaks):

sample | fastq_1 |fastq_2 | replicate | chipseq
-- | ------ |------ | ------ | ------
S1 | /data/SAMPLE1_1_R1.fastq.gz | /data/SAMPLE1_1_R2.fastq.gz | 1 | /data/SAMPLE1.narrowPeak
S1 | /data/SAMPLE1_2_R1.fastq.gz | /data/SAMPLE1_2_R2.fastq.gz | 2 | /data/SAMPLE1.narrowPeak
S2 | /data/SAMPLE2_1_R1.fastq.gz | /data/SAMPLE2_1_R2.fastq.gz | 1 | /data/SAMPLE2.narrowPeak
S2 | /data/SAMPLE2_2_R1.fastq.gz | /data/SAMPLE2_2_R2.fastq.gz | 2 | /data/SAMPLE2.narrowPeak

**Note** -
1) Remember, the pipeline requires chromosome names in the "chrX" format (e.g., chr1, chr14, chr21) in the narrowpeak file.
2) Ensure peak files follow this naming convention and the BED6+4 format.


#### Step 4.
To run, use the command inside the container use: 
```
/opt/nextflow run main.nf --design design.csv
```

If you want to call chipseq peaks first, you can use main_chipseq.nf with exactly the same parameters. The input should be exactly the same for all replicates in a given sample.
```
/opt/nextflow run main_chipseq.nf --design design.csv
```

#### Step 5.
The parameters of the pipeline can be found in the following table. All of them are optional: 

Parameter | Description | Default |
-- | ------ |------ |
--ref | Reference genome for the analysis. | /workspaces/hichip-nf-pipeline/ref/Homo_sapiens_assembly38.fasta
--outdir | Folder with the final results. | results
--design | .csv file containing information about samples and replicates. | /workspaces/hichip-nf-pipeline/design_high.csv
--chrom_sizes | Sizes of chromosomes for the specific reference genome. | /workspaces/hichip-nf-pipeline/hg38.chrom.sizes
--threads | Threads are to be used in each task. | 4
--mem | Memory to use (in GB) for all samtools tasks (per-sample - e.g., 4 samples with 4 threads with 4GB would consume 64GB of memory). | 16
--mapq | MAPQ for MAPS. | 30
--peak_quality | Quality parameter (q-value (minimum FDR) cutoff) for MACS3. | 0.05
--genome_size | Genome size string for MACS3. | hs

#### Step 6.
For Post-processing and figure recreation, please follow the scripts in the folder [post_processing](https://github.com/SFGLab/nf-hichip/tree/main/post_processing)

-------
## Citation
If you use nf-HiChIP in your research (the idea, the algorithm, the analysis scripts, or the supplemental data), please give us a star on the GitHub repo page and cite our paper as follows:    

- Preprint bioRxiv : 
Jodkowska, K., Parteka-Tojek, Z., Agarwal, A., Denkiewicz, M., Korsak, S., Chiliński, M., Banecki, K., & Plewczynski, D. (2024). Improved cohesin HiChIP protocol and bioinformatic analysis for robust detection of chromatin loops and stripes. In bioRxiv (p. 2024.05.16.594268). https://doi.org/10.1101/2024.05.16.594268
-------
