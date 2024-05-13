# nf-HiChIP Pipeline
Pipeline for processing HiChIP data.

Docker image available: https://hub.docker.com/repository/docker/mateuszchilinski/hichip-nf-pipeline/general

The pipeline is presented on the following figure:

![Figure 1](nf_HiChIP.pdf)

To use the pipeline, you need to have nextflow installed.

To run, use command: 

```nextflow run main.nf --design design.csv```

Example design.csv file:

sample | fastq_1 |fastq_2 | replicate |
-- | ------ |------ | ------ |
S1 | /data/SAMPLE1_1_R1.fastq.gz | /data/SAMPLE1_1_R2.fastq.gz | 1
S1 | /data/SAMPLE1_2_R1.fastq.gz | /data/SAMPLE1_2_R2.fastq.gz | 2
S2 | /data/SAMPLE2_1_R1.fastq.gz | /data/SAMPLE2_1_R2.fastq.gz | 1
S2 | /data/SAMPLE2_2_R1.fastq.gz | /data/SAMPLE2_2_R2.fastq.gz | 2

The parameters to the pipeline can be found in the following table. All of them are optional.


Parameter | Description | Default |
-- | ------ |------ |
--ref | Reference genome for the analysis. | /workspaces/hichip-nf-pipeline/ref/Homo_sapiens_assembly38.fasta
--outdir | Folder with the final results. | results
--design | Csv file containing information about samples and replicates. | /workspaces/hichip-nf-pipeline/design_high.csv
--chrom_sizes | Sizes of chromosomes for the specific reference genome. | /workspaces/hichip-nf-pipeline/hg38.chrom.sizes
--threads | Threads to use in each task. | 4
--mem | Memory to use (in GB) for sorting task. | 4
--mapq | MAPQ for MAPS. | 30
--peak_quality | Quality parameter for MACS3. | 0.05
--genome_size | Genome size string for MACS3. | hs 

