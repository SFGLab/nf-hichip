params.ref = "/workspaces/hichip-nf-pipeline/ref/Homo_sapiens_assembly38.fasta"
params.outdir = "results"
params.design = "/workspaces/hichip-nf-pipeline/design_high.csv"
params.threads = 8
params.mem = 4
params.mapq = 30
params.peak_quality = 0.05
params.genome_size = "hs"

all_chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
all_chromosomes_space = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
all_chromosomes_num = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"


workflow {
    files = Channel.fromPath(params.design).splitCsv(header: true).map {row -> tuple(row.sample, row.replicate, row.fastq_1, row.fastq_2) }
    Mapping(files)
    RemoveNotAligned(Mapping.out.sample, Mapping.out.replicate, Mapping.out.bam)
    MappingQualityFilter(Mapping.out.sample, Mapping.out.replicate, RemoveNotAligned.out.bam)
    RemoveDuplicates(Mapping.out.sample, Mapping.out.replicate, MappingQualityFilter.out.bam)
    CreateBigwig(Mapping.out.sample, Mapping.out.replicate, RemoveDuplicates.out.bam)
    CallPeaks(Mapping.out.sample, Mapping.out.replicate, RemoveDuplicates.out.bam)
    RunMapsSingleReplicate(Mapping.out.sample, Mapping.out.replicate, Mapping.out.fastq1, Mapping.out.fastq2, CallPeaks.out.narrowPeak)
}

process Mapping {
    tag "Mapping files"
 
    input:
    tuple val(sample), val(replicate), path(fastq1), path(fastq2)
 
    output:
    path 'output.bam', emit: bam
    path fastq1, emit:fastq1
    path fastq2, emit:fastq2
    val sample, emit:sample
    val replicate, emit:replicate
    //val cram.simpleName, emit:sample
 
    script:
    """
    bwa mem -M -v 0 -t ${params.threads} ${params.ref} ${fastq1} ${fastq2} | samtools view -bh - > output.bam
    """
}

process RemoveNotAligned {

    input:
    val sample
    val replicate
    path(mapped_bam)
    
    output:
    path 'output_RemoveNotAlignedReads.bam', emit: bam

    script:
    """
    samtools view -F 0x04 -b ${mapped_bam} >  config.outnames["mapped"]] > output_RemoveNotAlignedReads.bam
    """
}

process MappingQualityFilter {

    input:
    val sample
    val replicate
    path(removed_bam)

    output:
    path 'output_quality.bam', emit: bam

    script:
    """
    samtools view -q ${params.mapq} -t ${params.threads} -b ${removed_bam} > output_quality.bam
    """
}

process RemoveDuplicates {

    input:
    val sample
    val replicate
    path(quality_bam)

    output:
    path 'output_final.bam', emit: bam

    script:
    """
    samtools sort -n -t ${params.threads} -m 1G ${quality_bam} -o - | samtools fixmate --threads ${params.threads} - - | samtools rmdup -S - output_final.bam
    """
}

process CreateBigwig {

    input:
    val sample
    val replicate
    path(final_bam)

    output:
    path 'output.bigWig', emit: bigwig

    script:
    """
    samtools sort -t ${params.threads} -m 16G ${final_bam} -o output_final_sorted.bam
    samtools index output_final_sorted.bam
    bamCoverage -b output_final_sorted.bam -o output.bigWig -p ${params.threads}
    """
}

process CallPeaks {

    input:
    val sample
    val replicate
    path(final_bam)

    output:
    path "${sample}_peaks.narrowPeak", emit: narrowPeak

    script:
    """
    macs3 callpeak --nomodel -q ${params.peak_quality} -B -t ${final_bam} -n ${sample} -g ${params.genome_size} -f BAMPE
    """
}

process RunMapsSingleReplicate {

    input:
    val sample
    val replicate
    path(fastq1)
    path(fastq2)
    path(narrowPeak)

    publishDir "final_output/"
    
    output:
    path "MAPS_output/${sample}_current/${sample}.5k.2.sig3Dinteractions.bedpe", emit: bedpe

    script:
    """
    export DATASET_NUMBER=1
    export DATASET_NAME=${sample}
    export fastq1=${fastq1}
    export fastq2=${fastq2}
    export OUTDIR=.
    export MACS_OUTPUT=${narrowPeak}
    export BWA_INDEX=${params.ref}
    export MAPQ=${params.mapq}
    export THREADS=${params.threads}
    /workspaces/hichip-nf-pipeline/tasks/run_maps.sh > maps.txt
    """
}