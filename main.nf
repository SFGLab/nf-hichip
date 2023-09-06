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
    // Channel.fromPath(params.design).splitCsv(header: true).map(row -> tuple(row.sample, row.replicate)).groupTuple().view()
    Mapping(files)
    RemoveNotAligned(Mapping.out.sample, Mapping.out.replicate, Mapping.out.bam)
    MappingQualityFilter(Mapping.out.sample, Mapping.out.replicate, RemoveNotAligned.out.bam)
    RemoveDuplicates(Mapping.out.sample, Mapping.out.replicate, MappingQualityFilter.out.bam)
    CreateBigwig(Mapping.out.sample, Mapping.out.replicate, RemoveDuplicates.out.bam)
    CallPeaks(Mapping.out.sample, Mapping.out.replicate, RemoveDuplicates.out.bam)
    RunMapsSingleReplicate(Mapping.out.sample, Mapping.out.replicate, Mapping.out.fastq1, Mapping.out.fastq2, CallPeaks.out.narrowPeak)
    RunMapsMultipleReplicate(RunMapsSingleReplicate.out.info.groupTuple(), CallPeaks.out.narrowPeak)
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
    path "${sample}_${replicate}_peaks.narrowPeak", emit: narrowPeak
    
    publishDir "final_output/"

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
    tuple val(sample), val(replicate), emit: info
    path "MAPS_output/"
    path "feather_output/"

    script:
    """
    export DATASET_NUMBER=1
    export DATASET_NAME=${sample}_${replicate}
    export fastq1=${fastq1}
    export fastq2=${fastq2}
    export OUTDIR=.
    export MACS_OUTPUT=${narrowPeak}
    export BWA_INDEX=${params.ref}
    export MAPQ=${params.mapq}
    export THREADS=${params.threads}
    mkdir MAPS_output/
    mkdir feather_output/
    mkdir MAPS_output/${sample}_${replicate}_current/
    touch MAPS_output/${sample}_${replicate}_current/${sample}_${replicate}.5k.2.sig3Dinteractions.bedpe
    workspaces/hichip-nf-pipeline/tasks/run_maps.sh > maps.txt
    """
}


process RunMapsMultipleReplicate {
    debug true

    input:
    val sample
    path(narrowPeak)

    publishDir "final_output/"

    output:
    stdout
    val sample
    file "done.txt"

    script:
    """
    #!/usr/bin/env python3
    import ast
    import os
    from plumbum import local

    sample_text = "$sample"
    sample_text = sample_text.replace("[", '["', 1)
    sample_text = sample_text.replace(",", '",', 1)
    sample = ast.literal_eval(sample_text)
    sample[1].sort()
    print(sample)
    str2 = 'local.env(DATASET_NUMBER=%s, DATASET_NAME="%s", OUTDIR=".", MACS_OUTPUT="/workspaces/hichip-nf-pipeline/final_output/%s", BWA_INDEX="%s", MAPQ="%s", THREADS="%s"' % (sample[1][-1], sample[0], sample[0]+"_"+str(sample[1][-1])+"_peaks.narrowPeak", "$params.ref", "$params.mapq", "$params.threads")
    for i in sample[1]:
        str2 += ', DATASET%s="/workspaces/hichip-nf-pipeline/final_output/feather_output/%s_%s_current/"' % (i, sample[0], i)
    
    str2 += ")"
    print(str2)
    with eval(str2):
        run_maps = local["/workspaces/hichip-nf-pipeline/tasks/run_maps.sh"]
        (run_maps > "done.txt")()
    """
}
