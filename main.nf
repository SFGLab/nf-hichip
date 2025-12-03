// Authors: Mateusz Chiliński (nextflow version) & Zofia Tojek (original version)

params.ref = "/app/ref/Homo_sapiens_assembly38.fasta"
params.ref_short = "hg38"
params.outdir = "results"
params.design = "/app/design_high.csv"
params.chrom_sizes = "/app/hg38.chrom.sizes"
params.threads = 8
params.mem = 4
params.mapq = 30
params.peak_quality = 0.05
params.genome_size = "hs"

all_chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
all_chromosomes_space = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
all_chromosomes_num = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

workflow {
    files = Channel.fromPath(params.design)
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }


    merged_fastq = files
        .map {
            meta, fastq1, fastq2 ->
                def meta_clone = meta.clone()
                meta_clone.remove('replicate')
                meta_clone.id = meta_clone.id
                meta_clone.chipseq = meta_clone.chipseq
                [ meta_clone, fastq1, fastq2 ]
        }
        .groupTuple(by: [0])

    MergeFiles(merged_fastq)
    Mapping(MergeFiles.out.sample, MergeFiles.out.fastq1, MergeFiles.out.fastq2)
    FilterQuality(Mapping.out.sample, Mapping.out.bam)
    RemoveDuplicates(Mapping.out.sample, FilterQuality.out.bam)
    CreateBigwig(Mapping.out.sample, RemoveDuplicates.out.bam)
    CallPeaks(Mapping.out.sample, MergeFiles.out.chipseq, RemoveDuplicates.out.bam)
    RunMapsSingleReplicate(Mapping.out.sample, MergeFiles.out.fastq1, MergeFiles.out.fastq2, CallPeaks.out.narrowPeak)
    CallHiCHeatMap(Mapping.out.sample, RunMapsSingleReplicate.out.input_file)
}

process MergeFiles {
    tag "Merging files"

    input:
    tuple val(meta), path(fastq1, stageAs: "?/*"), path(fastq2, stageAs: "?/*")

    output:
    val(meta.id), emit: sample
    val(meta.chipseq), emit: chipseq
    path("${meta.id}_sample_R1.fastq"), emit: fastq1
    path("${meta.id}_sample_R2.fastq"), emit: fastq2
 
    script:
    """
    if test ${fastq1[0]} = 1
    then
        cat 1/${fastq1[1]} > ${meta.id}_sample_R1.fastq
        cat 1/${fastq2[1]} > ${meta.id}_sample_R2.fastq
    else
        cat ${fastq1.join(' ')} > ${meta.id}_sample_R1.fastq
        cat ${fastq2.join(' ')} > ${meta.id}_sample_R2.fastq
    fi
    """
}

process Mapping {
    tag "Mapping files"

    input:
    val sample
    path fastq1
    path fastq2
 
    output:
    val(sample), emit: sample
    path("${sample}_output.bam"), emit: bam
 
    script:
    """
    set -euo pipefail

    bwa mem -M -v 0 -t ${params.threads} ${params.ref} ${fastq1} ${fastq2} \
        | samtools view -@ ${params.threads} -bh - \
        > ${sample}_output.bam
    """
}

process FilterQuality {

    input:
    val sample
    path mapped_bam
    
    output:
    path "${sample}_output_filtered.bam", emit: bam

    script:
    """
    set -euo pipefail

    # Remove unmapped reads; apply MAPQ filter
    samtools view \
        -@ ${params.threads} \
        -F 0x04 \
        -q ${params.mapq} \
        -b ${mapped_bam} \
        > ${sample}_output_filtered.bam
    """
}

process RemoveDuplicates {

    input:
    val sample
    path bam

    output:
    path "${sample}_dedup.bam", emit: bam
    publishDir "final_output/bam/"

    script:
    """
    set -euo pipefail

    # Compute per-thread memory chunk (in GB)
    # Distribute total mem across threads
    CHUNK_MEM=\$(( ${params.mem} / (${params.threads}) ))
    if [ "\$CHUNK_MEM" -lt 1 ]; then CHUNK_MEM=1; fi

    echo "Using samtools sort -m \${CHUNK_MEM}G per thread with ${params.threads} threads (total mem=${params.mem}G)"

    # 1) Name-sort for fixmate
    samtools sort -n -@ ${params.threads} -m \${CHUNK_MEM}G \
        -o ${sample}.nsort.bam \
        ${bam}

    # 2) Add mate information on name-sorted BAM
    samtools fixmate -m -@ ${params.threads} \
        ${sample}.nsort.bam \
        ${sample}.fixmate.bam

    # 3) Coordinate-sort for duplicate marking
    samtools sort -@ ${params.threads} -m \${CHUNK_MEM}G \
        -o ${sample}.possort.bam \
        ${sample}.fixmate.bam

    # 4) Mark/remove duplicates
    samtools markdup -r -@ ${params.threads} \
        ${sample}.possort.bam \
        ${sample}_dedup.bam

    # 5) Cleanup
    rm -f ${sample}.nsort.bam ${sample}.fixmate.bam ${sample}.possort.bam
    """
}

process CreateBigwig {

    input:
    val sample
    path final_bam

    output:
    path "${sample}_output.bigWig", emit: bigwig
    publishDir "final_output/coverage/"

    script:
    """
    set -euo pipefail

    # final_bam is already coordinate-sorted (from RemoveDuplicates)  
    samtools index -@ ${params.threads} ${final_bam}

    bamCoverage \
        -p ${params.threads} \
        -b ${final_bam} \
        -o ${sample}_output.bigWig
    """
}

process CallPeaks {

    input:
    val sample
    val chipseq
    path(final_bam)

    output:
    path "${sample}_peaks.narrowPeak", emit: narrowPeak
    
    publishDir "final_output/peaks/"

    script:
    """
    if test ${chipseq} = None
    then
        macs3 callpeak --nomodel -q ${params.peak_quality} -B -t ${final_bam} -n ${sample} -g ${params.genome_size} -f BAMPE
    else
        cp ${chipseq} ${sample}_peaks.narrowPeak
    fi
    
    """
}

process RunMapsSingleReplicate {

    input:
    val sample
    path(fastq1)
    path(fastq2)
    path(narrowPeak)

    publishDir "final_output/loops/", pattern: '*.bedpe'

    output:
    val(sample), emit: info
    path "${sample}.bedpe"
    path "${sample}.hic.input", emit: input_file

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
    /app/tasks/run_maps.sh > ${sample}_maps.txt
    mv MAPS_output/${sample}_current/${sample}.5k.2.sig3Dinteractions.bedpe .
    mv ${sample}.5k.2.sig3Dinteractions.bedpe ${sample}.bedpe
    mv feather_output/${sample}_current/${sample}.hic.input .
    """
}

process CallHiCHeatMap {

    input:
    val sample
    path(input_file)

    publishDir "final_output/hic/"

    output:
    path "${sample}.hic"

    script:
    """
    /opt/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar pre ${input_file} ${sample}.hic ${params.ref_short}
    """
}

def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.replicate           = row.replicate
    meta.chipseq           = row.chipseq

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    fastq_meta = [ meta, file(row.fastq_1), file(row.fastq_2) ]
    
    return fastq_meta
}
