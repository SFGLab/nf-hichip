// Authors: Mateusz Chiliński (nextflow version) & Zofia Tojek (original version)

params.ref = "/app/ref/Homo_sapiens_assembly38.fasta"
params.outdir = "results"
params.design = "/app/design_high.csv"
params.chrom_sizes = "/app/hg38.chrom.sizes"
params.threads = 4
params.mem = 16
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
            meta, fastq1, fastq2, input1, input2 ->
                def meta_clone = meta.clone()
                meta_clone.remove('replicate')
                meta_clone.id = meta_clone.id
                [ meta_clone, fastq1, fastq2, input1, input2 ]
        }
        .groupTuple(by: [0])

    MergeFiles(merged_fastq)
    Mapping(MergeFiles.out.sample, MergeFiles.out.fastq1, MergeFiles.out.fastq2, MergeFiles.out.input1, MergeFiles.out.input2)
    FilterQuality(Mapping.out.sample, Mapping.out.bam, Mapping.out.bam_input)
    RemoveDuplicates(Mapping.out.sample, FilterQuality.out.bam, FilterQuality.out.bam_input)
    CreateBigwig(Mapping.out.sample, RemoveDuplicates.out.bam, RemoveDuplicates.out.bam_input)
    CallPeaks(Mapping.out.sample, RemoveDuplicates.out.bam_input)
    RunMapsSingleReplicate(Mapping.out.sample,MergeFiles.out.fastq1,MergeFiles.out.fastq2,CallPeaks.out.narrowPeak)
}

process MergeFiles {
    tag "Merging files"

    input:
    tuple val(meta), path(fastq1, stageAs: "?/*"), path(fastq2, stageAs: "?/*"), path(input1, stageAs: "?/*"), path(input2, stageAs: "?/*")


    output:
    val(meta.id), emit: sample
    path("${meta.id}_sample_R1.fastq"), emit: fastq1
    path("${meta.id}_sample_R2.fastq"), emit: fastq2
    path("${meta.id}_input_sample_R1.fastq"), emit: input1
    path("${meta.id}_input_sample_R2.fastq"), emit: input2
 
    script:
    """
    if test ${fastq1[0]} = 1
    then
        cat 1/${fastq1[1]} > ${meta.id}_sample_R1.fastq
        cat 1/${fastq2[1]} > ${meta.id}_sample_R2.fastq
        cat 1/${input1[1]} > ${meta.id}_input_sample_R1.fastq
        cat 1/${input2[1]} > ${meta.id}_input_sample_R2.fastq
    else
    cat ${fastq1.join(' ')} > ${meta.id}_sample_R1.fastq
    cat ${fastq2.join(' ')} > ${meta.id}_sample_R2.fastq
    cat ${input1.join(' ')} > ${meta.id}_input_sample_R1.fastq
    cat ${input2.join(' ')} > ${meta.id}_input_sample_R2.fastq
    fi
    """
}

process Mapping {
    tag "Mapping files"

    input:
    val(sample)
    path(fastq1)
    path(fastq2)
    path(input1)
    path(input2)
 
    output:
    val(sample), emit: sample
    path("${sample}_output.bam"), emit: bam
    path("${sample}_input_output.bam"), emit: bam_input
 
    script:
    """
    bwa mem -M -v 0 -t ${params.threads} ${params.ref} ${fastq1} ${fastq2} | samtools view -@ ${params.threads} -bh - > ${sample}_output.bam
    bwa mem -M -v 0 -t ${params.threads} ${params.ref} ${input1} ${input2} | samtools view -@ ${params.threads} -bh - > ${sample}_input_output.bam
    """
}

process FilterQuality {

    input:
    val sample
    path(mapped_bam)
    path(bam_input)
    
    output:
    path "${sample}_output_filtered.bam", emit: bam
    path "${sample}_input_output_filtered.bam", emit: bam_input

    script:
    """
    samtools view -@ ${params.threads} -F 0x04 -b ${mapped_bam} > ${sample}_output_removed_not_aligned.bam
    samtools view -@ ${params.threads} -F 0x04 -b ${bam_input} > ${sample}_input_output_removed_not_aligned.bam
    samtools view -@ ${params.threads} -q ${params.mapq} -b ${sample}_output_removed_not_aligned.bam > ${sample}_output_filtered.bam
    samtools view -@ ${params.threads} -q ${params.mapq} -b ${sample}_input_output_removed_not_aligned.bam > ${sample}_input_output_filtered.bam
    """
}

process RemoveDuplicates {

    cpus params.threads
    memory "${params.mem} GB"

    input:
    val sample
    path bam
    path input_bam

    output:
    path "${sample}_dedup.bam", emit: bam
    path "${sample}_input_dedup.bam", emit: bam_input

    publishDir "final_output/bam/", mode: 'copy'

    script:
    """
    set -euo pipefail

    THREADS=${task.cpus}
    TOTAL_MEM_GB=${params.mem}

    CHUNK_MEM=\$(( TOTAL_MEM_GB / THREADS ))
    if [ "\$CHUNK_MEM" -lt 1 ]; then
        CHUNK_MEM=1
    fi

    echo "RemoveDuplicates using:"
    echo "  Threads: \${THREADS}"
    echo "  Total memory: \${TOTAL_MEM_GB}G"
    echo "  samtools sort memory per thread: \${CHUNK_MEM}G"

    echo "Processing ChIP/sample BAM: ${bam}"

    samtools sort -n -@ \${THREADS} -m \${CHUNK_MEM}G \
        -T ${sample}.tmp.nsort \
        -o ${sample}.nsort.bam \
        ${bam}

    samtools fixmate -m -@ \${THREADS} \
        ${sample}.nsort.bam \
        ${sample}.fixmate.bam

    samtools sort -@ \${THREADS} -m \${CHUNK_MEM}G \
        -T ${sample}.tmp.possort \
        -o ${sample}.possort.bam \
        ${sample}.fixmate.bam

    samtools markdup -r -@ \${THREADS} \
        ${sample}.possort.bam \
        ${sample}_dedup.bam


    echo "Processing input/control BAM: ${input_bam}"

    samtools sort -n -@ \${THREADS} -m \${CHUNK_MEM}G \
        -T ${sample}_input.tmp.nsort \
        -o ${sample}_input.nsort.bam \
        ${input_bam}

    samtools fixmate -m -@ \${THREADS} \
        ${sample}_input.nsort.bam \
        ${sample}_input.fixmate.bam

    samtools sort -@ \${THREADS} -m \${CHUNK_MEM}G \
        -T ${sample}_input.tmp.possort \
        -o ${sample}_input.possort.bam \
        ${sample}_input.fixmate.bam

    samtools markdup -r -@ \${THREADS} \
        ${sample}_input.possort.bam \
        ${sample}_input_dedup.bam


    rm -f ${sample}.nsort.bam ${sample}.fixmate.bam ${sample}.possort.bam
    rm -f ${sample}_input.nsort.bam ${sample}_input.fixmate.bam ${sample}_input.possort.bam
    rm -f ${sample}.tmp.nsort* ${sample}.tmp.possort*
    rm -f ${sample}_input.tmp.nsort* ${sample}_input.tmp.possort*
    """
}

process CreateBigwig {

    cpus params.threads
    memory "${params.mem} GB"

    input:
    val sample
    path(final_bam)
    path(input_final_bam)

    output:
    path "${sample}_output.bigWig", emit: bigwig
    path "${sample}_input_output.bigWig", emit: bigwig_input

    publishDir "final_output/coverage/", mode: 'copy'

    script:
    """
    set -euo pipefail

    THREADS=${task.cpus}

    samtools index -@ \${THREADS} ${final_bam}
    bamCoverage -p \${THREADS} \
        -b ${final_bam} \
        -o ${sample}_output.bigWig

    samtools index -@ \${THREADS} ${input_final_bam}
    bamCoverage -p \${THREADS} \
        -b ${input_final_bam} \
        -o ${sample}_input_output.bigWig
    """
}

process CallPeaks {

    input:
    val sample
    path chipseq_bam
    
    output:
    path "${sample}_peaks.narrowPeak", emit: narrowPeak
    
    publishDir "final_output/peaks/"

    script:
    """
    macs3 callpeak \
        --nomodel \
        -q ${params.peak_quality} \
        -B \
        -t ${chipseq_bam} \
        -n ${sample} \
        -g ${params.genome_size} \
        -f BAMPE
    """
}

process RunMapsSingleReplicate {

    tag "MAPS loop calling: ${sample}"

    cpus params.threads
    memory "${params.mem} GB"

    input:
    val sample
    path fastq1
    path fastq2
    path narrowPeak

    output:
    path "${sample}.bedpe", emit: loops
    path "${sample}_maps.txt", emit: maps_log

    publishDir "final_output/loops/", mode: 'copy'

    script:
    """
    set -euo pipefail

    export DATASET_NUMBER=1
    export DATASET_NAME=${sample}
    export fastq1=${fastq1}
    export fastq2=${fastq2}
    export OUTDIR=.
    export MACS_OUTPUT=${narrowPeak}
    export BWA_INDEX=${params.ref}
    export MAPQ=${params.mapq}
    export THREADS=${task.cpus}

    /app/tasks/run_maps.sh > ${sample}_maps.txt 2>&1

    if [ ! -f MAPS_output/${sample}_current/${sample}.5k.2.sig3Dinteractions.bedpe ]; then
        echo "ERROR: MAPS loop output not found:"
        echo "MAPS_output/${sample}_current/${sample}.5k.2.sig3Dinteractions.bedpe"
        echo "Check ${sample}_maps.txt for details."
        exit 1
    fi

    cp MAPS_output/${sample}_current/${sample}.5k.2.sig3Dinteractions.bedpe ${sample}.bedpe
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
    if (!file(row.input_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Input Read 1 FastQ file does not exist!\n${row.input_1}"
    }
    if (!file(row.input_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Input Read 2 FastQ file does not exist!\n${row.input_2}"
    }
    fastq_meta = [ meta, file(row.fastq_1), file(row.fastq_2), file(row.input_1), file(row.input_2) ]
    
    return fastq_meta
}
