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
    CallPeaks(Mapping.out.sample, RemoveDuplicates.out.bam, RemoveDuplicates.out.bam_input)
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
        cat ${input1[1]} > ${meta.id}_input_sample_R1.fastq
        cat ${input2[1]} > ${meta.id}_input_sample_R2.fastq
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

    input:
    val sample
    path(bam)
    path(input_bam)

    output:
    path "${sample}_dedup.bam", emit: bam
    path "${sample}_input_dedup.bam", emit: bam_input

    script:
    """
    samtools sort -n -@ ${params.threads} -m ${params.mem}G ${bam} -o - | samtools fixmate -@ ${params.threads} - - | samtools rmdup -S - ${sample}_dedup.bam
    samtools sort -n -@ ${params.threads} -m ${params.mem}G ${input_bam} -o - | samtools fixmate -@ ${params.threads} - - | samtools rmdup -S - ${sample}_input_dedup.bam
    """
}

process CreateBigwig {

    input:
    val sample
    path(final_bam)
    path(input_final_bam)

    output:
    path "${sample}_output.bigWig", emit: bigwig
    path "${sample}_input_output.bigWig", emit: bigwig_input
    publishDir "final_output/coverage/"

    script:
    """
    samtools sort -@ ${params.threads} -m ${params.mem}G ${final_bam} -o ${sample}_output_final_sorted.bam
    samtools index -@ ${params.threads} ${sample}_output_final_sorted.bam
    bamCoverage -p ${params.threads} -b ${sample}_output_final_sorted.bam -o ${sample}_output.bigWig
    samtools sort -@ ${params.threads} -m ${params.mem}G ${input_final_bam} -o ${sample}_input_output_final_sorted.bam
    samtools index -@ ${params.threads} ${sample}_input_output_final_sorted.bam
    bamCoverage -p ${params.threads} -b ${sample}_input_output_final_sorted.bam -o ${sample}_input_output.bigWig
    """
}

process CallPeaks {

    input:
    val sample
    path(final_bam)
    path(input_final_bam)
    
    output:
    path "${sample}_peaks.narrowPeak", emit: narrowPeak
    
    publishDir "final_output/peaks/"

    script:
    """
    macs3 callpeak --nomodel -q ${params.peak_quality} -B -t ${final_bam} -c ${input_final_bam} -n ${sample} -g ${params.genome_size} -f BAMPE
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
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.input_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    fastq_meta = [ meta, file(row.fastq_1), file(row.fastq_2), file(row.input_1), file(row.input_2) ]
    
    return fastq_meta
}
