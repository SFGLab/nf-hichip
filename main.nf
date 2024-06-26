// Authors: Mateusz Chiliński (nextflow version) & Zofia Tojek (original version)

params.ref = "/app/hichip-nf-pipeline/ref/Homo_sapiens_assembly38.fasta"
params.outdir = "results"
params.design = "/app/hichip-nf-pipeline/design_high.csv"
params.chrom_sizes = "/app/hichip-nf-pipeline/hg38.chrom.sizes"
params.threads = 4
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
    ParsePairtools(Mapping.out.sample, FilterQuality.out.bam)
    RemoveDuplicates(Mapping.out.sample, ParsePairtools.out.pairsam)
    MakeFinalBam(Mapping.out.sample, RemoveDuplicates.out.pairsam)
    CreateBigwig(Mapping.out.sample, MakeFinalBam.out.bam)
    CallPeaks(Mapping.out.sample, MergeFiles.out.chipseq, MakeFinalBam.out.bam)
    RunMapsSingleReplicate(Mapping.out.sample, MergeFiles.out.fastq1, MergeFiles.out.fastq2, CallPeaks.out.narrowPeak)
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
    val(sample)
    path(fastq1)
    path(fastq2)
 
    output:
    val(sample), emit: sample
    path("${sample}_output.bam"), emit: bam
 
    script:
    """
    bwa mem -M -v 0 -t ${params.threads} ${params.ref} ${fastq1} ${fastq2} | samtools view -bh - > ${sample}_output.bam
    """
}

process FilterQuality {

    input:
    val sample
    path(mapped_bam)
    
    output:
    path "${sample}_output_filtered.bam", emit: bam

    script:
    """
    samtools view -q 30 -t ${params.threads} -b ${mapped_bam} > ${sample}_output_filtered.bam
    """
}

process ParsePairtools {

    input:
    val sample
    path(quality_bam)
    
    output:
    path "${sample}_paired.pairsam", emit: pairsam

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir temp
    samtools view -h ${quality_bam} | pairtools parse -c ${params.chrom_sizes} --add-columns mapq | pairtools sort --nproc 5 --memory 8G --tmpdir temp/ --output ${sample}_paired.pairsam
    """
}

process RemoveDuplicates {

    input:
    val sample
    path(pairsam)

    output:
    path "${sample}_dedup.pairsam", emit: pairsam

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    pairtools dedup --mark-dups --output-stats ${sample}_dedup.stats --output ${sample}_dedup.pairsam ${pairsam}
    """
}

process MakeFinalBam {

    input:
    val sample
    path(dedup_pairsam)

    output:
    path "${sample}_output_final.bam", emit: bam

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    pairtools split --output-sam ${sample}_output_final.bam ${dedup_pairsam}
    """
}

process CreateBigwig {

    input:
    val sample
    path(final_bam)

    output:
    path "${sample}_output.bigWig", emit: bigwig
    publishDir "final_output/coverage/"

    script:
    """
    samtools sort -t ${params.threads} -m ${params.mem}G ${final_bam} -o ${sample}_output_final_sorted.bam
    samtools index ${sample}_output_final_sorted.bam
    bamCoverage -b ${sample}_output_final_sorted.bam -o ${sample}_output.bigWig -p ${params.threads}
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

    publishDir "final_output/loops/"

    output:
    val(sample), emit: info
    path "${sample}.5k.2.sig3Dinteractions.bedpe"

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
    /app/hichip-nf-pipeline/tasks/run_maps.sh > ${sample}_maps.txt
    mv MAPS_output/${sample}_current/${sample}.5k.2.sig3Dinteractions.bedpe .
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
