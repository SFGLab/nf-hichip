#!/usr/bin/env bash
# for hichip
#bwa mem -SP5M -v 0 -t${THREADS} ${REFERENCE} ${R1} ${R2} | samtools view -bhS - > ${OUTNAME_MAPPED}
# dor chip-seq
bwa mem -M -v 0 -t${THREADS} ${REFERENCE} ${R1} ${R2} | samtools view -bh - > ${OUTNAME_MAPPED}