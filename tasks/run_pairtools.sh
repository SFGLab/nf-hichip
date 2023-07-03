#!/usr/bin/env bash

samtools view -h $MAPPED | pairtools parse -c $CHROMOSOMES --add-columns mapq | pairtools sort --nproc 5 --tmpdir ${TMPDIR} --memory 8G --output ${SORTED_PAIRS_PATH}