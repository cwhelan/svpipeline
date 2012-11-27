#!/bin/bash

set -e
set -u

HAP1_REF=$1
HAP2_REF=$2
OUTPUT_PREFIX=$3

# dwgsim params:
#
# -e and -E: 0.02 error rate
# -d, -s, -1, -2: 100bp reads from 300 bp fragments with stddev of 30
# -C: 30X coverage
# -r: 0.0010 mutation rate (default)
# -R: no indels (indels already applied to reference)
# -H: haploid mode

dwgsim -e 0.02-0.02 -E 0.02-0.02 -d 100 -s 30 -1 100 -2 100 -C 10 -r 0.0010 -R 0 -H $HAP1_REF $HAP1_REF
dwgsim -e 0.02-0.02 -E 0.02-0.02 -d 100 -s 30 -1 100 -2 100 -C 10 -r 0.0010 -R 0 -H $HAP2_REF $HAP2_REF

cat ${HAP1_REF}.bwa.read1.fastq ${HAP2_REF}.bwa.read1.fastq | gzip -c > ${OUTPUT_PREFIX}.read1.fastq.gz
cat ${HAP1_REF}.bwa.read2.fastq ${HAP2_REF}.bwa.read2.fastq | gzip -c > ${OUTPUT_PREFIX}.read2.fastq.gz

rm ${HAP1_REF}.bwa.read1.fastq
rm ${HAP1_REF}.bwa.read2.fastq
rm ${HAP1_REF}.bfastq.fastq

rm ${HAP2_REF}.bwa.read1.fastq
rm ${HAP2_REF}.bwa.read2.fastq
rm ${HAP2_REF}.bfastq.fastq
