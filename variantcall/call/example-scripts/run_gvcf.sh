#! /bin/bash

./vc -n 8 -g \
    -r file:/mnt/ref_genome/HanXRQr1.0-20151230.fa \
    -b file:/mnt/ref_genome/HanXRQr1.0-20151230_allTEs_ubc.bed \
    -i file:/mnt/variants/DBGBGBS_GB331.bam \
    -gatk4
