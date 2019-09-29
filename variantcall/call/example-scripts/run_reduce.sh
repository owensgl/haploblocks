#! /bin/bash

./reduce \
    -r file:/mnt/ref_genome/HanXRQr1.0-20151230.fa \
    -L file:/mnt/variants/reduce_sample_list \
    -cat -gatk4
