#! /bin/bash

./align \
    -r file:/mnt/data/ref_genome/HanXRQr1.0-20151230.fa \
    -i file:/mnt/data/align.job \
    -x file:/mnt/data/ref_genome/PE-Marco.fa \
    -lossy -d 1
