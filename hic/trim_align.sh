input=$1 #Input without the R1.fastq.gz part
/project/rpp-rbruskie/gowens/HiC/homer/bin/homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 -q 15 ${input}_R1.fastq.gz

mv ${input}_R1.fastq.gz.trimmed ${input}_R1.trimmed.fastq
pigz -p 30 ${input}_R1.trimmed.fastq

/home/gowens/bin/NextGenMap-0.5.4/bin/ngm-0.5.4/ngm -r /scratch/gowens/ref/Ha412HOv2.0-20181130.onlychr.fasta -t 30 -q ${input}_R1.trimmed.fastq.gz -o ${input}_R1.sam; 

/project/rpp-rbruskie/gowens/HiC/homer/bin/homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 -q 15 ${input}_R2.fastq.gz

mv ${input}_R2.fastq.gz.trimmed ${input}_R2.trimmed.fastq
pigz -p 30 ${input}_R2.trimmed.fastq

/home/gowens/bin/NextGenMap-0.5.4/bin/ngm-0.5.4/ngm -r /scratch/gowens/ref/Ha412HOv2.0-20181130.onlychr.fasta -t 30 -q ${input}_R2.trimmed.fastq.gz -o ${input}_R2.sam;
