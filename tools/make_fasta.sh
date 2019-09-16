sample=$4 #Passed from parallel
samplefile=$1
ref=$2
script_bin=$3
bam=$(cat $samplefile | grep $sample | cut -f 1)

gatk HaplotypeCaller -R $ref -I $bam --intervals gene.bed -ERC BP_RESOLUTION -O $sample.g.vcf

cat $sample.g.vcf | perl $script_bin/gvcf2fasta_nogaps.pl > $sample.fa
