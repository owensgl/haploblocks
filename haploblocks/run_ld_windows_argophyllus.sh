chr=$1

species="Argophyllus"
prefix="tranche90.snp.gwas.90.bi.remappedHa412HO"

#bcftools view --threads 2 -r Ha412HOChr${chr} -i 'MAF[0]>0.05' -O z $species.${prefix}.beagle.vcf.gz > $species.${prefix}.beagle.Chr${chr}.tmp.vcf.gz
#tabix -p vcf $species.${prefix}.beagle.Chr${chr}.tmp.vcf.gz
#/home/owens/bin/emeraLD/bin/emeraLD --phased -i  $species.${prefix}.beagle.Chr${chr}.tmp.vcf.gz   --stdout --window 1000000000 | perl /home/owens/bin/reformat/emerald2windowldcounts.pl | gzip >$species.${prefix}.Chr${chr}.windows.ld.gz
bcftools view -O v -r Ha412HOChr$chr $species.${prefix}.vcf.gz | vcftools --vcf -  --maf 0.05 --thin 100  -c --geno-r2 | perl /home/owens/bin/reformat/emerald2windowldcounts.pl | gzip >$species.${prefix}.thin100.maf5.Chr${chr}.windows.ld.gz
#rm $species.${prefix}.beagle.Chr${chr}.tmp.vcf.gz
