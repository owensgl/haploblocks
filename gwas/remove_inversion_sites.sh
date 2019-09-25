full_species="petiolaris"
capital_species="Petiolaris"
group_tag="petfal"
inversion_version="v3"

#Remove the main inversion SNPs
bcftools view --threads 5 -T ^/home/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.${inversion_version}.inversions.regions.${full_species}.bed ../${full_species}/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.vcf.gz -O z > $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp.vcf.gz
tabix -p vcf $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp.vcf.gz

zcat $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp.vcf.gz | head -n 90000 | grep CHR | cut -f 10- | tr '\t' '\n' > $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.samplelist.txt
#Add inversion genotypes to the VCF to use correlation to remove more SNPs
vcftools --vcf /home/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/$full_species/Ha412HO_inv.${inversion_version}.pcasites.vcf --keep $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.samplelist.txt --recode --out $full_species/Ha412HO_inv.${inversion_version}.pcasites.$full_species

/home/owens/bin/vcftools_0.1.12a/perl/vcf-shuffle-cols -t $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp.vcf.gz $full_species/Ha412HO_inv.${inversion_version}.pcasites.$full_species.recode.vcf  | bgzip  > $full_species/Ha412HO_inv.${inversion_version}.pcasites.$full_species.shuffled.vcf.gz
tabix -p vcf $full_species/Ha412HO_inv.${inversion_version}.pcasites.$full_species.shuffled.vcf.gz
bcftools concat --threads 6 -O z -a $full_species/Ha412HO_inv.${inversion_version}.pcasites.$full_species.shuffled.vcf.gz  $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp.vcf.gz >  $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp2.vcf.gz

vcftools --gzvcf $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp2.vcf.gz --geno-r2-positions /home/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/$full_species/Ha412HO_inv.${inversion_version}.pcasites.vcf --out $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp2


awk '($6 > 0.5)' $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp2.list.geno.ld | grep -v "nan" | cut -f 3,4 | sed '1d' | sort | uniq > $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp2.removesites.txt


bcftools view -T ^$full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp2.removesites.txt -O z $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp.vcf.gz > $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv.ldfilter.vcf.gz
tabix -p vcf $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv.ldfilter.vcf.gz

rm $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp.vcf.gz*
rm $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.${inversion_version}noinv_tmp2.vcf.gz*

#####################3
#Repeat with beagle imputed data

#Remove the main inversion SNPs
bcftools view --threads 5 -T ^/home/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.${inversion_version}.inversions.regions.${full_species}.bed ../${full_species}/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.vcf.gz -O z > $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp.vcf.gz
tabix -p vcf $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp.vcf.gz

#Add inversion genotypes to the VCF to use correlation to remove more SNPs

bcftools concat --threads 6 -O z -a $full_species/Ha412HO_inv.${inversion_version}.pcasites.$full_species.shuffled.vcf.gz  $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp.vcf.gz >  $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp2.vcf.gz

vcftools --gzvcf $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp2.vcf.gz --geno-r2-positions /home/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/$full_species/Ha412HO_inv.${inversion_version}.pcasites.vcf --out $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp2


awk '($6 > 0.5)' $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp2.list.geno.ld | grep -v "nan"| cut -f 3,4 | sed '1d' | sort | uniq > $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp2.removesites.txt


bcftools view -T ^$full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp2.removesites.txt -O z $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp.vcf.gz > $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv.ldfilter.vcf.gz
tabix -p vcf $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv.ldfilter.vcf.gz

zcat $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv.ldfilter.vcf.gz | perl /home/owens/bin/reformat/vcf2tped.pl $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv.ldfilter
rm $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp.vcf.gz*
rm $full_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv_tmp2.vcf.gz*
