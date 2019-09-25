full_species="petiolaris"
capital_species="Petiolaris"
short_species="petfal"
group_tag="petfal"
inversion_version="v3"

vcftools --vcf /home/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/$full_species/Ha412HO_inv.${inversion_version}.pcasites.vcf --keep $short_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.samplelist.txt --recode --out $short_species/Ha412HO_inv.${inversion_version}.pcasites.$short_species

/home/owens/bin/vcftools_0.1.12a/perl/vcf-shuffle-cols -t $short_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.vcf.gz $short_species/Ha412HO_inv.${inversion_version}.pcasites.$short_species.recode.vcf  | bgzip  > $short_species/Ha412HO_inv.${inversion_version}.pcasites.$short_species.shuffled.vcf.gz
tabix -p vcf $short_species/Ha412HO_inv.${inversion_version}.pcasites.$short_species.shuffled.vcf.gz

zcat $short_species/Ha412HO_inv.${inversion_version}.pcasites.$short_species.shuffled.vcf.gz | perl /home/owens/bin/reformat/vcf2tped.pl $short_species/Ha412HO_inv.${inversion_version}.pcasites.$short_species.shuffled


zcat $short_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.vcf.gz | perl /home/owens/bin/reformat/vcf2tped.pl $short_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle

zcat $short_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv.ldfilter.vcf.gz | perl /home/owens/bin/reformat/vcf2tped.pl $short_species/${capital_species}.tranche90.snp.${group_tag}.90.bi.remappedHa412HO.beagle.${inversion_version}noinv.ldfilter


