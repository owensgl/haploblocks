#6199548 7199548 8699260 9699260 0 Ha412HOChr01 annuus Annuus env neg1
line=$1
start1=$(echo $line | cut -f 1 -d " " )
end1=$(echo $line | cut -f 2 -d " ")
start2=$(echo $line | cut -f 3 -d " ")
end2=$(echo $line | cut -f 4 -d " ")
chosen_genotype=$(echo $line | cut -f 5 -d " ")
chr=$(echo $line | cut -f 6 -d " ")
mds=$(echo $line | cut -f 10 -d " ")
species=$(echo $line | cut -f 7 -d " ")
samplelist="Annuus.tranche90.snp.env.90.bi.samplelist.txt"
vcf="Annuus.tranche90.snp.env.90.bi.remappedHa412HO.vcf.gz"

#echo $chr $mds 

cat /home/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/$species/Ha412HO_inv.v3.pcasites.$chr.$mds.genotypes.txt | grep -f $samplelist | grep $chosen_genotype$ | cut -f 1 > tmp.$chr.$mds.samplelist.txt

bcftools view -S tmp.$chr.$mds.samplelist.txt -r $chr $vcf| vcftools --vcf -  --maf 0.05 --thin 100  -c --geno-r2  | perl /home/owens/bin/reformat/emerald2windowldcounts.pl | gzip > $species.$chr.$mds.withinhaplo.ld.txt.gz

rm tmp.samplelist.$chr.$mds.txt



