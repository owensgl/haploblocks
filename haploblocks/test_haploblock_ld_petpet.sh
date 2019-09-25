#6199548 7199548 8699260 9699260 0 Ha412HOChr01 annuus Annuus env neg1
while read line
do
start1=$(echo $line | cut -f 1 -d " " )
end1=$(echo $line | cut -f 2 -d " ")
start2=$(echo $line | cut -f 3 -d " ")
end2=$(echo $line | cut -f 4 -d " ")
chosen_genotype=$(echo $line | cut -f 5 -d " ")
if [ "$chosen_genotype" == "NA" ]
then
	continue
fi
chr=$(echo $line | cut -f 6 -d " ")
mds=$(echo $line | cut -f 10 -d " ")
species=$(echo $line | cut -f 7 -d " ")
samplelist="Petiolaris.tranche90.snp.petpet.90.bi.samplelist.txt"
vcf="Petiolaris.tranche90.snp.petpet.90.bi.remappedHa412HO.vcf.gz"


cat /home/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/$species/Ha412HO_inv.v3.pcasites.$chr.$mds.genotypes.txt | grep -f $samplelist | grep $chosen_genotype$ | cut -f 1 > tmp.samplelist.txt
bcftools view -S tmp.samplelist.txt -r $chr $vcf| vcftools --vcf -  --maf 0.05 --thin 100  -c --geno-r2  | perl /home/owens/bin/reformat/emerald2windowldcounts.pl | gzip > $species.$chr.$mds.withinhaplo.ld.txt.gz

echo $chr $mds 
done < petpet.haploblockld.info
