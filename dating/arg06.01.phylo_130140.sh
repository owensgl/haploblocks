## Run haploblock phylogeny using SNPs within. 
## run as: bash arg06.01.phylo.sh
## Kaichi Huang 2019 Apr

inv_n=15
sv_name="arg06.01"

wild_gwas_git="/scratch/gowens/wild_gwas/wild_gwas_2018"
inversion_version="v3"

inv_page="$wild_gwas_git/MDS_outliers/Ha412HO/all/Ha412HO_inv.$inversion_version.txt"
species=$(cat $inv_page | sed -n ${inv_n}p | cut -f 1)
chr_n=$(cat $inv_page | sed -n ${inv_n}p | cut -f 4)
chr_padded=$(printf %02d $chr_n)
chr="Ha412HOChr$chr_padded"
direction=$(cat $inv_page | sed -n ${inv_n}p | cut -f 6)
mds_n=$(cat $inv_page | sed -n ${inv_n}p | cut -f 5)
mds=${direction}$mds_n

cat $wild_gwas_git/MDS_outliers/Ha412HO/$species/Ha412HO_inv.$inversion_version.pcasites.$chr.$mds.genotypes.txt | awk '($2 != "NA") && ($6 == 2) && ($4 > 0.85)' | shuf | head -n 10 | cut -f 1 > $sv_name.group_2.txt
cat $wild_gwas_git/MDS_outliers/Ha412HO/$species/Ha412HO_inv.$inversion_version.pcasites.$chr.$mds.genotypes.txt | awk '($2 != "NA") && ($6 == 0) && ($4 < 0.15)' | shuf | head -n 10 | cut -f 1 > $sv_name.group_0.txt

cat $sv_name.group_0.txt $sv_name.group_2.txt ann_list > $sv_name.bcftools.list
awk '{print $1"\t0"}' $sv_name.group_0.txt >> $sv_name.sample_info.tsv
awk '{print $1"\t2"}' $sv_name.group_2.txt >> $sv_name.sample_info.tsv
awk '{print $1"\tann"}' ann_list >> $sv_name.sample_info.tsv

bcftools view -r Ha412HOChr06:130156472-155128857 -S $sv_name.bcftools.list --force-samples -O z -o $sv_name.vcf.gz /moonriseNFS/wild_gwas/Annuus.tranche90.snp.annarg.90.bi.remappedHa412HO.vcf.gz
zcat $sv_name.vcf.gz | perl vcf2fasta.pl > $sv_name.fasta
~/programs/iqtree-1.6.10-Linux/bin/iqtree -s $sv_name.fasta -st DNA -m GTR+ASC -nt 20 -bb 1000 -pre $sv_name