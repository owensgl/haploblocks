## This script takes an inversion, pulls out inversion regions in Ha412, gets the genes in that region, then pulls out orthologous genes in XRQ. 
## cedar:/home/hkchi/scratch/BEAST/Partition/
## run as: bash run_inv_phylogeny.sh <lineNum>
## modified from G. Owens
## Kaichi Huang 2019 Apr

inv_n=$1 # line number from the inversion sheet, 2-19

working_directory="/home/hkchi/scratch/BEAST/Partition"
wild_gwas_git="/scratch/gowens/wild_gwas/wild_gwas_2018/"
inversion_version="jan09" # {jan09, mar27, v3}
outgroups="GRO_2043,DIV_1956"
xrq_ref="/scratch/gowens/ref/HanXRQr1.0-20151230.fa"

cd $working_directory
inv_page="$wild_gwas_git/MDS_outliers/Ha412HO/all/Ha412HO_inv.$inversion_version.txt"
species=$(cat $inv_page | sed -n ${inv_n}p | cut -f 1)
chr_n=$(cat $inv_page | sed -n ${inv_n}p | cut -f 4)
chr_padded=$(printf %02d $chr_n)
chr="Ha412HOChr$chr_padded"
direction=$(cat $inv_page | sed -n ${inv_n}p | cut -f 6)
mds_n=$(cat $inv_page | sed -n ${inv_n}p | cut -f 5)
mds=${direction}$mds_n
Ha412_genes="$wild_gwas_git/resources/Han412-HO_gene_filt.gff3.gz"

#Samples from other species are pre selected and stored in these files
other_samples="$wild_gwas_git/resources/$species.physamples.txt"

mkdir -p $species
cd $species
mkdir -p "$chr.${mds}"
cd "$chr.${mds}"

if [ ! -f group_2.txt ]
then
#Get list of samples appropriate for each inversion test
cat $wild_gwas_git/MDS_outliers/Ha412HO/$species/Ha412HO_inv.$inversion_version.pcasites.$chr.$mds.genotypes.txt | awk '($2 != "NA") && ($6 == 2) && ($4 > 0.85)' | shuf | head -n 5 | cut -f 1 > group_2.txt
cat $wild_gwas_git/MDS_outliers/Ha412HO/$species/Ha412HO_inv.$inversion_version.pcasites.$chr.$mds.genotypes.txt | awk '($2 != "NA") && ($6 == 0) && ($4 < 0.15)' | shuf | head -n 5 | cut -f 1 > group_0.txt
cat group_0.txt group_2.txt $other_samples > $chr.$mds.samplelist.txt

#Get HA412 genes in inversion region
cat $wild_gwas_git/MDS_outliers/Ha412HO/all/Ha412HO_inv.$inversion_version.inversions.regions.v1.txt | grep $species | grep $chr | grep $mds | cut -f 2-4 > inversion.regions.txt
n_regions=$(wc -l inversion.regions.txt | cut -f 1 -d " ")
rm inv.ha412.genes.txt
for n in `seq $n_regions`
do
	region_chr=$(cat inversion.regions.txt | sed -n ${n}p | cut -f 3)
	region_start=$(cat inversion.regions.txt | sed -n ${n}p | cut -f 1)
	region_end=$(cat inversion.regions.txt | sed -n ${n}p | cut -f 2)
	zcat $Ha412_genes | perl $wild_gwas_git/perl_tools/pull_out_genes.pl $region_chr $region_start $region_end >> inv.ha412.genes.txt
done

#Get orthologous XRQ genes
cat $wild_gwas_git/resources/Ha412HO_HanXRQv1.1to1ortho.txt | grep -f inv.ha412.genes.txt | cut -f 2 | tr -d ' '  > inv.xrq.genes.txt
cat $wild_gwas_git/resources/HanXRQr1.0-20151230-EGN-r1.1.genelocations.txt | grep -f inv.xrq.genes.txt > inv.xrq.genes.locations.txt
fi

#Make all the fasta files
# most time-consuming step
inv_genes="inv.xrq.genes.locations.txt"
module load gatk
n_genes=$(wc -l $inv_genes | cut -f 1 -d " ")
for i in `seq $n_genes`
do
	if [ ! -s "inv.gene_$i.fasta" ]
	then
		cat $inv_genes | sed -n ${i}p | cut -f 2- > gene.bed
		cat $chr.$mds.samplelist.txt | parallel -j 20 bash $wild_gwas_git/perl_tools/make_fasta.sh $wild_gwas_git/resources/sample_info_file_all_samples_2018_12_07.tsv $xrq_ref $wild_gwas_git/perl_tools
		for f in *.fa; do (cat "${f}"; echo); done > inv.gene_$i.fasta
		rm *.g.vcf*
		rm *.fa
	fi
done

# Filter by missing sites of each gene and generate a list of good genes for phylogeny
module load perl
perl $working_directory/filter_phylo_genes.pl inv inv.phylo_genes.list
