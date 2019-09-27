## This script takes a gene list, generates XML input and run coalescence with BEAST, then pulls out MCC tree. 
## cedar:/home/hkchi/scratch/BEAST/Partition/
## run as: bash beast_run.sh <lineNum>
## Kaichi Huang 2019 Apr

inv_n=$1 # line number from the inversion sheet, 2-19

working_directory="/home/hkchi/scratch/BEAST/Partition"
wild_gwas_git="/scratch/gowens/wild_gwas/wild_gwas_2018/"
inversion_version="jan09" # {jan09, mar27, v3}
inv_page="$wild_gwas_git/MDS_outliers/Ha412HO/all/Ha412HO_inv.$inversion_version.txt"
species=$(cat $inv_page | sed -n ${inv_n}p | cut -f 1)
chr_n=$(cat $inv_page | sed -n ${inv_n}p | cut -f 4)
chr_padded=$(printf %02d $chr_n)
chr="Ha412HOChr$chr_padded"
direction=$(cat $inv_page | sed -n ${inv_n}p | cut -f 6)
mds_n=$(cat $inv_page | sed -n ${inv_n}p | cut -f 5)
mds=${direction}$mds_n

cd "$working_directory/$species/$chr.${mds}"

n_filter_gene=`wc -l inv.phylo_genes.list | cut -f 1 -d " "`
if [ "$n_filter_gene" -ge 100 ]
then
	n_filter_gene=100
fi

# Generate XMLs
cat inv.phylo_genes.list | shuf | head -n $n_filter_gene | sed 's/^/inv.&/g' | sed 's/$/&.fasta/g' | perl $working_directory/fastas2BEAST1xml.pl $chr $mds inv > $working_directory/$species/$chr.$mds.inv.xml

# Run BEASTv1.10.4 on a 48-core node
cd ..
module load beagle-lib/3.1.1
export LD_LIBRARY_PATH=/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/beagle-lib/3.1.1/lib:$LD_LIBRARY_PATH
~/program/BEASTv1.10.4/bin/beast -threads 48 $working_directory/$species/$chr.$mds.inv.xml

# Assuming convergence after 10M runs, generate MCC tree with TreeAnnotator
~/program/BEASTv1.10.4/bin/treeannotator -burnin 100000 -heights median $working_directory/$species/$chr.$mds.inv.trees $working_directory/$species/$chr.$mds.inv.mcc.tre

# Download and visualize TREs manually in FigTree
