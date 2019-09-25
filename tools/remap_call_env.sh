module load bcftools
module load tabix
module load bwa
module load samtools

species="Annuus"
prefix="tranche90.snp.env.90.bi"
#Merge all chromosomes including 00c
if [ ! -f $species.$prefix.tmp.vcf.gz ]
then
	bcftools concat --threads 10 -O z -a Han*.$prefix.vcf.gz > $species.$prefix.tmp.vcf.gz
fi

#Remap to Ha412HO
if [ ! -f $species.$prefix.remappedHa412.vcf.gz ]
then
	perl /home/gowens/bin/xrqpos2ha412pos_bwa.pl $species.$prefix.tmp.vcf.gz $species.$prefix
	tabix -p vcf $species.$prefix.remappedHa412.vcf.gz
fi

#Remove extra contigs
if [ ! -f $species.$prefix.remappedHa412HO.vcf.gz ]
then
	bcftools view --threads 10 -R ../Ha412HO.chrlist.bed -O z $species.$prefix.remappedHa412.vcf.gz  > $species.$prefix.remappedHa412HO.vcf.gz
	tabix -p vcf $species.$prefix.remappedHa412HO.vcf.gz
fi
#Phase and impute
if [ ! -f $species.$prefix.remappedHa412HO.beagle.vcf.gz ]
then
	module load java

	java -Xmx120g -jar /home/gowens/bin/beagle.10Jun18.811.jar gt=$species.$prefix.remappedHa412HO.vcf.gz out=$species.$prefix.remappedHa412HO.beagle impute=true nthreads=30
	tabix -p vcf $species.$prefix.remappedHa412HO.beagle.vcf.gz
fi
