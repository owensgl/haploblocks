species="petfal";
capital_species="Petiolaris"
prefix="${capital_species}.tranche90.snp.petfal.90.bi.remappedHa412HO"
INV_file="Ha412HO_inv.v3.pcasites.${species}.shuffled"
#cov_file="Annuus.tranche99.snp.gwas.90.bi.remappedHa412HO.jan09noinv.ldr0p2.4pcas.cov"
for phenotype in `ls /home/owens/bin/wild_gwas_2018/gwas/$species | grep -v znorm.txt | sed s/.txt//g `
do
	if [ -f $species/results/${INV_file}.${phenotype}.ps.gz ]
	then
		continue
	fi
        echo "Running $species with FULL GENOME, without inversion kinf..."
        /home/owens/bin/emmax-intel64 -v -d 10 -t $species/${prefix}.beagle -p /home/owens/bin/wild_gwas_2018/gwas/$species/${phenotype}.txt -k $species/${prefix}.beagle.v3noinv.ldfilter.aIBS.kinf -o $species/results/${prefix}.${phenotype}.beagle.fullgenome.v3noinv.ldfilter -c $species/${prefix}.v3noinv.ldfilter.ldr0p2.PC1-3.cov.txt
        gzip -f $species/results/${prefix}.${phenotype}.beagle.fullgenome.v3noinv.ldfilter.ps
	#Run on snps with full genome PCA and kinship
	echo "Running $species with FULL GENOME..."
	/home/owens/bin/emmax-intel64 -v -d 10 -t $species/${prefix}.beagle -p /home/owens/bin/wild_gwas_2018/gwas/$species/${phenotype}.txt -k $species/${prefix}.beagle.aIBS.kinf -o $species/results/${prefix}.${phenotype}.beagle -c $species/${prefix}.ldr0p2.PC1-3.cov.txt
	#Run on inversions with inversion remove genome PCA and kinship
	echo "Running $species with INVERSIONS..."
	/home/owens/bin/emmax-intel64 -v -d 10 -t $species/$INV_file -p /home/owens/bin/wild_gwas_2018/gwas/$species/${phenotype}.txt -k $species/${prefix}.beagle.v3noinv.ldfilter.aIBS.kinf -o $species/results/${INV_file}.$phenotype -c $species/${prefix}.v3noinv.ldfilter.ldr0p2.PC1-3.cov.txt
	#Run on SNPs with inversion removed genome PCA and kinship
	echo "Running $species with INV CONTROLLED SNPS..."
	/home/owens/bin/emmax-intel64 -v -d 10 -t $species/${prefix}.beagle -p /home/owens/bin/wild_gwas_2018/gwas/$species/${phenotype}.txt -k $species/${prefix}.beagle.v3noinv.ldfilter.aIBS.kinf -o $species/results/${prefix}.${phenotype}.beagle.v3noinv.ldfilter -c $species/${prefix}.v3noinv.ldfilter.ldr0p2.PC1-3.cov.txt
	#Gzip the files
	gzip -f $species/results/${prefix}.${phenotype}.beagle.ps
	gzip -f $species/results/${prefix}.${phenotype}.beagle.v3noinv.ldfilter.ps
	gzip -f $species/results/${INV_file}.${phenotype}.ps
done
