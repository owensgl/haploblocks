sample=$1
name=$(echo $sample | cut -f 5-6 -d "_")
#homer/bin/analyzeHiC $sample -res 1000000 -simpleNorm > $sample.1MB.simple.matrix.txt
cat $sample.1MB.simple.matrix.txt | perl hicmatrix2table.pl $name | pigz -p 5 >  $sample.1MB.simple.table.txt.gz 
