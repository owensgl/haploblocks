
#!/usr/bin/perl
use warnings;
use strict;

#this script will filter a snp table based on the minor allele frequency, the number of individuals with genotypes called, the observed frequency of heterozygotes

#the SNP table
my $list = $ARGV[0];
open LIST, $list;

#summary of the SNP allele frequencies
open OUT, ">$list.summary";
#filtered SNP table
open OUT2, ">$list.table";

my $allele_cut=0.03; #minor allele frequency minimum
my $geno_cut = 5; #miminum number of individuals with genotype information
my $het_cut = .7; #max frequency of heterozygotes

my $sites=0; #count up the total number of variable sites

#read in the calls
while (<LIST>){

	#allele counts
	my $a_count=();
	my $t_count=();
	my $c_count=();
	my $g_count=();

	#allele total
	my $alleles=();

	#genotype counts
	my $aa_count=();
	my $tt_count=();
	my $cc_count=();
	my $gg_count=();
	my $nn_count=();
	my $ag_count=();
	my $at_count=();
	my $ac_count=();
	my $cg_count=();
	my $ct_count=();
	my $gt_count=();
	my $count=();
	
	#total number of genotypes
	my $geno=();
	
	#frequency of heterozygotes
	my $het=();

	my $join=();
	my @split=();
	if (/#/){
		next;
	}
	if (/^$/){
		next;
	}
	chomp;
	
	@split=split;
	my @temp=();
	push (@temp, $split[0], $split[1]);
	my $position=();
	$position=join "\t", @temp;
	#remove the contig and position information
	shift(@split);
	shift(@split); 
        $join=join"", @split;
        #@split=split "", $join;

	#get the number of each genotype
        my $A=$join;
        $A =~ s/T|C|G|N|X|-|K|R|W|M|S|Y//g;
        my $T=$join;
        $T =~ s/A|C|G|N|X|-|K|R|W|M|S|Y//g;
        my $C=$join;
        $C =~ s/T|A|G|N|X|-|K|R|W|M|S|Y//g;
        my $G=$join;       
		$G =~ s/A|T|C|N|X|-|K|R|W|M|S|Y//g;
        my $N=$join;       
		$N =~ s/A|T|C|G|X|-|K|R|W|M|S|Y//g;
        my $R=$join;
        $R =~ s/A|T|C|G|X|-|K|N|W|M|S|Y//g;
        my $W=$join;
        $W =~ s/A|T|C|G|X|-|K|N|R|M|S|Y//g;
        my $M=$join;
        $M =~ s/A|T|C|G|X|-|K|N|W|R|S|Y//g;
        my $S=$join;
        $S =~ s/A|T|C|G|X|-|K|N|W|M|R|Y//g;
        my $Y=$join;
        $Y =~ s/A|T|C|G|X|-|K|N|W|M|S|R//g;
        my $K=$join;
        $K =~ s/A|T|C|G|X|-|N|Y|W|M|S|R//g;
       	my $test=$join;
		$test =~ s/A|T|C|G|X|-|N|Y|W|M|S|R|K//g; 
	#make sure there are no other characters in the table
	if (length($test)>0){
		print "error unrecognized character $test in contig $position\n";
		next;
	}

		$aa_count=length($A);
        $tt_count=length($T);
        $cc_count=length($C);
        $gg_count=length($G);
        $nn_count=length($N);
        $ag_count=length($R);
        $at_count=length($W);
        $ac_count=length($M);
        $cg_count=length($S);
        $ct_count=length($Y);
        $gt_count=length($K);


	#total the number of genotypes
	$geno=$aa_count+$tt_count+$cc_count+$gg_count+$ag_count+$at_count+$ac_count+$cg_count+$ct_count+$gt_count;
	my $homo=$aa_count+$tt_count+$cc_count+$gg_count;

	#the frequncy of heterozygotes
	if($geno > 0){
		$het=($ag_count+$at_count+$ac_count+$cg_count+$ct_count+$gt_count)/$geno;
	}

	#array of genotype counts
	my @geno_count=();
        push @geno_count, $aa_count, $tt_count, $cc_count, $gg_count, $ag_count, $at_count, $ac_count, $cg_count, $ct_count, $gt_count;
        my %hash=();

	#count the number of alleles
        $a_count=2*$aa_count+$ag_count+$at_count+$ac_count;
        $t_count=2*$tt_count+$at_count+$ct_count+$gt_count;
        $c_count=2*$cc_count+$ac_count+$cg_count+$ct_count;
        $g_count=2*$gg_count+$ag_count+$cg_count+$gt_count;

	#count the total number of alleles
	$alleles=$a_count+$t_count+$c_count+$g_count;

	#push onto allele hash key = base value = count
	$hash{"A"}=$a_count;
	$hash{"T"}=$t_count;
	$hash{"C"}=$c_count;
	$hash{"G"}=$g_count;
	my $key;

	#sort allele counts
	my @sorted = sort {$hash{$b} <=> $hash{$a} }keys %hash;
	my @top= sort { $b <=> $a} @geno_count;

	#determine if locus is variable
	my $var=();
	if ($homo == ($top[0]/2)){
			$var=1;
    }else{$var=0;}
	#print out each variable site, removing high frequency heterozygotes and SNPs that have too few individuals with genotypes called
		if ($var == 0  && $het <= $het_cut && $geno >= $geno_cut) {
			if ($hash{$sorted[0]} < $alleles){
				#determine major and minor allele frequencies
				my $mj=$hash{$sorted[0]}/$alleles;
				my $mn=$hash{$sorted[1]}/$alleles;
				@temp=();
				my $join= join "\t", @temp;
				#print the summary information for each locus if the minor allele frequency is greater than the threshold and the SNP calls to a table
				if ($mn > $allele_cut){
					print OUT "$position\t$alleles\t$geno\t$top[0]\t$sorted[0]\t$hash{$sorted[0]}\t$mj\t$sorted[1]\t$hash{$sorted[1]}\t$mn\t$het\n";
					print OUT2 "$_\n";
				}
			}
		}

	#count up the total number of sites
	if($alleles>1){
		$sites=$sites+1 ;
	}	
}

print "total number of sites $sites \n";
