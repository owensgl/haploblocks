#!/usr/bin/perl

use warnings;
use strict;

unless (@ARGV == 3) {die;}

my $in = $ARGV[0]; #merged bsnp table loci (rows) pop (columns)
my $pop = $ARGV[1];# list of population with ind \t pop
my $out = $ARGV[2]; #name of out file
my @pop=();
my $loc;
my %pop;
my $c=1;
my %h;
my %counter;
my %samples;
my @samples;
my @loc;
my $locNum;
my @a;
my $pop_count;
my %pop_no;
my @ind;
my $cn=1;
my $num=1;
open POP, $pop;
while (<POP>){
	chomp;
	unless (/^\s*$/ | /^ID/){
		@a = split;
		$pop{$cn}=$a[1];
		$pop_no{$a[1]}=1;
		$cn++;
}
}
close POP;
open OUT3, ">$out.pop_order";
my $no=1;
foreach my $p (sort keys %pop_no){
		print OUT3 "$no $p\n";
		push(@pop, $p);
		$no++;
}

foreach my $p (sort keys %pop){
		#print "pop $p $pop{$p}\n";
}

my $pop_no=scalar keys %pop_no;
print "there are $pop_no populations\n";


open IN, $in;
open OUT, ">$out";
open OUT2, ">$out.loci";

#read in snp table
while (<IN>){
	chomp; #delete hard returns
	if(/^#/){ #skip any header
		next;
	}	
	@a = split;
	$locNum++;
	my @l=();
	my $l=shift (@a);
	my $l2=shift(@a);	
	push(@l, $l, $l2);
    $loc = join'__', @l; #name with genomic contig and genomic position
	$pop_count=1;		
	foreach my $i (@a){	
		if($i =~ /A|T|G|C/ ){		
#			print "population $pop{$pop_count}\n";
			
			$counter{$loc}{$pop{$pop_count}}{$i}++; #hash of hash of a hash locus -> pop -> allele count
			$counter{$loc}{$pop{$pop_count}}{$i}++;			
		}if($i =~ /K|W|Y/){
			$counter{$loc}{$pop{$pop_count}}{'T'}++;
		}if($i =~ /K|R|S/){
			$counter{$loc}{$pop{$pop_count}}{'G'}++;
		}if($i =~ /R|W|M/){
			$counter{$loc}{$pop{$pop_count}}{'A'}++;
		}if ($i =~ /M|S|Y/){
			$counter{$loc}{$pop{$pop_count}}{'C'}++;
		}
		$pop_count=$pop_count+1; 
	}
		foreach my $k (keys %counter){
				my %alleles=();
				for (my $i=0; $i < @pop  ; $i++){
						foreach my $k2 (keys %{$counter{$k}{$pop[$i]}}){		
								#print "$k2\n";
								$alleles{$k2}=1;#get the alleles for each loci
					}
				}
				foreach my $key (keys %alleles){
					#	print "alleles $key\n";
				}
	

		for (my $j=0; $j < @pop; $j++){
			foreach my $key (keys %alleles){
				unless(exists $counter{$k}{$pop[$j]}{$key}){
					$counter{$k}{$pop[$j]}{$key}=0;#insert 0 for each population where that allele does not exist
					#print "MISS $k $key $counter{$k}{$pop[$j]}{$key}\n";
				}		
			}	
		}		
	
	}

                my @temp1=();
                my $allele_count=();
                my @temp=();
		#print out allele count and number in each population
		foreach my $key (sort keys %counter){
#		print OUT "\n";
				$num++;
				@temp1=();
				$allele_count=();
				@temp=();
				foreach my $key2 (sort keys %{$counter{$key}}){
		#				print "pop counter $key2 \t";
						$c=1;		
						foreach my $key3 (sort keys %{$counter{$key}{$key2}}){
						#		print "counter $counter{$key}{$key2}{$key3}\t allele key $key3 pop $key2\n";	
								$allele_count=scalar keys %{$counter{$key}{$key2}};
								if($c==1){
										push(@temp, $counter{$key}{$key2}{$key3});
								}if($c==2){
										push(@temp1, $counter{$key}{$key2}{$key3});
								}
									$c++;
						}
				}
				if($allele_count == 2){#do not print if there are more than 2 alleles
						print OUT2 "$num\t$key\n";#print out the loci order for those that have 2 alleles
						foreach (@temp){
						print OUT "$_\t";
						}
						print OUT "\n";
						foreach (@temp1){
								print OUT "$_\t";
						}
						print OUT "\n";
				}				
	
	
		}
%counter=();
}
