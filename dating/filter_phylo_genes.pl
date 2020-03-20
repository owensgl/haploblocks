#!/bin/perl
# Filter gene fasta based on averge missing rate across samples, and generate a gene list. 
# run as: perl filter_phylo_genes.pl inv inv.phylo_genes.list
## Kaichi Huang 2019 Apr

use List::Util qw(sum);  
use List::MoreUtils qw(none);

my $type = $ARGV[0];
my $n = `ls $type.gene_*.fasta | wc -l | cut -f 1 -d " "`;

my @gene_list;

foreach my $i (1..$n) {
	my @missing;
	open GENE, "$type.gene_$i.fasta";
	while(<GENE>){
		chomp;
		if(/^>/){
			next;
		} else {
			my $seq = $_;
			my $count = ($seq=~tr/N/N/);
			my $missing = $count / length($seq);
			push @missing, $missing;
		}
	}
	close GENE;
	if ((sum(@missing)/($#missing+1) < 0.3) or (none {$_>0.5} @missing)) {
		push @gene_list, "gene_$i";
	}
}

open OUT, ">$ARGV[1]"; 
foreach my $gene (@gene_list) {
	print OUT "$gene\n";
}
close OUT;