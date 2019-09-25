#!/bin/perl;
use warnings;
use strict;

my $first_line;
my $line_counter =1;
while(<STDIN>){
	chomp;
	my @a = split(/\t/,$_);
	my %hash;
	foreach my $i (4..$#a){
		$hash{$a[$i]}++;
	}
	if ($first_line){
		unless($hash{2}){
			$hash{2} = 0;
		}
		unless($hash{1}){
			$hash{1} = 0;
		}
		my $freq = $hash{2}/($hash{1}+$hash{2});
		print "\n$a[0]\t$a[3]\t$freq"; 
	}else{
                unless($hash{2}){
                        $hash{2} = 0;
                }
                unless($hash{1}){
                        $hash{1} = 0;
                }
		$first_line++;
		my $freq = $hash{2}/($hash{1}+$hash{2});
		print "$a[0]\t$a[3]\t$freq";
	}
	$line_counter++;
	if ($line_counter % 100000 == 0){
		print STDERR "Printing $a[0] $a[3]...\n";
	}
}
