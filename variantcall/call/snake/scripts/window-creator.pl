#!/bin/perl
use warnings;
use strict;
#This script lets you put in a number and it will pick a genomic window for printing vcf windows in GATK. Each window will not overlap contigs.

my $n = $ARGV[0] -1; #The window number
my $window_size = 1000000;
my %sizes;
my @chromosomes;
#Pipe in the header lines of the reference genome.
while(<STDIN>){
  chomp;
  my @a = split(/ /,$_);
  my $chrom = $a[0];
  $chrom =~ s/>//;
  my $len = $a[1];
  $len =~ s/len=//;
  $sizes{$chrom} = $len;
  push(@chromosomes,$chrom);
}
my $current_chrom = shift(@chromosomes);
my $window_start = 1;
my $window_end = $window_start + $window_size - 1;
if ($n){
  foreach my $i (1..$n){
    $window_start += $window_size;
    $window_end = $window_start + $window_size - 1;
    if ($window_end > $sizes{$current_chrom}){
      $window_end = $sizes{$current_chrom};
    }
    if ($window_start > $sizes{$current_chrom}){
      #Move to next chromosome;
      unless (@chromosomes){
	      die("Ran out of chromosomes.");
      }
      $current_chrom = shift(@chromosomes);
      $window_start = 1;
      $window_end = $window_start + $window_size - 1;
      if ($window_end > $sizes{$current_chrom}){
	 $window_end = $sizes{$current_chrom};
      }
    }
  }
}
print "$current_chrom:$window_start-$window_end";

