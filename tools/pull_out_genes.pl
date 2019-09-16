#!/bin/perl
use warnings;
use strict;
#Picks genes in regions
my $chr = $ARGV[0];
my $start = $ARGV[1];
my $end = $ARGV[2];

while(<STDIN>){
  chomp;
  if ($_ =~ m/^#/){next;}
  my @a = split(/\t/,$_);
  if ($a[2] ne "gene"){next;}
  if (($a[3] >= $start) and ($a[4] <= $end) and ($a[0] eq $chr)){
    my @fields = split(/;/,$a[8]);
    my $gene = $fields[0];
    $gene=~ s/ID=gene://g;
    print "$gene\n";
  }
}
