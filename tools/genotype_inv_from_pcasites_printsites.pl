#!/bin/perl
use strict;
use warnings;

#This script takes the list of MDS PCA outliers and  records what alleles are for each group. It then uses a vcf file to record how many of each type of allele each sample has. 
my $mds_fst_file = $ARGV[0]; #List of fst outliers
my $gzvcf = $ARGV[1];
#my $min_freq = 0.95;
my $Ha412_remap = "TRUE"; #Set to TRUE to output the remapped positions

open (MDS, "gunzip -c $mds_fst_file |");
my %allele_1;
my %allele_2;
my %ha412_chr;
my %ha412_pos;
while(<MDS>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[2];
  my $pos = $a[3];
  my $ha412_chr = $a[0];
  my $ha412_pos = $a[1];
  $ha412_chr{$chr.$pos}= $ha412_chr;
  $ha412_pos{$chr.$pos}= $ha412_pos;
  my $inv_0_allele = $a[7];
  my $inv_2_allele = $a[8];
  $allele_1{$chr}{$pos} = $inv_0_allele;
  $allele_2{$chr}{$pos} = $inv_2_allele;
}
close MDS;
open(IN, "gunzip -c $gzvcf |");

my %sample;
my %counts;
my $site_count = 0;
print "XRQchr\tXRQpos\tchr\tpos\tsample\tgenotype\tdepth";
while(<IN>){
  chomp;
  if ($_ =~ m/^##/){
    next;
  }
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#/){
    print STDERR "Loading vcf file now...\n";
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
    }
    next;
  }
  my $chr = $a[0];
  my $pos = $a[1];
  $site_count++;
  if ($site_count % 100000 == 0){print STDERR "$chr $pos processed...\n";}
  unless ($allele_1{$chr}{$pos}){next;}
  my $ref = $a[3];
  my $alt = $a[4];
  my @alts = split(/,/,$alt);
    #Check to make sure that the alleles match;
    #Find out which allele in this site matches with the inversion alleles
    my $allele_1 = "NA";
    my $allele_2 = "NA";
    if ($ref eq $allele_1{$chr}{$pos}){
      $allele_1 = 0;
    }elsif ($ref eq $allele_2{$chr}{$pos}){
      $allele_2 = 0;
    }
    foreach my $x (0..$#alts){
       if ($alts[$x] eq $allele_1{$chr}{$pos}){
         $allele_1 = ($x+1);
       }elsif ($alts[$x] eq $allele_2{$chr}{$pos}){
         $allele_2 = ($x+1);
       }
     }
    my @info = split(/:/,$a[8]);
    unless($info[1]){next;}
    my $dp_n;
    if ($info[1] eq "DP"){
      $dp_n = 1;
    }elsif ($info[2] eq "DP"){
      $dp_n = 2;
    }else{
      print STDERR "No DP, skipping line\n";
      next;
    }
    foreach my $i (9..$#a){
      if (($a[$i] eq './.') or ($a[$i] eq '.')){next;}
      my @fields = split(/:/,$a[$i]);
      my $genotype = $fields[0];
      my $dp = $fields[$dp_n];
      unless($dp){next;}
      if ($dp == 0){next;}
      if (($genotype eq '.') or ($genotype eq './.')){next;}
      if (($genotype eq "$allele_1\/$allele_2") or  ($genotype eq "$allele_2\/$allele_1")){
	$counts{$sample{$i}}{1}++;
        print "\n$chr\t$pos\t$ha412_chr{$chr.$pos}\t$ha412_pos{$chr.$pos}\t$sample{$i}\t01\t$dp";
      }elsif ($genotype eq "$allele_1\/$allele_1"){
	print "\n$chr\t$pos\t$ha412_chr{$chr.$pos}\t$ha412_pos{$chr.$pos}\t$sample{$i}\t00\t$dp";      
      }elsif ($genotype eq "$allele_2\/$allele_2"){
	print "\n$chr\t$pos\t$ha412_chr{$chr.$pos}\t$ha412_pos{$chr.$pos}\t$sample{$i}\t11\t$dp";
      }
    
  }
}  
