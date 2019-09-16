#!/bin/perl
use strict;
use warnings;

#This script takes a VCF, a list of outgroup samples, and outputs only sites where the outgroups are fixed for a single allele. It then polarizes all SNPs as A (ancestral) or B (derived). Requires biallelic sites only.
#In this version the outgroup is defined in a separate file

my $outfile = $ARGV[0]; #e.g. perennial_alleles.20190619.txt.gz
my $vcffile = $ARGV[1]; #The VCF to be processed
#First run through the VCF to load all sites that need processing

my %used_sites;
my $linecounter = 1;
open(VCF, "gunzip -c $vcffile |");
while(<VCF>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#/){
    next;
  }
  my $info = $a[7];
  my @infos = split(/\;/,$info);
  my $XRQ_info = $infos[$#infos];
  $XRQ_info =~ s/XRQ=//g;
  my $location = $XRQ_info;
  $used_sites{"$location"}++;
  if ($linecounter % 1000000 == 0){print STDERR "Initial vcf processing $location...\n";}
}
close VCF;
print STDERR "Loaded all sites\n";
my %ancestral_state;
open(OUTGROUP, "gunzip -c $outfile |");
$linecounter = 1;
while(<OUTGROUP>){
  chomp;
  if ($. == 1){
    next;
  }
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $ref = $a[2];
  my $count = $a[3];
  if ($count < 2){next;}
  $linecounter++;
  if ($linecounter % 1000000 == 0){print STDERR "Processing outgroup file $chr $pos...\n";}
  unless ($used_sites{"$chr.$pos"}){next;}
  $ancestral_state{"$chr.$pos"} = $ref;
}
close OUTGROUP;
print STDERR "Finished processing outgroup file\n";
open(VCF2, "gunzip -c $vcffile |");
my %sample;
while(<VCF2>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#CHR/){
    print "chr\tpos";
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
      print "\t$a[$i]";
    }
    next;
  }
  if ($_ =~ m/^#/){next;}
  my $info = $a[7];
  my @infos = split(/\;/,$info);
  my $XRQ_info = $infos[$#infos];
  $XRQ_info =~ s/XRQ=//g;
  my $location = $XRQ_info;
  my $ref = $a[3];
  my $alts = $a[4];
  unless ($ancestral_state{"$location"}){print STDERR "Ancestral state not found\t$location\n";next;}
  unless (($ancestral_state{"$location"} eq $ref) or (($ancestral_state{"$location"} eq $alts))){print STDERR "Ancestral allele not found\n";next;}
  my $outgroup_allele;
  if ($ref eq $ancestral_state{"$location"}){
    $outgroup_allele = 0;
  }else{
    $outgroup_allele = 1;
  }

  my @alts = split(/,/,$alts);
  if ($alts[1]){next;} #Skips sites with more than two alleles
  
  print "\n$a[0]\t$a[1]";
  #Now to figure out if the samples have derived or ancestral
  foreach my $i (9..$#a){
    if (($a[$i] eq '.') or ($a[$i] eq '.:0,0')){
      print "\tN";
    }else{
      my @infos = split(/:/,$a[$i]);
      if ($infos[0] eq './.'){
        print "\tN";
        next;
      }
      my @genos = split(/\//,$infos[0]);
      if (($genos[0] == $outgroup_allele ) and ($genos[1] == $outgroup_allele )){
        print "\tA";
      }elsif (($genos[0] != $outgroup_allele ) and ($genos[1] != $outgroup_allele)){
        print "\tB";
      }else{
        print "\tH";
      }
    }
  }
}
close VCF2;
