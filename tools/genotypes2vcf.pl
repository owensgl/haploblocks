#!/bin/perl
use warnings;
use strict;
#This script creates VCFs for inversions from .genotypes.txt files with triangle genotypes. 

my $info_file = "/home/owens/bin/wild_gwas_2018/sample_info_apr_2018.tsv";
my $date_chosen = "v3";
my $species_name = "argophyllus";
my $species_abbreviation = "Arg";
my $mdsoutlier_directory="/home/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/".$species_name;
my $vcf_header="##fileformat=VCFv4.2\n##source=Ha412HO_inv.$date_chosen\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
my $vcf_sample_line="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

open INFO, $info_file;
my %species;
my @samples;
while(<INFO>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $sample = $a[2];
  my $species = $a[3];
  $species{$sample} = $species;
  if ($species =~ m/$species_abbreviation/){
    $vcf_sample_line .= "\t$sample";
    push(@samples, $sample);
  }
}

my @files = `ls $mdsoutlier_directory | grep genotypes.txt`;
my $out_file = $mdsoutlier_directory."/Ha412HO_inv.".$date_chosen.".pcasites.vcf";
open my $vcf_out, '>',$out_file;

print $vcf_out "$vcf_header\n$vcf_sample_line";
my %genotypes;
my %mds_numbers;
print STDERR "@files\n";
foreach my $file (@files){
  my @file_infos = split(/\./,$file);
  my $chr = $file_infos[3];
  my $mds = $file_infos[4];
  my $mds_n = $file_infos[4];
  $mds_n =~ s/neg//g;
  $mds_n =~ s/pos//g;
  $mds_n =~ s/syn//g;
  $mds_numbers{$chr}{$mds_n} = $mds;
  open (GENOTYPES, "cat $mdsoutlier_directory/$file | ");
  while(<GENOTYPES>){
    chomp;
    if ($. == 1){next;}
    my @a = split(/\t/,$_);
    my $sample = $a[0];
    my $genotype = $a[5]; #Triangle genotype
    unless ($species{$sample}){next;}
    if ($species{$sample} =~ m/$species_abbreviation/){
      $genotypes{$chr}{$mds_n}{$sample} = $genotype;
    }
  }
}
foreach my $chr (sort keys %genotypes){
  foreach my $mds (sort keys %{$genotypes{$chr}}){
    print $vcf_out "\n$chr\t$mds\t$mds_numbers{$chr}{$mds}\t";
    print $vcf_out "A\tT\t.\t.\t.\tGT";
    foreach my $sample (@samples){
      unless(exists($genotypes{$chr}{$mds}{$sample})){
        print $vcf_out "\t./.";
	next;
      }
      if ($genotypes{$chr}{$mds}{$sample} == 0){
        print $vcf_out "\t0/0";
      }elsif ($genotypes{$chr}{$mds}{$sample} == 1){
        print $vcf_out "\t0/1";
      }elsif ($genotypes{$chr}{$mds}{$sample} == 2){
        print $vcf_out "\t1/1";
      }
    }
  }
}

