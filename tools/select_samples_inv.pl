#!/bin/perl
use strict;
use warnings;
use List::Util 'shuffle';

my $species_used = $ARGV[0];
my $genotype_file = $ARGV[1];
my $samplelist = "/home/owens/bin/wild_gwas_2018/resources/sample_info_file_all_samples_wildspecies.tsv";

my %sample_geno;
open LIST, $samplelist;

my %ok_list;
while(<LIST>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $name = $a[2];
  my $species = $a[3];
  if ($species eq $species_used){
    $ok_list{$name}++;
  }
}
close LIST;

my %backup_list;
open GENO, $genotype_file;
while(<GENO>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $name = $a[0];
  my $geno = $a[5];
  if ($geno == 1){next;}
  if (($a[3] > 0.1) and ($a[3] < 0.9)){
    if ($ok_list{$name}){
      $backup_list{$name} = $geno;
      next;
    }
  }
  if ($ok_list{$name}){
    $sample_geno{$name} = $geno;
  }
}
close GENO;

my $zero_print = 0;
my $two_print = 0;
foreach my $name ( shuffle keys %sample_geno ){
  if (($sample_geno{$name} == 0) and ($zero_print < 5)){
    print "$name\t0\t$genotype_file\n";
    $zero_print++;
  }
  if (($sample_geno{$name} == 2) and ($two_print < 5)){
    print "$name\t2\t$genotype_file\n";
    $two_print++;
  }
}
#Backup files that are admixed if no pure types are available
foreach my $name ( shuffle keys %backup_list ){
  if (($backup_list{$name} == 0) and ($zero_print < 5)){
    print "$name\t0\t$genotype_file\n";
    $zero_print++;
  }
  if (($backup_list{$name} == 2) and ($two_print < 5)){
    print "$name\t2\t$genotype_file\n";
    $two_print++;
  }
}





if ($zero_print < 5){
  print STDERR "$species_used\t$genotype_file\tprinted $zero_print/5 for 0\n"
}
if ($two_print < 5){
  print STDERR "$species_used\t$genotype_file\tprinted $two_print/5 for 2\n"
}
