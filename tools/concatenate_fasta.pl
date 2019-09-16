#!/bin/perl
use warnings;
use strict;
#Concatenate fasta and make partition file

my $partion_name = "partions.txt";
my %samples;
my %seqs;
open(my $partition_file, '>', $partion_name);
my $current_partion_start = 0;
#Feed it a list of fasta files to concatenate
while(<STDIN>){
  chomp;
  my $file = $_;
  open FILE, $file;
  my $seq_length = 0;
  my $sample;
  while(<FILE>){
    chomp;
    next if /^\s*$/;
    if ($_ =~ m/>/){
      $sample = $_;
      $samples{$sample}++;
    }else{
      $seqs{$sample} .= $_;
      $seq_length = length($_);
    }
  }
  close FILE;
  if ($seq_length == 0){next;}
  my $partion_start = $current_partion_start + 1;
  my $partion_end = $current_partion_start + $seq_length;
  print $partition_file "DNA, $file = $partion_start-$partion_end\n";
  $current_partion_start = $partion_end;
}
close $partition_file;
foreach my $sample (sort keys %samples){
  print "$sample\n$seqs{$sample}\n";
}
