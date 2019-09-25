#!/bin/perl
#Turns a matrix into a table
use strict;
use warnings;

my $interchromosomal = "TRUE"; #Does it print out interchromosomal reactions. TRUE OR FALSE
my %chrs;
my %chr_list;
my %pos;
my %data;
my $genotype = $ARGV[0]; #Genotype column

unless($genotype){
	print STDERR "No genotype given\n";
	exit;
}
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. == 1){
    foreach my $i (2..$#a){
      my @info = split(/-/,$a[$i]);
      $chrs{$i} = $info[0];
      $chr_list{$info[0]}++;
      $pos{$info[0]}{$i} = $info[1];
    }
    next;
  }
  my @info = split(/-/,$a[0]);
  my $row_chr = $info[0];
  my $row_pos = $info[1];
  foreach my $i (2..$#a){
    $data{$row_chr}{$chrs{$i}}{$row_pos}{$pos{$chrs{$i}}{$i}} = $a[$i];
  }
}

my $first_line;
if ($interchromosomal eq "FALSE"){
  foreach my $chr (sort keys %chr_list){
    foreach my $pos1 (sort values %{$pos{$chr}}){
      foreach my $pos2 (sort values %{$pos{$chr}}){
        if (exists($data{$chr}{$chr}{$pos1}{$pos2})){
          if ($first_line){
            print "\n$chr\t$pos1\t$chr\t$pos2\t$data{$chr}{$chr}{$pos1}{$pos2}\t$genotype";
          }else{
            $first_line++;
            print "$chr\t$pos1\t$chr\t$pos2\t$data{$chr}{$chr}{$pos1}{$pos2}\t$genotype";
	  }
        }
      }
    }
  }
}else{
  foreach my $chr1 (sort keys %chr_list){
    foreach my $chr2 (sort keys %chr_list){
      foreach my $pos1 (sort values %{$pos{$chr1}}){
        foreach my $pos2 (sort values %{$pos{$chr2}}){
          if (exists($data{$chr1}{$chr2}{$pos1}{$pos2})){
            if ($first_line){
              print "\n$chr1\t$pos1\t$chr2\t$pos2\t$data{$chr1}{$chr2}{$pos1}{$pos2}\t$genotype";
            }else{
              $first_line++;
              print "$chr1\t$pos1\t$chr2\t$pos2\t$data{$chr1}{$chr2}{$pos1}{$pos2}\t$genotype";
            }
          }
        }
      }
    }
  }
}
