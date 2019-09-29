#!/usr/bin/env perl
# -*- cperl-indent-level: 2; perl-indent-level: 2; -*-
use warnings;
use strict;
use Time::HiRes qw( gettimeofday );

# Greg Owens 2018 This script takes a GATK4 vcf and pulls out site
# specific stats. It only outputs the top 3 alternate allele
# frequencies, even if there are more alternate alleles. It doesn't
# work allele specific annotations.

# Contributors: Jean-Sebastien Legare 2018
#
#
# Usage:
#   This script receives VCF data on stdin, and outputs tabular information on
#   each site of the file. output has a header.
#
#   CHR POS TYPE ALTLENGTH QUAL HET_PERCENT INFO+ AF0 [AF1] [AF2]
#
#   STATS is: AN BaseQRankSum ClippingRankSum DP ExcessHet FS InbreedingCoeff MQ MQRankSum QD ReadPosRankSum SOR
#

use constant NA => "+nan";
use constant COL_CHR => 0;
use constant COL_POS => 1;
use constant COL_REF => 3;
use constant COL_ALT => 4;
use constant COL_QUAL => 5;
use constant COL_FILTER => 6;
use constant COL_INFO => 7;

my %col_map = (
  "chrom"  => COL_CHR,
  "pos"    => COL_POS,
  "ref"    => COL_REF,
  "alt"    => COL_ALT,
  "qual"   => COL_QUAL,
  "filter" => COL_FILTER,
  "info"   => COL_INFO,
);

my %s;
my $flag_filter_multi = 0; #If true, removes non-biallelic sites.
my $flag_count_het = 0; #If true, it will also output observed heterozygosity
my $start_ts = gettimeofday();

my $last_count = 0;
my $last_ts = $start_ts;
sub show_progress { # count, chr, pos
  my $now = gettimeofday();
  if (($now - $last_ts) < 5.0) {
    return;
  }
  my $avg_speed = ($_[0] - $last_count) / ($now - $last_ts + 0.000001);
  $last_count = $_[0];
  $last_ts = $now;
  print STDERR sprintf("[%10.3f] vcf2stats progress: %s\t%s\t(%.3f sites/min)\n", ($now - $start_ts), $_[1], $_[2], $avg_speed * 60.0);
}

# Defines columns to select, as well as their output order
my @info_cols = ("AN", "BaseQRankSum", "ClippingRankSum", "DP", "ExcessHet", "FS", "InbreedingCoeff", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR");
foreach (@info_cols) {
    $s{$_} = NA;
}

print "chr\tpos\ttype\talt_length\tqual\tfilter\thet_percent";
foreach my $info_col (@info_cols) {
  print "\t$info_col";
}
foreach my $i (0..2){
  print "\tAF_$i";
}
print "\n";

my $count = 0;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^##/){next;}
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#/){
    # header sanity check
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  PET0316 ...
    (my $hdr = lc($_)) =~ s/^#+//g;
    my @hdr_cols = split(/\s+/, $hdr);
    foreach my $hdri (0 .. $#hdr_cols) {
      my $hdr_chk = lc($hdr_cols[$hdri]);
      if (exists($col_map{$hdr_chk})) {
	if ($hdri != $col_map{$hdr_chk}) {
	  die("Expected column $hdr_chk in position $col_map{$hdr_chk}. Found in position $hdri");
	}
      }
    }
    next;
  } else {
    $count++;
    my $chr = $a[COL_CHR];
    my $pos = $a[COL_POS];
    if ($count % 10000 == 0) { show_progress($count, $chr, $pos); }
    my $ref = $a[COL_REF];
    my $alt = $a[COL_ALT];
    my $filter = $a[COL_FILTER];

    my $type = "SNP";
    my @alts = split(/,/,$alt);
    if ($flag_filter_multi){
      if ($alts[1]){next;}
    }
    my $alt_length = scalar(@alts);
    if (length($ref) > 1){
      $type = "Indel";
    }
    foreach my $i (0..$#alts){
      if (length($alts[$i]) > 1){
        $type = "Indel";
      }
    }
    my $qual = $a[COL_QUAL];

    # reset info fields
    foreach my $info_col (@info_cols) {
	$s{$info_col} = NA;
    }

    my %AF;
    my $info = $a[COL_INFO];
    my @infos = split(/;/,$info);
    foreach my $infos (@infos){
      my @parts = split(/=/,$infos);
      if ($s{$parts[0]}){
        $s{$parts[0]} = $parts[1];
      }
      if ($parts[0] eq "AF"){
        my @AFs = split(/,/,$parts[1]);
        foreach my $i (0..$#alts){
          $AF{$i} = $AFs[$i];
        }
      }
    }
    my $het = 0;
    my $total = 0;
    my $het_percent = NA;

    if ($flag_count_het) {
      foreach my $i (9..$#a){
	my @calls = split(/:/,$a[$i]);
	if ($calls[0] eq './.'){next;}
	if ($calls[0] eq '.'){next;}
	my @genotypes = split(/\//,$calls[0]);
	if ($genotypes[0] ne $genotypes[1]){
	  $het++;
	}
	$total++;
      }
      if ($total >= 1) {
	$het_percent = $het*100/$total;
      }
    }
    print "$chr\t$pos\t$type\t$alt_length\t$qual\t$filter\t$het_percent";

    foreach my $info_col (@info_cols){
      print "\t$s{$info_col}";
    }

    my @AF_sorted =  (sort {$b <=> $a} values %AF);
    foreach my $i (0..2){
      if ($AF_sorted[$i]){
        print "\t$AF_sorted[$i]";
      }else{
        print "\t" . NA;
      }
    }
    print "\n";
  }
}
