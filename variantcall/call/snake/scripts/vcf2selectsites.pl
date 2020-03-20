#!/usr/bin/env perl
# -*- cperl-indent-level: 2; perl-indent-level: 2; -*-

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


# VERSION 2

=head1 SYNOPSIS

  vcf2selectsites.pl [--help] [--header|--noheader] SNPOUT INDELOUT SITEFILE [SITEFILE]*

=head1 DESCRIPTION

  Select VCF sites from stdin matching entries provided in one or more
  sitefiles.  The site list is specified as a positional argument and
  contains at least the columns "chr pos type" identifying the sites
  that should be kept. the records with type "snp" (case insensitive)
  will be sent to SNPOUT, and those with "indel" sent to INDELOUT.

  By default it is assumed that site files will each have a one line
  header on the first line.

=cut

use constant COL_CHR => 0;
use constant COL_POS => 1;
use constant COL_TYPE => 2;


my $has_header = 1;
my $help = 0;
GetOptions(
  'header!' => \$has_header,
  'help'    => \$help,
);
if ($help) { pod2usage(-verbose => 2, -exit => 1); }

my $snpout = shift(@ARGV);
my $indelout = shift(@ARGV);

if (!$snpout || !$indelout) {
  pod2usage({ -message => "Missing SNPOUT or INDELOUT", -verbose => 1, -exit => 2});
}

#site list file should have a header.
my $site_file = shift(@ARGV);

if (!$site_file) {
  pod2usage({ -message => "At least one site file must be provided.", -verbose => 1, -exit => 2});
}

open(SNPOUT, ">", $snpout) or die("Can't open file $snpout");
open(INDELOUT, ">", $indelout) or die ("Can't open file $indelout");

my %snp_sites;
my %indel_sites;

my %col_map = (
  "chr"        => COL_CHR,
  "chrom"      => COL_CHR,
  "chromosome" => COL_CHR,
  "pos"        => COL_POS,
  "position"   => COL_POS,
  "type"       => COL_TYPE,
  "typ"        => COL_TYPE,
);


while (defined $site_file) {
    open(LIST, $site_file) or die("Can't open file $site_file");
    while(<LIST>){
	chomp;
	if ($. == 1){
	  my @a = split(/\t/,$_);
	  if (0+@a < 3) { die("Expected 3 columns"); }
	  if ($has_header) {
	    foreach my $hdri (0 .. $#a) {
	      my $hdr_chk = lc($a[$hdri]);
	      if (exists($col_map{$hdr_chk})) {
		if ($hdri != $col_map{$hdr_chk}) {
		  die("Expected column $hdr_chk in position $col_map{$hdr_chk}. Found in position $hdri");
		}
	      }
	    }
	    next;
	  }
	}

	if ($_ =~ m/#/){next;}
	my @a = split(/\t/,$_);
	if (lc($a[COL_TYPE]) eq "snp") {
	  $snp_sites{$a[COL_CHR]."_".$a[COL_POS]}++;
	} elsif (lc($a[COL_TYPE]) eq "indel") {
	  $indel_sites{$a[COL_CHR]."_".$a[COL_POS]}++;
	} else {
	  die("Unknown variant type '$a[COL_TYPE]' in $site_file line $.: $_");
	}
    }
    close(LIST);
    $site_file = shift(@ARGV);
}

while(<STDIN>){
  chomp;
  if ($. == 1){
    print SNPOUT "$_\n";
    print INDELOUT "$_\n";
    next;
  }
  if ($_ =~ m/^##/) {
    print SNPOUT "$_\n";
    print INDELOUT "$_\n";
  } elsif ($_ =~ m/^#/) {
    # header
    my @a = split(/\t/,$_);
    if (lc($a[COL_CHR]) !~ /^#?chrom$/ || lc($a[COL_POS]) ne "pos") {
      die("Expected CHROM, POS in VCF. Found ".$a[COL_CHR].", ".$a[COL_POS]);
    }
    print SNPOUT "$_\n";
    print INDELOUT "$_\n";
  } else {
    my @a = split(/\t/,$_);
    my $chr = $a[COL_CHR];
    my $pos = $a[COL_POS];
    if ($snp_sites{"${chr}_$pos"}){
      print SNPOUT "$_\n";
    }
    if ($indel_sites{"${chr}_$pos"}) {
      print INDELOUT "$_\n";
    }
  }
}
