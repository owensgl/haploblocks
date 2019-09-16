#!/bin/perl
use warnings;
use strict;

#This script is for creating tree files that are tested for each gene. 
my $species = $ARGV[0]; #String for species name.
my $group_1 = $ARGV[1]; #File name of inversion group 1 list
my $group_2 = $ARGV[2]; #File name of inversion group 2 list
my $other_samples = $ARGV[3]; #File name of other samples;
my $sampleinfo = $ARGV[4]; #Sample info file

my $g1;
my $g2;
#Load in inversion samples
open G1, $group_1;
while(<G1>){
	chomp;
	if ($. == 1){
		$g1 = $_;
	}else{
		$g1 .= ",$_";
	}
}
close G1;
open G2, $group_2;
while(<G2>){
	chomp;
	if ($. == 1){
		$g2 = $_;
	}else{
		$g2 .= ",$_";
	}
}
close G2;
#Load in species identity
my %info;
open INFO, $sampleinfo;
while(<INFO>){
	chomp;
	if ($. == 1){next;}
	my @a = split(/\t/,$_);
	my $sample = $a[2];
	my $species = $a[3];
	#Put all perennials in a single group;
	$species =~ s/Gig/Per/;
	$species =~ s/Div/Per/;
	$species =~ s/Dec/Per/;
	$species =~ s/Gro/Per/;
	$info{$sample} = $species;
}
close INFO;

my %species_groups;
#Load in other samples
open OTHER, $other_samples;
while(<OTHER>){
	chomp;
	my $species = $info{$_};
	if ($species_groups{$species}){
		$species_groups{$species} .= ",$_";
	}else{
		$species_groups{$species} = $_;
	}
}
#Print out the actual phylogenies we want to test.
if ($species eq "Ann"){
	print "((((($g1),($g2)),$species_groups{'Arg'}),($species_groups{'Pet'})),($species_groups{'Per'}));\n";
	print "((((($g2),($species_groups{'Arg'})),$g1),($species_groups{'Pet'})),($species_groups{'Per'}));\n";
	print "((((($g2),($species_groups{'Arg'})),$species_groups{'Pet'}),($g1)),($species_groups{'Per'}));\n";
	print "((((($g1),($species_groups{'Arg'})),$g2),($species_groups{'Pet'})),($species_groups{'Per'}));\n";
	print "((((($g1),($species_groups{'Arg'})),$species_groups{'Pet'}),($g2)),($species_groups{'Per'}));";
}elsif ($species eq "Arg"){
	print "((((($g1),($g2)),$species_groups{'Ann'}),($species_groups{'Pet'})),($species_groups{'Per'}));\n";
        print "((((($g2),($species_groups{'Ann'})),$g1),($species_groups{'Pet'})),($species_groups{'Per'}));\n";
        print "((((($g2),($species_groups{'Ann'})),$species_groups{'Pet'}),($g1)),($species_groups{'Per'}));\n";
        print "((((($g1),($species_groups{'Ann'})),$g2),($species_groups{'Pet'})),($species_groups{'Per'}));\n";
        print "((((($g1),($species_groups{'Ann'})),$species_groups{'Pet'}),($g2)),($species_groups{'Per'}));";
}
