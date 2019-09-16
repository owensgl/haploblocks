#!/bin/perl
#This script takes a vcf file and outputs a fasta. 
use strict;
use warnings;
use POSIX;
my $IUPAC = "TRUE";
my $subsample = 1; #Proportion of sites to subsample down.
my $raxml_filter = "TRUE"; #This filters sites where the alternate allele is only represented by heterozysgotes, making the site effectively invariant.
#This script will turn a vcf into a fasta.
#Set this to FALSE if you want to keep invariant sites.
my $remove_missing_data = "FALSE";
my $outprefix = $ARGV[0];
my(%table) = (
        'AC' => 'M',
        'CA' => 'M',
        'AG' => 'R',
        'GA' => 'R',
        'AT' => 'W',
        'TA' => 'W',
        'CG' => 'S',
        'GC' => 'S',
        'CT' => 'Y',
        'TC' => 'Y',
        'TG' => 'K',
        'GT' => 'K',
        'AA' => 'A',
        'CC' => 'C',
        'GG' => 'G',
        'TT' => 'T',
        'NN' => 'N',
	'A*' => 'A',
	'*A' => 'A',
	'C*' => 'C',
	'*C' => 'C',
	'G*' => 'G',
	'*G' => 'G',
	'T*' => 'T',
	'*T' => 'T',
	'**' => '-'
);

my %ind;
my %sites;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^##/){next;}
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#/){
    foreach my $i (9..$#a){
      $ind{$i} = $a[$i];
    }
  }else{
    my $rand = rand(1);
    if ($rand > $subsample){next;}
    my $chr = $a[0];
    my $pos = $a[1];
    my $ref = $a[3];
#    if (length($ref) > 1){next;}
    my $alt = $a[4];
    my @alts = split(/,/,$alt);
    foreach my $alt (@alts){
	    if (length($alt) > 1){goto NEXTLINE;}
    }
    my %homocounts;
    my $missingcounts;
    my $realcounts;
    foreach my $i(9..$#a){
      my @tmp = split(/:/,$a[$i]);
      my $call = $tmp[0];
      my $current_call;
      if ($call eq "."){
        $current_call .= "NN";
	$missingcounts++;
	goto PICKALLELE;
      }
      my @bases = split(/\//,$call);
      if ($bases[0] eq $bases[1]){
	if ($bases[0] eq "0"){
          $homocounts{"ref"}++;
	}elsif (($bases[0] eq "1") and ($alts[0] ne '*')){$homocounts{0}++;
	}elsif (($bases[0] eq "2") and ($alts[1] ne '*')){$homocounts{1}++;
	}elsif (($bases[0] eq "3") and ($alts[2] ne '*')){$homocounts{2}++;
	}elsif (($bases[0] eq "4") and ($alts[3] ne '*')){$homocounts{3}++;}	
      }foreach my $j (0..1){
        if ($bases[$j] eq "0"){
          $current_call .= $ref;
        }elsif($bases[$j] eq "1"){
          $current_call .= $alts[0];
        }elsif($bases[$j] eq "2"){
          $current_call .= $alts[1];
        }elsif($bases[$j] eq "3"){
          $current_call .= $alts[2]; 
        }elsif($bases[$j] eq "4"){
          $current_call .= $alts[3];		
        }elsif($bases[$j] eq '.'){
          $current_call .= "N";
	  $missingcounts++;
	}
       
      }
      PICKALLELE:
      if ($IUPAC eq "TRUE"){
	my $code;
	if ($current_call =~ m/NON_REF/){
          $code = "N";
        }else{
          $code = $table{$current_call};
        }
	unless($code){print STDERR "$current_call\n";}
        $sites{$ind{$i}} .= $code;
	if ($code ne '-'){$realcounts++;}
	
      }else{
        my @calls = split(//,$current_call);
        my $rand = int(rand(2));
	if ($calls[$rand] =~ m/NON_REF/){$calls[$rand] = "N";}
	$sites{$ind{$i}} .= $calls[$rand];
	if ($calls[$rand] ne '-'){$realcounts++;}
      }
    }
    unless($realcounts){
      chop(%sites);next;
    }
    if ($raxml_filter eq "TRUE"){
      if ((scalar keys %homocounts) < 2){
	chop(%sites);
	next;
      }
    }
    if ($remove_missing_data eq "TRUE"){
      if ($missingcounts){
	chop(%sites);
	next;
      }
    }
  }
  NEXTLINE:
}
&print_fasta;

sub print_fasta {
  my @a = sort values %ind;
  unless($sites{$a[0]}){exit;}
  foreach my $i (0..$#a){
    print ">$a[$i]\n";
    print "$sites{$a[$i]}\n";
  }
}

