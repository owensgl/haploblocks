#!/usr/bin/env perl
# -*- cperl-indent-level: 2; perl-indent-level: 2; -*-

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

# VERSION 1

=head1 SYNOPSIS

  xrqpos2ha412pos_bwa.pl [--help] INPUTVCFGZ OUTPUTPREFIX

=head1 DESCRIPTION


  This script takes the input vcf file, uses samtools to pull out the
  XRQ regions (for each SNP in the input, a small window
  around the position in the XRQ reference is selected). These regions are blasted on the HA412 genome (with bwa),
  and pulls out hits with enough mapping quality.

  A side effect of this mapping process is that sites which map to the
  same location on Ha412 will be eliminated (due to resulting low
  quality mappings).

  Arguments:

     INPUTVCFGZ   path to a .vcf.gz file containing SNPs

     OUTPUTPREFIX basename used as a prefix for temporary and final output
                  files. output files are written in the current directory.

=head1 OUTPUT

    OUTPUTPREFIX.remappedHa412.vcf.gz (+ .tbi)

       A VCF file with all input positions remapped to Ha412 positions
       identified via Blast.

    OUTPUTPREFIX.Ha412conversionstats.txt

       Consult this file to get statistics on quality of remapping.
=cut

my $help = 0;
GetOptions(
  'help'    => \$help,
);
if ($help) { pod2usage(-verbose => 2, -exit => 1); }

my $input_file = shift(@ARGV);
my $tmp_prefix = shift(@ARGV);
if (!$input_file || !$tmp_prefix) {
  pod2usage({ -message => "Missing arguments", -verbose => 1, -exit => 2});
}

sub nofail {
  my ($cmd) = @_;
  print STDERR "CMD: $cmd\n";
  system($cmd);
  if ($? == -1) {
    print STDERR "CMD:   failed to execute: $!\n";
  }
  elsif ($? & 127) {
    printf STDERR "CMD:   child died with signal %d, %s coredump\n",
	($? & 127),  ($? & 128) ? 'with' : 'without';
  }
  else {
    printf STDERR "CMD:   child exited with value %d\n", $? >> 8;
  }
  return 0;
}
  
#SLURM COMMANDS to load before running this script:
#system("module load bcftools");
#system("module load bwa");
#system("module load bedtools");

my $ha412_ref = $ENV{'HA412_REF'} // "/scratch/gowens/ref/Ha412HOv2.0-20181130.fasta";
my $xrq_ref   = $ENV{'XRQ_REF'}   // "/scratch/gowens/ref/HanXRQr1.0-20151230.fa";
#my $ha412_ref = "/home/owens/ref/Ha412HOv2.0-20181130.fasta";
#my $xrq_ref = "/home/owens/ref/HanXRQr1.0-20151230.fa";

print STDERR "Using reference: HA412 $ha412_ref\n";
print STDERR "Using reference: XRQ   $xrq_ref\n";

my $bwa = "bwa";
my $bcftools = "bcftools";
my $bases_surrounding=100; #Number of bases before and after the target site for blasting.
my $counter = 0;
my $ncores_bwa = 6; #Cores for BWA
my $min_mq = 40;


print STDERR "Reading contig listfrom fasta...";
my @contig_list;
open FASTA, $ha412_ref;
while(<FASTA>){
  chomp;
  if ($_ =~ m/\>/){
    my $contig = $_;
    $contig =~ s/\>//g;
    $contig =~ s/ /, /g;
    $contig =~ s/len/length/g;
    push(@contig_list,$contig);
  }
}
print STDERR "..done.\n";

#Turn SNPs into fastq of surrounding region
open(INPUT, "gunzip -c $input_file |");
open my $bed, '>',  "$tmp_prefix.bed";
print STDERR "Creating BED file of sites ...\n";
while(<INPUT>){
  chomp;
  if ($_ =~ m/^#/){next;}
  my @a = split(/\t/,$_);
  $counter++;
  my $chr = $a[0];
  if ($chr eq "chr"){next;}
  my $pos = $a[1];
  my $start = $pos - $bases_surrounding+1;  #+1 for bed
  my $end = $pos + $bases_surrounding+2; #+2 for bed
  if ($start < 0){next;}
  print $bed "$chr\t$start\t$end\n";
}
close INPUT;
close $bed;

if ($counter == 0) {
  die("Invalid input file given.");
}

open my $fastq, '>', "$tmp_prefix.fastq";
print STDERR "Using bedtools to create fasta of surrounding sequence..\n";
nofail("bedtools getfasta -fi $xrq_ref -bed $tmp_prefix.bed -fo $tmp_prefix.fasta");
print STDERR "Converting fasta to fastq for alignment...\n";
open FASTA, "$tmp_prefix.fasta";
my $even_odd_line = 1;
my $previous_chr;
my $previous_start;
my $previous_end;
my $previous_pos;
my $first_line;
while(<FASTA>){
  chomp;
  if ($even_odd_line % 2 == 1){ #It's the header line
    my $header = $_;
    $header =~ s/>//g;
    my @info = split(/:/,$header);
    my @start_end = split(/-/,$info[1]);
    $previous_chr = $info[0];
    $previous_start = $start_end[0];
    $previous_end = $start_end[1];
    $previous_pos = $previous_start + $bases_surrounding - 1;
  }else{
    my $seq = $_;
    my $fastq_entry;
    if($first_line){
      $fastq_entry .= "\n";
    }else{
      $first_line++;
    }
    $fastq_entry .= "\@$previous_chr:$previous_pos:$previous_start:$previous_end\n";
    $fastq_entry .= "$seq\n";
    $fastq_entry .= "+\n";
    my $qual_length = ((2*$bases_surrounding)+1);
    foreach my $i (1..$qual_length){
      $fastq_entry .= "H";
    }
    print $fastq "$fastq_entry";
    if ($counter % 10000 == 0){print STDERR "Processed $previous_chr\t$previous_pos\n";}
  }
  $even_odd_line++;
}

close $fastq;
close FASTA;

print STDERR "Aligning SNP positions using BWA...\n";
#Align to HA412 and pull out SNP position from mapping.
nofail("$bwa mem -t $ncores_bwa $ha412_ref $tmp_prefix.fastq > $tmp_prefix.sam");
open SAM, "$tmp_prefix.sam";
$counter=0;
print STDERR "Pulling SNP positions from SAM file...\n";
my %loc_hash;
while(<SAM>){
  chomp;
  if ($_ =~ m/^@/){next;}
  $counter++;
  my @a = split(/\t/,$_);
  my @info = split(/:/,$a[0]);
  my $orig_chr = $info[0];
  my $orig_snp = $info[1];
  my $orig_start = $info[2];
  my $orig_end = $info[3];
  my $flag = $a[1];
  my $target_chr = $a[2];
  my $target_pos = $a[3];
  my $mq = $a[4];
  my $code = $a[5];
  if ($flag >= 2048){next;} #Skip secondary alignments
  if ($mq < $min_mq){next;} #Skip alignments with too poor mapping quality.
  #Try to figure out if there are soft clipping at the start for determining exact alignment location
  my @softclip = split(/S/, $code);
  my $clip_value = 0;
  if ($softclip[0] =~  /^[0-9]+$/ ){
    $clip_value = $softclip[0];
  }
  my $target_snp = $target_pos + $bases_surrounding - $clip_value;
  $loc_hash{$orig_chr.$orig_snp}{"chr"} = $target_chr;
  $loc_hash{$orig_chr.$orig_snp}{"pos"} = $target_snp;
  if ($counter % 10000 == 0){print STDERR "Processed sam $orig_chr\t$orig_snp\n"}
}

close SAM;
print STDERR "Outputting remapped VCF file...\n";
open(INPUT, "gunzip -c $input_file |");
open my $unsorted_vcf, '>', "$tmp_prefix.unsorted.vcf";
my @unaligned_array;
my %alignment_stats;
my %used_sites;
$counter=0;
my $first_contig_line;
my $printed_contig_lines;
while(<INPUT>){
  chomp;
  $counter++;
  if ($. == 1){print $unsorted_vcf "$_";next;}
  if ($_ =~ m/^##contig/){
    unless($first_contig_line){
      print $unsorted_vcf "\n##INFO=<ID=XRQ,Number=1,Type=String,Description='XRQ genome location'";
      foreach my $contig (@contig_list){
        print $unsorted_vcf "\n##contig=<ID=$contig>";
      }
      print $unsorted_vcf "\n##source=Remapped to $ha412_ref using $bwa";
      $printed_contig_lines=1;
      $first_contig_line++;
    }
    next;
  }
  if ($_ =~ m/^#CHROM/){
    unless ($printed_contig_lines){
      print $unsorted_vcf "\n##INFO=<ID=XRQ,Number=1,Type=String,Description='XRQ genome location'";
      foreach my $contig (@contig_list){
        print $unsorted_vcf "\n##contig=<ID=$contig>";
      }
      print $unsorted_vcf "\n##source=Remapped to $ha412_ref using $bwa";
      $printed_contig_lines=1;
    }
  }
  if ($_ =~ m/^#/){
    print $unsorted_vcf "\n$_";
    next;
  }
  
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  if ($loc_hash{$chr.$pos}{"chr"}){
    until (!$used_sites{$loc_hash{$chr.$pos}{'chr'}.$loc_hash{$chr.$pos}{'pos'}}){ #If the specific location was already used, shift the position by one.
       $loc_hash{$chr.$pos}{'pos'}++;
     }
     print $unsorted_vcf "\n$loc_hash{$chr.$pos}{'chr'}\t$loc_hash{$chr.$pos}{'pos'}";
     foreach my $i (2..6){
       print $unsorted_vcf "\t$a[$i]";
     }
     if ($a[7] ne '.'){
       print $unsorted_vcf "\t$a[7];";
     }else{
       print $unsorted_vcf "\t";
     } #Don't print info if its only .
     print $unsorted_vcf "XRQ=$chr.$pos";
     foreach my $i (8..$#a){
       print $unsorted_vcf "\t$a[$i]";
     }
     $alignment_stats{'aligned'}++;
     $used_sites{$loc_hash{$chr.$pos}{'chr'}.$loc_hash{$chr.$pos}{'pos'}}++;
  }else{
     $alignment_stats{'unaligned'}++;
     push(@unaligned_array,"$chr\t$pos");
  }
  if ($counter % 10000 == 0){print STDERR "Processed remapped vcf $chr\t$pos\n"}
}
close INPUT;
close $unsorted_vcf;

print STDERR "Outputting remapping stats...\n";
open my $alignment_stats, '>', "$tmp_prefix.Ha412conversionstats.txt";
unless($alignment_stats{'unaligned'}){
  $alignment_stats{'unaligned'} = 0;
}
unless($alignment_stats{'aligned'}){
  $alignment_stats{'aligned'} = 0;
}
my $percent_aligned = 100 * ($alignment_stats{'aligned'} / ($alignment_stats{'aligned'} + $alignment_stats{'unaligned'}));

print $alignment_stats "Percent_snps_remapped=$percent_aligned\n";
print $alignment_stats "SNPS_remapped=$alignment_stats{'aligned'}\n";
print $alignment_stats "SNPS_not_remapped=$alignment_stats{'unaligned'}\n";
print $alignment_stats "Unmapped sites";
foreach my $site (@unaligned_array){
  print $alignment_stats "\n$site";
}
my $sorted_vcf = "$tmp_prefix.remappedHa412.vcf.gz";
#Sort using BCFtools
nofail("mkdir ./tmp_$tmp_prefix");
nofail("$bcftools sort -T ./tmp_$tmp_prefix $tmp_prefix.unsorted.vcf -O z > $sorted_vcf");
nofail("tabix -p vcf $sorted_vcf");
#Clean up temporary files
#system("rm $tmp_prefix.unsorted.vcf");
#system("rm $tmp_prefix.fasta");
#system("rm $tmp_prefix.bed");
#system("rm $tmp_prefix.fastq");
#system("rm $tmp_prefix.sam"); 

