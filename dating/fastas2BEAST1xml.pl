#!/bin/perl
# Take a FASTA file and convert it to a BEASTv1 xml
# run as: cat fasta.inv.list | perl fastas2BEAST1xml.sub.pl $chr $mds inv [sample_list] > $chr.$mds.inv.xml
## Kaichi Huang 2019 Apr

use warnings;
use strict;

my $chr = $ARGV[0];
my $mds = $ARGV[1];
my $inv = $ARGV[2];
my $sample_list = $ARGV[3]; # optional sample_list for subsetting
my $chainLength = 1000000; # MCMC chain length
my $screenLog = 1000; # MCMC log to screen
my $fileLog = 100; # MCMC log to file (and tree)

my @taxa;
my %seqs;
my @sample_list; # optional sample_list for subsetting

my $i = 1;
#Feed it a list of fasta files to concatenate
while (<STDIN>) {
	chomp;
	my $file = $_;
	open FILE, $file;
	my $sample;
	$seqs{$i} = {};
	while(<FILE>){
		chomp;
		next if /^\s*$/;
		if ($_ =~ m/>/) {
			$_ =~ s/>//g;
			$sample = $_;
			if ($i == 1) {
				push @taxa, $sample;
			}
		} else {
			$seqs{$i}{$sample} = $_;
		}
	}
	close FILE;
	$i = $i + 1;
}

# subsetting
if ($sample_list) {
	open LIST, $sample_list;
	while(<LIST>){
		chomp;
		push @sample_list, $_;
	}
	close LIST;
	my %sample_list = map{$_ => 1} @sample_list;
	@taxa = grep {$sample_list{$_}} @taxa;
}

my $n_gene = keys %seqs;

print "<beast version=\"1.10.4\">\n";

print "<taxa id=\"taxa\">\n";
foreach my $sample (@taxa) {
	print "<taxon id=\"$sample\"/>\n";
}
print "</taxa>\n";

foreach my $j (sort {$a<=>$b} keys %seqs) {
	print "<alignment id=\"alignment$j\" dataType=\"nucleotide\">\n";
	foreach my $sample (sort keys %{$seqs{$j}}) {
		if (grep {$_ eq $sample} @taxa) {
			print "<sequence>\n";
			print "<taxon idref=\"$sample\"/>\n";
			print "$seqs{$j}{$sample}\n";
			print "</sequence>\n";
		}
	}
	print "</alignment>\n";
}

for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<patterns id=\"gene$a.patterns\" from=\"1\" strip=\"false\">
<alignment idref=\"alignment$a\"/>
</patterns>\n";
}

print "<constantSize id=\"constant\" units=\"years\">
<populationSize>
<parameter id=\"constant.popSize\" value=\"0.04\" lower=\"0.0\"/>
</populationSize>
</constantSize>

<coalescentSimulator id=\"startingTree\">
<taxa idref=\"taxa\"/>
<constantSize idref=\"constant\"/>
</coalescentSimulator>

<treeModel id=\"treeModel\">
<coalescentTree idref=\"startingTree\"/>
<rootHeight>
<parameter id=\"treeModel.rootHeight\"/>
</rootHeight>
<nodeHeights internalNodes=\"true\">
<parameter id=\"treeModel.internalNodeHeights\"/>
</nodeHeights>
<nodeHeights internalNodes=\"true\" rootNode=\"true\">
<parameter id=\"treeModel.allInternalNodeHeights\"/>
</nodeHeights>
</treeModel>

<treeLengthStatistic id=\"treeLength\">
<treeModel idref=\"treeModel\"/>
</treeLengthStatistic>

<coalescentLikelihood id=\"coalescent\">
<model>
<constantSize idref=\"constant\"/>
</model>
<populationTree>
<treeModel idref=\"treeModel\"/>
</populationTree>
</coalescentLikelihood>\n";

for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<strictClockBranchRates id=\"gene$a.branchRates\">
<rate>
<parameter id=\"gene$a.clock.rate\" value=\"1.0\"/>
</rate>
</strictClockBranchRates>
<rateStatistic id=\"gene$a.meanRate\" name=\"gene$a.meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">
<treeModel idref=\"treeModel\"/>
<strictClockBranchRates idref=\"gene$a.branchRates\"/>
</rateStatistic>\n";
}

for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<HKYModel id=\"gene$a.hky\">
<frequencies>
<frequencyModel dataType=\"nucleotide\">
<alignment idref=\"alignment$a\"/>
<frequencies>
<parameter id=\"gene$a.frequencies\" dimension=\"4\"/>
</frequencies>
</frequencyModel>
</frequencies>
<kappa>
<parameter id=\"gene$a.kappa\" value=\"2.0\" lower=\"0.0\"/>
</kappa>
</HKYModel>
<siteModel id=\"gene$a.siteModel\">
<substitutionModel>
<HKYModel idref=\"gene$a.hky\"/>
</substitutionModel>
<gammaShape gammaCategories=\"4\">
<parameter id=\"gene$a.alpha\" value=\"0.5\" lower=\"0.0\"/>
</gammaShape>
</siteModel>
<statistic id=\"gene$a.mu\" name=\"mu\">
<siteModel idref=\"gene$a.siteModel\"/>
</statistic>\n";
}

for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<treeDataLikelihood id=\"gene$a.treeLikelihood\" useAmbiguities=\"false\">
<partition>
<patterns idref=\"gene$a.patterns\"/>
<siteModel idref=\"gene$a.siteModel\"/>
</partition>
<treeModel idref=\"treeModel\"/>
<strictClockBranchRates idref=\"gene$a.branchRates\"/>
</treeDataLikelihood>\n";
}

print "<operators id=\"operators\" optimizationSchedule=\"default\">\n";
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<scaleOperator scaleFactor=\"0.75\" weight=\"1\">
<parameter idref=\"gene$a.kappa\"/>
</scaleOperator>
<scaleOperator scaleFactor=\"0.75\" weight=\"1\">
<parameter idref=\"gene$a.alpha\"/>
</scaleOperator>\n";
}
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<scaleOperator scaleFactor=\"0.75\" scaleAll=\"true\" ignoreBounds=\"true\" weight=\"3\">
<parameter idref=\"treeModel.allInternalNodeHeights\"/>
</scaleOperator>\n";
}
print "<subtreeSlide size=\"1.0\" gaussian=\"true\" weight=\"30\">
<treeModel idref=\"treeModel\"/>
</subtreeSlide>
<narrowExchange weight=\"30\">
<treeModel idref=\"treeModel\"/>
</narrowExchange>
<wideExchange weight=\"3\">
<treeModel idref=\"treeModel\"/>
</wideExchange>
<wilsonBalding weight=\"3\">
<treeModel idref=\"treeModel\"/>
</wilsonBalding>
<scaleOperator scaleFactor=\"0.75\" weight=\"3\">
<parameter idref=\"treeModel.rootHeight\"/>
</scaleOperator>
<uniformOperator weight=\"30\">
<parameter idref=\"treeModel.internalNodeHeights\"/>
</uniformOperator>
<scaleOperator scaleFactor=\"0.75\" weight=\"3\">
<parameter idref=\"constant.popSize\"/>
</scaleOperator>\n";
print "</operators>\n";

print "<mcmc id=\"mcmc\" chainLength=\"$chainLength\" autoOptimize=\"true\" operatorAnalysis=\"$chr.$mds.$inv.ops\">
<joint id=\"joint\">
<prior id=\"prior\">\n";
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<logNormalPrior mu=\"1.0\" sigma=\"1.25\" offset=\"0.0\">
<parameter idref=\"gene$a.kappa\"/>
</logNormalPrior>
<exponentialPrior mean=\"0.5\" offset=\"0.0\">
<parameter idref=\"gene$a.alpha\"/>
</exponentialPrior>\n";
}
print "<gammaPrior shape=\"10.0\" scale=\"0.004\" offset=\"0.0\">
<parameter idref=\"constant.popSize\"/>
</gammaPrior>
<coalescentLikelihood idref=\"coalescent\"/>\n";
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<strictClockBranchRates idref=\"gene$a.branchRates\"/>\n";
}
print "</prior>
<likelihood id=\"likelihood\">\n";
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<treeDataLikelihood idref=\"gene$a.treeLikelihood\"/>\n";
}
print "</likelihood>
</joint>
<operators idref=\"operators\"/>\n";
print "<log id=\"screenLog\" logEvery=\"$screenLog\">
<column label=\"Joint\" dp=\"4\" width=\"12\">
<joint idref=\"joint\"/>
</column>
<column label=\"Prior\" dp=\"4\" width=\"12\">
<prior idref=\"prior\"/>
</column>
<column label=\"Likelihood\" dp=\"4\" width=\"12\">
<likelihood idref=\"likelihood\"/>
</column>
<column label=\"rootHeight\" sf=\"6\" width=\"12\">
<parameter idref=\"treeModel.rootHeight\"/>
</column>
</log>\n";
print "<log id=\"fileLog\" logEvery=\"$fileLog\" fileName=\"$chr.$mds.$inv.log\" overwrite=\"false\">
<joint idref=\"joint\"/>
<prior idref=\"prior\"/>
<likelihood idref=\"likelihood\"/>
<parameter idref=\"treeModel.rootHeight\"/>
<treeLengthStatistic idref=\"treeLength\"/>
<parameter idref=\"constant.popSize\"/>\n";
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<parameter idref=\"gene$a.kappa\"/>
<parameter idref=\"gene$a.alpha\"/>\n";
}
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<parameter idref=\"gene$a.clock.rate\"/>\n";
}
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<rateStatistic idref=\"gene$a.meanRate\"/>\n";
}
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<treeDataLikelihood idref=\"gene$a.treeLikelihood\"/>\n";
}
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<strictClockBranchRates idref=\"gene$a.branchRates\"/>\n";
}
print "<coalescentLikelihood idref=\"coalescent\"/>
</log>\n";
print "<logTree id=\"treeFileLog\" logEvery=\"$fileLog\" nexusFormat=\"true\" fileName=\"$chr.$mds.$inv.trees\" sortTranslationTable=\"true\">
<treeModel idref=\"treeModel\"/>\n";
for (my $a=1; $a<=$n_gene; $a=$a+1){
	print "<trait name=\"rate\" tag=\"gene$a.rate\">
<strictClockBranchRates idref=\"gene$a.branchRates\"/>
</trait>\n";
}
print "<joint idref=\"joint\"/>
</logTree>
</mcmc>\n";

print "<report>
<property name=\"timer\">
<mcmc idref=\"mcmc\"/>
</property>
</report>\n";

print "</beast>\n";
