# Tools used in various parts.

* concatenate_fasta.pl

  Concatenate multiple genes fasta files and retain partition information for BEAST
* genotype_inv_from_pcasites_printsites.pl

  Extracts haploblock diagnostic sites from gvcf files and saves in a tidy fashion.
* genotypes2vcf.pl

  Creates a pseudo-vcf for haploblock genotypes

* gvcf2fasta_nogaps.pl

  Extracts genotypes from a gcvf to turn into a fasta file.
  
* make_fasta.sh

  Wrapper script for creating gvcf and turning it into a fasta for a list of genes.
* pull_out_genes.pl

  Extracts genes within regions. 
* remap_call_env.sh

  Example of wrapper script for remapping XRQ vcf to HA412v2 positions.
* remove_inversion_sites.sh

  Removes sites within and highly correlated to haploblocks.
* select_samples_inv.pl

  Extracts a random selection of samples that are homozygous for a particular haploblock.
* tped2allelefreq.pl

  Extracts allele frequencies from tped files.

* vcf2DerivedAncestral_type2.pl
  
  Polarizes vcf file to derived and ancestral based on outgroup.
* vcf2fasta_basic.pl

  Converts a vcf file into a basic fasta.
* xrqpos2ha412pos_bwa.pl
  
  Extracts surrounding bases for each locus and remaps to Ha412v2 using BWA. This creates a new VCF with the old location saved in the INFO column.
