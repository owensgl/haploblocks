# Phylogeny and dating for haploblocks

* run_inv_phylogeny*.sh
   
   Generate gene list for phylogeny/dating
* filter_phylo_genes.pl
   
   Filter gene FASTAs based on averge missing rate across samples; used in run_inv_phylogeny*.sh
* beast_run*.sh
   
   Run coalescence with BEASTv1
* fastas2BEAST1xml.pl
  
  Convert FASTA to a BEASTv1 XML; used in beast_run*.sh
* arg06.01.phylo_130140.sh
  
  Run phylogeny using SNPs in haploblock arg06.01 
* arg06.01.phylo_130140.R
  
  Plot phylogeny of arg06.01 using ggtree
* ann05.01.beast_example.r
  
  Plot BEAST tree of ann05.01 using ggtree
