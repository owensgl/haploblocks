# Genome environment associations

* env_gea_genomewide.R
   
   Printing of all GEA with and without haploblock regions removed. Selecting top candidate genes.
   
* environment_variables.txt

   All environmental variables with general category

* print_all_gea.R
    
    Printing all haploblock environment associations summary
    
* run_baypass_core_model.pbs

   Running the core model of BayPass for soil GEA

* run_baypass_cov_model.pbs

   Running the covariate model of BayPass for soil GEA
   
* vcf2vertical_sunflower_modifies.pl

   Filters for vcf file parameters: Quality, Mapping Quality, Depth of Coverage and Genotype Quality and generate genotype file
   
* snp.coverage.pl

   Filters the genotype files from vcf2vertical for minor allele frequency, calling rate and heterozygosity
   
* Baypass_input_file_prep.pl

   Converts the genotype file into baypass input file (allele count)
   
* BayPass_job_arrays_Graham_cedar.sh

   Runs baypass 5k genotype files as a job array on compute canada servers (Graham and Cedar) for climate GEA
   
* baypass_outlier_plot.R

  Summerizes BayPass outputs, calulates outlier percentage and plots the results
   
* haploblocks_input_format.R

   Converts haploblocks vcf files (012 format) to baypass input files and allele frequency
