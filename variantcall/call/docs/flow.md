
Pipeline Steps
==============


Input data:

- Reference Genome HanXRQr1.0-20151230 (2.9G on disk).

- Illumina configuration (markers)

- Transposable Elements (maps of contigs containing transposable elements)

- Sequence Data (FastQ)  
  2400 individuals. Sequencing calibrated for avg coverage of 8.  
  Some individuals are sequenced twice (or more) to get more coverage.  
  Each sequence is roughly 3-4GB of fastq (gzip'd). Some plants will require maybe twice that when another run is added.

Step 1: Align each sequencing run

- ngmAlign. takes about 30min and under 80GB ram to align 3GB of sequence data.
  32 core (Intel(R) Xeon(R) CPU E5-2683 v4 @ 2.10GHz).

    - Pre processing: trimming raw sequences for low scores
    - Pre processing: removing illumina markers
    - Post processing: convert from sam to bam
    - Post processing: concatenate partial bams
    - Post processing: sort mapped reads. sorted based on mapping position in reference
    - Post processing: mark duplicates (illumina pcr artefacts)
    - Post processing: reheader bam with sample metadata
    - Post processing: generate bam index

Step 2: Make sure produced bams meet the provided quality requirements

    - inspect coverage obtained. if more reads are needed, resequence
      plant and run additional step 1's. Sequencing is a random process, and
      some sections of the genome may be under-represented.

    - we discovered a labelling error at this stage. is FOO.BAM and
      BAR.BAM the same plant?

    - Run a BAM validation tool.

Step 3: Merge all the created bams

    - Map each bam in step 2 to a sample name to use in the cohort

    - Aligned files from step 2 are merged if they have the same name.  
      This is a N-way merge of sorted bams. (works like merge sort.)

    - Merged BAM is on average 10GB  (this is per plant).  
      min=5.75GB max=29.82GB n=587 avg=9.58GB sum=5.6TB std=2.34GB


Step 4: Generate Haplotypes per-sample (HaplotypeCaller -ERC GVCF)

    - Produce GVCF file per-sample.  
      https://www.biorxiv.org/content/biorxiv/early/2017/11/14/201178.1.full.pdf

    - Perform local realignment of "active regions" (hot spots of variation)

    - For each position, determine likelihood of a read going on one
      haplotype vs the other. Then determine the raw genotype likelihood of
      each variant.

Step 5: Combine all of step 4 using GenotypeGVCFs

    - Exponential runtime in number of samples. But massively paralellizable.

    - At each position of the input gVCFs, this tool will combine all
      spanning records, produce correct genotype likelihoods,
      re-genotype the newly merged record, and then re-annotate it.

    - Produces large VCFs of SNPs and Indels (but they need to be
      filtered.)

    - Variant calling is done after this step completes. Filtering starts.

Step 6: Decide on model for true variants

    - This step consists in building a model of what a true variant should be.

    - We re-run all of step 5 on a small subset cohort. We pick the N samples
      which have been sampled the best to make that cohort.

    - We perform a hard filter on the variants of the raw vcf and call that
      the gold-set.

Step 7: Run the variant recalibration based on model.

- We combine two different runs of the variant calling.

        1) raw vcf obtained over entire sample set
        2) hard-filtered gold set vcf used as a source of true snps.

- Run VariantRecalibrator

    From: https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.8.1/org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator.php


    > Build a recalibration model to score variant quality for filtering purposes
This tool performs the first pass in a two-stage process called Variant Quality Score Recalibration (VQSR). Specifically, it builds the model that will be used in the second step to actually filter variants. This model attempts to describe the relationship between variant annotations (such as QD, MQ and ReadPosRankSum, for example) and the probability that a variant is a true genetic variant versus a sequencing or data processing artifact. It is developed adaptively based on "true sites" provided as input, typically HapMap sites and those sites found to be polymorphic on the Omni 2.5M SNP chip array (in humans). This adaptive error model can then be applied to both known and novel variation discovered in the call set of interest to evaluate the probability that each call is real. The result is a score called the VQSLOD that gets added to the INFO field of each variant. This score is the log odds of being a true variant versus being false under the trained Gaussian mixture model.

   - This offers a series of cutoff VQSLOD values which allow you to
     flexibly choose your quantity of false positives and false negatives.
     One metric you have to determine whether you've included enough is the
     Ti/Tv ratio. But there are other ways to determine how close you got:

     https://software.broadinstitute.org/gatk/documentation/article.php?id=6308


Step 8: Assign each variant to a tranche

     Mark each variant in the callset with the tranche it falls into. The tranche is simply a range to VQSLOD values.

Step 9: Apply single pass filters over dataset repeatedly.

Step 10: GWAS

