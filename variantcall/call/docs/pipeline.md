[Back to start](./README.md)

Variant-Calling Pipeline
========================

   The pipeline operations is programmed into `snake/Snakefile`. The program defines
   a series of `targets` (or commands) that can be produced, which we outline here. Each
   phase of the pipeline has its own set of targets. The order of the targets reflect the
   steps of the pipeline:

   | Target | Description |
   | ------ | ----------- |
   | *all_bam*    | Aligning each sample sequence against the reference (ngm alignment, removal of duplicates, markers trimmed). This will download the sequence files if necessary. |
   | *all_merged* | Merge bams from the same sample into a single bam. The notion of "same sample" is determined based on read-group information. |
   | *all_bam_stats* | Produce statistics on produced merged bam files. i.e. `samtools stats` |
   | *all_bases_mapped* | Produce a listing of sample names in order of # of bases mapped, descending |
   | *all_gvcf* | Produce per-sample variant calls (with GATK HaplotypeCaller). |
   | *all_gendb* | Combine variant calls into efficient database format (GenomicsDB) |
   | *all_vcf_chrom* | Perform jointcalling on samples from gendb database (produces one unfiltered (i.e. RAW) VCF) per chromosome |
   | *all_vcf_chrom_sites_only* | Produces a smaller version of the unfiltered VCF where genotypic information has been removed. |
   | *all_vcf_goldset* | Produce a hard-filtered version of the raw VCF subset (which we call a "gold set"). The parameters for filtering are configurable. |
   | *all_variant_models* | Produce a recalibration model computed over one VCF with RAW SNPs, and using another hard-filtered VCF goldset |
   | *all_vcf_chrom_filtered* | Apply the recalibrated model over the set of raw SNPs. Produces a the final VCF for this pipeline |

   The pipeline command template is the following:

   ```bash
   <pipeline_cmd> [<pipeline_opt>...] <target1> [<target2> [<target3> [... <targetN>]]]
   ```

   Which will produce all the demanded targets in parallel. `<target>` would be either:

   - one of the command names presented in the table above
   - a path to a file that needs to be created

   For `<pipeline_cmd>`, there are three choices offered, which will depend on your environment.

   ```bash
   cd "$REPODIR"/snake/

   # option1: produce TARGET by using slurm and sbatch (compute-canada)
   # PREFIX is a short name which will be prepended to all submitted jobs.
   ./compute-canada-go --run-prefix PREFIX TARGET

   : -- OR --

   # option2: produce TARGET on compute canada, but directly on the
   # current node (which could be a login node, or an `salloc` node).
   # the --local-only flag instructs the script to avoid sbatch calls.
   ./compute-canada-go --local-only -j NTHREADS TARGET

   : -- OR --
   # option3: compile the tools in a non-compute-canada environment. e.g. your own
   #          multicore server.
   ./local-go -j NTHREADS TARGET
   ```

   Options 2 and 3 work similarly, in the sense that the jobs run on
the same node as the controller -- no new worker nodes are
allocated. The difference between them is that some paths will change
when creating singularity containers -- on compute-canada, paths
`/scratch`, `/localscratch`, and `/project` need to be bind-mounted
into the worker nodes.

   Option 1 is the common case when running with on a cluster. It
might be beneficial, in cases where a shorter turnaround is needed, to
invoke rules directly on an already allocated node, so option 2 is
given as an example above. In option1, the argument `--run-prefix
PREFIX` will cause all submitted job names to be prefixed with
`PREFIX`. This allows jobs to be killed in a single swoop, based on their names.

**Parameters and options**

   Any parameter provided to the pipeline commands will be passed to
the underlying scheduler (snakemake), so it is in one's best interest to familiarize oneself
with its capabilities.

A useful option, for instance, `--restart-times N`, instructs the scheduler to
reschedule any job that has failed, up to N times. Each resubmission
provides an opportunity to allocate additional RAM or time to the job
(an equation which is defined inside the `Snakefile`.

Another option which you will want to abuse is `--config CFG`, which allows you to
override parameters in the config with those given on the command line. The following command will generate VCF files (one per
chromosome), comparing the cohort set defined by samples whose names are
listed in `my_cohort.txt` (one name per line, single column)

```bash
./compute-canada-go --run-prefix K1C6_ --config "samplenames='./samples/my_cohort.txt'" all_chromosome_vcf --restart-times 2
```

Walkthrough -- A typical variant call
=====================================

We will walk you through the set of commands to invoke to produce a VCF for
a variant call:

1. Define a list of sample names to include in the cohort.

   The first step is to create a plain text file
   (e.g. `my_cohort.txt`) where the names of all the samples to
   include are provided, one per line.  Each name provided should
   match one of the entries in your yaml database of samples (see
   [Importing Samples](./importing_samples.md)).
   
   Ex:
   ```text
   # my_cohort.txt
   ANN0001
   PET1001
   ARG1234
   # lines can start with `#` which turns them into comments

   # empty lines are allowed too...
   
   BLAH555
   ```

1. Compute BAM statistics

   We tell the pipeline to align the samples, merge them, and provide us with a listing of samples, ordered by number of bases
   mapped. This is indicative of the depth of the samples. We will use the samples with the most bases mapped to calculate a goldset.

   ```bash
   ./compute-canada-go --run-prefix GOLD --config "samplenames='./my_cohort.txt'" all_bases_mapped
   ```

   Note that we didn't have to run the alignment and merging steps explicitly because the target `all_bases_mapped` depends on
   them, so they are computed automatically in the process.

3. Form a subset of samples to make a gold set

   The previous step will output a file containing statistics on each sample. We pick the first 20 samples, and form a new
   cohort:

   ```bash
   #gold_cohort.txt
   ANN1000
   ...
   ANN1019
   ```  
   _Note: use your 20 highest-quality samples_

4. Make a contig list

   We want to select high quality SNPs from chromosomal DNA only, so before we build our hard-filtered VCF,
   we make a file containing the contigs of interest (windows over which the VCF will be computed):

   ```text
   # contig-list-gold.txt
   HanXRQChr01 1 1000000 153905722
   HanXRQChr01 1000001 2000000 153905722
   ...
   HanXRQChr00c1511 1 641 641
   HanXRQChr00c1512 1 518 518
   #HanXRQMT 1 301004 301004
   #HanXRQCP 1 151101 151101
   ```

   The syntax for a contig list is: `<chromosome> <start_bp> <end_bp>
   <chrom_len>`. Here, we've included the whole genome, except that
   we've commented out the chloroplasmic and mitochondrial sections.
   You wouldn't normally generate this file by hand: you can use
   `snake/scripts/generate_contigs.sh REF_GENOME.fa` to produce it.

5. Hard filter a gold set.

   Now, form your goldset VCF file:

   ```bash
   ./compute-canada-go --contiglist data/gwas/contig-list-gold.txt \
                    --run-prefix GOLDFILTER \
                    --config "samplenames='./gold_cohort.txt'" \
                    --config vcf_batch_size_kbp=50000 \
                    all_vcf_goldset
   ```

   The `vcf_batch_size_kbp` config flag controls how big of a window
   will be computed at once per job. It only affects the fan-out level
   of the job. Picking a number too low will increase the total number
   of jobs. Picking a number too high will make jobs run longer (and
   possibly run over the time limit). We have found 50000 to be a good
   amount for 20 samples.

   This will produce 2 files: goldset.snps.tgz, goldset.indels.tgz (plus 2 index files). Which
   we will use in later steps.

6. Compute a raw VCF over the entire cohort.

   We can compute a set of RAW snps over the entire cohort. This does not yet use the goldset produced earlier.
   
   ```bash
   ./compute-canada-go --run-prefix RAW1 --config "samplenames='./my_cohort.txt'" --config vcf_batch_size_kbp=12000 all_vcf_chrom
   ```

   Note that we are providing a `vcf_batch_size_kbp` that is lower
   than when we computed the goldset. This is because we have more
   samples to compute, so we have to reduce the amount of work per job
   to fit the time allocation. The number you need to provide involves
   some guesswork to provide optimal performance. Aim in the ballpark
   of 5000 for cohorts of more than 1200 samples, and 50000 for
   cohorts with less than 100 samples. The work needed is exponential in the number
   of samples in the cohort.

7. Filter SNPs

   As the final step, use the goldset computed earlier to produce a filtered (VQSR) set of
   SNPs. you will have to modify paths to account for the locations where the files have been
   generated by the previous steps.

   ```bash
   ./compute-canada-go --contiglist data/gwas/contig-list.txt \
                    --run-prefix FILTER1 \
                    --config "samplenames='./samples/my_cohort.txt'" \
                    --config "goldset_snps='data/gwas/gold/goldset.snps.vcf.gz'" \
                    --config "goldset_indels='data/gwas/gold/goldset.indels.vcf.gz'" \
                    --config vcf_batch_size_kbp=12000 \
                    all_vcf_chrom_filtered
   ```

   This will recalibrate a model using the goldset over the entire dataset, and apply a round of VQSR
   on the raw snps produced previously.

   Note here that the numeric configuration parameters provided to
   this command must match those entered earlier
   (esp. `vcf_batch_size_kbp=12000`). This is because the scheduler
   can only reuse previously computed results if the parameters match
   those of a subsequent invocation.


8. Final round of filtering for GWAS

   TODO describe how to apply a custom filter to make the dataset ready for GWAS.

Job and resource Management
===============

TODO get status, get ETA, start/pause/stop/resume. how to adapt resource requirements for each job (ram, disk). how to control parallelism and threads.

