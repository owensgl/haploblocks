[Back to start](./README.md)

VC Tool
==============

Directory: `./vc`

Performs per-sample variant calling against a reference.

You provide an input BAM file, as well as a reference, and it will produce a g.vcf file as well as an index file.

Multiple of those .g.vcf outputs can be collated together and genotyped into a .vcf file.

Internally, it will call GATK's HaplotypeCaller in GVCF mode. It takes care of parallelizing the computation across the
regions of interest in the genome (encoded into a BED file), and gathering the multiple work units into a single .g.vcf.

See docs on GATK4's [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.3.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php).
