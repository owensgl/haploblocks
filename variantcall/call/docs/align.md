[Back to start](./README.md)

Alignment Tool
==============

Performs alignment. 

`Usage: align JOBFILE REFGENOME`

Connects tools:

* ngm NextGenMap

  used to align sam files

* trimmomatic

  used to trim low quality samples in the lossy version of the aligner.

* picard
* bamutils
* sambamba, samtools

  the default is samtools, but it can use sambamba instead

  used to convert sam files to bam files, sort bam files, and concatenate bam files

  

Jobfile
-------

A "jobfile" is a json file specifying multiple parameters.

* name
* sample1 (path to sample), hash1 (path to file containing md5sum of sample1)
* sample2 (path to sample), hash2 (path to file containing md5sum of sample2)


Pipeline
---------

* lossy pipeline

  1. trim low quality samples. 

  1. align paired and unpaired files from trimmomatic, with ngm. produces 3 files: paired, unpaired1, unpaired2

  1. convert sam to bam (samtools)

  1. concatenate bams (samtools)

  1. sort concatenated bam file (samtools sort)

  1. if culling is enabled, index the bam file (samtools index)