[Back to start](./README.md)

Bed Library
==============

Tools for merging vcf files through concatenation. Abstraction over:

* `bcftools concat X Y Z -o OUTPT` (freebayes)
* `gatk3 CatVariants` (gatk3)
* `gatk4 GatherVcfs -I X -I Y -I Z -O OUTPUT` (gatk4)
* `platypus mergeVariants --intervalFiles=X,Y,Z --output=OUTPUT` (platypus)

Tools for naming and defining interval files.

* FIXME

