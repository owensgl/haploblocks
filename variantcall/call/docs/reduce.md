[Back to start](./README.md)

Reduce Tool
==============

Reduces samples. Final step in the pipeline.
	
Usage: `Usage: reduce -r REF -L LIST OPERATION CALLER [opts]*`

Each file in `LIST` is a URL to a VCF file.

The output filename is based on the basename of the first input file.

The merging is done with one of multiple options: `freebayes` merge, `platypus` merge, `gatk3` merge, or `gatk4` merge. These merge functions are implemented in the [bed](./bed.md) module
