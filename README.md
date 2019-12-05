## scSNV

scSNV is an alternative to Cell Ranger + velocyto for the alignment, quantification and SNV calling of 10X Single Cell RNA-seq data.
Example python scripts to annotate SNVs, identify barcode tags that represent cells and generate the data from the scSNV manuscript is available at http://github.com/GWW/scsnvpy/.

scSNV requires the HDF5 C and C++ Libraries for compilation.

Building instructions:


The scSNV is separated into 6 steps:

1. [index](#indexing) -- Generating a reusable genome index
1. [count](#count) -- Counting the barcode tags
1. [map](#map) -- Mapping the alignment tags and correcting the barcodes
1. [quant](#quantify) -- Deduplicating and quantifying the UMI tags
1. [collapse](#collapse) -- Collapsing the sorted bam file fragments
1. [pileup](#pileup) -- Piling up the collapsed bam file

### Indexing
### Count
### Map
### Quantify
### Collapse
### Pileup

### Example Commands

### Future Tasks

1. Parallelizing the scSNV mapping step on a computational cluster:

