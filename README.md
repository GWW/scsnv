```
scmap index -g genes.gtf -r genome.fa -l 100 out_prefix
    -g, --gtf
        GTF file (required)
    -r, --ref
        Genome fasta file (required)
    --retained-introns
        Keep retained introns, note alignments to these will be considered exonic
    -l, --min-length
        Minimum Transcript Length [100]
    -h, --help
        shows this help message

```

```
scmap count -k known_barcodes.txt -o barcodes.gz <fastq folder 1> <fastq folder 2> ...
    -k, --known-barcodes
        Known Barcode File
    -o, --out
        Barcode count output file
    -l, --library
        libary type (V2, V3)
    -h, --help
        shows this help message
```

```
scmap map -i <transcript index prefix> -g <genome bwa index> -b <barcode prefix> -o <out prefix> <fastq folder 1> ... <fastq folder N>
    -h, --help
        shows this help message
    -i, --index
        Transcript Index
    -g, --genome
        Genome BWA mem Index (optional)
    -b, --barcodes
        Barcode count prefix
    -t, --threads
        Number of threads
    -o, --output
        Output prefix
    --ei-ratio
        Minimum ratio of exonic / intergenic bases to be considered a cDNA alignment (Default: 0.75)
    -l, --library
        libary type (V2)
    --no-bam
        Disable writing the sorted bam files of the countable (ie. uniquely mapped reads)
    --bam-tmp
        Temporary directory to store sorted bam files (Default: {out_prefix}_tmp
    --bam-thread
        Number of output reads to buffer for each thread (Default: 50000)
    --bam-file
        Number of output reads per file (Default: 5000000)
    --bam-write
        Number of writer threads to use when emitting sorted bam files (Default 1)
```

```

scmap quant -i <transcript index prefix> -o <out prefix> <map prefix 1> ... <map prefix N>
    -h, --help
        shows this help message
    -i, --index
        Transcript Index
    -t, --threads
        Number of threads
    --no-bam
        Disable Writing data required to correct the UMI of bam files
    -c, --count-groups
        Gene Groups for cell quantification
    -m, --min-molecules
        Minimum number of barcode cDNA molecules
    -l, --library
        libary type (V2)
    -o, --output
        Output file prefix, ie. sample/summary
```
