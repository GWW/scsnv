## scSNV

scSNV is an alternative to Cell Ranger + velocyto for the alignment, quantification and SNV calling of 10X Single Cell RNA-seq data. Currently, scSNV supports 10X 5-prime and 3-prime 10X version 2 and 3 libraries.

Example python scripts to annotate SNVs, identify barcode tags that represent cells and generate the data from the scSNV manuscript is available at http://github.com/GWW/scsnvpy/.

For reproducibility the version of Cell Ranger used in our manuscript comparisons is available at [here](https://github.com/GWW/cellranger_211_mirror)


scSNV requires the HDF5 C and C++ Libraries for compilation.

## Building instructions:

```bash
git clone --recurse-submodules https://github.com/GWW/scsnv.git
cd scsnv
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

```bash
cd scsnvpy
python setup.py install
```

The only executable needed to use scsnv is located in scsnv/build/src/scsnv
The scsnvmisc python executable will be installed

## Directory Setup

scSNV expects the fastq files to be generated with the Cell Ranger mkfastq command.  For example:
```
./run1/Library1_S1_L001_R1_001.fastq.gz
./run1/Library1_S1_L001_R2_001.fastq.gz
./run1/Library1_S1_L001_I1_001.fastq.gz
```
We recommend organizing your data as follows and assume this organization for the manual
```
sample/run1/<fastq_files>
```
Where sample is the name you wish to use for your sample and run1 is the name you wish to use for the folder containing your fastq files

## Example alignment and pileup for one sample

This assumes a 10X V2 3'-library using the 737K barcode set from 10X genomics.  For V3 libraries you need to use the 3M set.
This assumes you have already indexed the genome fasta file with BWA.

```bash
#Build the index, only required once
scsnv index -g genes.gtf -r genome.fa index_prefix

#Count the number of barcodes
scsnv count -o sample/barcode  -k scsnv/data/737K-august-2016.txt -l V2 sample/run1

#Map the reads, quantify gene expression, and write the sorted mRNA-tag alignments
#See below regarding the gene_group file
scsnv map -l V2 -i index_prefix -g bwa_genome_index -b sample/barcode -t 24 --bam-write 4 -q 4 -c index_prefix/gene_groups.txt -o sample/ sample/run1

#Collapse the mRNA-tags into collapsed molecules
scsnv collapse -l V2 -r genome_fasta -i index_prefix -o sample/ --threads 4 --bam-write 8 -b sample/barcode_counts.txt.gz sample/merged.bam

#Find the optimal number of cells
scsnvmisc cells -o sample sample/summary.h5

#Pileup the reads from the collapsed molecules using a list of passed barcodes
scsnv pileup -l V2 -i index_prefix -r genome_fasta -o sample/pileup -p ./sample/passed_barcodes.txt.gz -t 4 -x 4 ./sample/collapsed.bam

#The pileup can be annotated and bi-allelic strand-specific SNVs can be called using the scsnvpy annotate command (See Below)

#Quantify SNV co-expression and collapsed molecule lengths
scsnv snvcounts -t -s sample/pileup_passed_snvs.txt.gz -b sample/passed_barcodes.txt.gz -i index_prefix -o sample/snv -l V2 sample/collapsed.bam
```

Note:

The pileup command will also work with the Cell Ranger BAM file or STAR Solo BAM file. It requires the -c flag to be set and a passed_barcodes.txt.gz file with the following format:

```
barcode
GCAATCAGTACTTGAC
AAATGCCGTCGCCATG
...
etc
```
The barcodes should not contain a -1 in them.  For the bam files these will automatically be removed from the CB tag. The annotate command will work on the output from these files.


##### Output Files:
| File        | Contents      |
| ---------------|:--------------|
| sample/pileup\_barcode_matrices.h5 | Pileup data matrices (See Below) |
| sample/pileup.txt.gz | Summary data for each SNV that passed the filtering criteria |

##### Pileup H5 File:
| Group        | Value      |
| ---------------|:--------------|
| barcodes | Cellular barcodes (used by cellular barcode indexes) |
| pos | Zero-based SNV position (used by SNV row indexes) |
| tid | SNV chromosome Name Index (used by SNV row indexes) |
| refs | SNV chromosome Names (used by SNV row indexes) |
| base_A, base_C, base_G, base_C | Sparse base count matrices 1. **barcode_ids**: Cell Barcode Indexes, 2. **snps**: SNV Row Indexes, 3. **plus**: Plus strand counts, 4. **minus**: Minus strand counts. |
| barcode_rates | Per cellular barcode base annotation counts for each nucleotide and strand |
| coverage | Sparse matrix of the base coverage rates for the plus and minus strand 1. chroms: Chromosome names, 2. **tid**: Chromosome ids 3. **pos**: Zero-based position 4. **plus**: Plus strand coverage, 5. **minus**: Minus strand coverage|

#### Annotate
##### Command:
```bash
scsnvmisc annotate -r ./edits/repeat_masker.txt.gz -d 1000GENOMES.txt.gz -e REDIportal.bed.gz sample/pileup
```
##### Arguments:
| Option        | Argument      | Function | Required |
| ---------------|:--------------|:---------|:---------|
| -e, --edits | path | BED file of RNA edits from REDIportal | No |
| -r, --repeats | path | Path to UCSC repeat masker annotations | No |
| -d, --g1000 | path | Path to 1000 Genomes txt file (see below to make) | No |
| -c, --capture | path | Exome capture region bed file, uses the location of these + flank | No |
| -f, --flank | path | Flank for capture regions (Default 100 bp) | No |
| -v, --vcf | path | VCF file to annotate with (uses the name field below), can be repeated for multiple VCF files | No |
| -n, --name | path | Name key to use if a SNV is found in this VCF file (one for each VCF file) | No |
| -u, --purity | float | Proportion of non-reference counts that must match a single allele (Default: 0.95 | No |
| -s, --strand-purity | float | Proportion of reads that must be derived from the same strand to call a stranded SNP (Default: 0.9) | No |

##### 1000 Genomes File Generation:
```bash
wget ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz
wget ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz.csi
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi

bcftools view -m2 -M2 -v snps -i "MAF>=0.01" -o 1000GENOMES-filtered.vcf.gz -O z 1000GENOMES-phase_3.vcf.gz
bcftools +fill-tags ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz -- -t AF | bcftools view -m2 -M2 -v snps -i "MAF>=0.01" -o ALL.chrMT-filtered.vcf.gz

zgrep --no-filename -E "^(#)" -v ALL.chrMT-filtered.vcf.gz | gzip > ALL.chrMT-filtered-no-header.vcf.gz
zgrep --no-filename -E "^(#)\1+" -v 1000GENOMES-filtered.vcf.gz ALL.chrMT-filtered-no-header.vcf.gz | cut -f 1,2,4,5 | gzip > 1000GENOMES.txt.gz
```

##### Output Files:
| File        | Contents      |
| ---------------|:--------------|
| sample/pileup_annotated.h5 | Serialized python flammkuchen file with annotated SNV data and per cell counts |
| sample/pileup_barcodes.txt.gz | Barcode molecule counts and base coverage |
| sample/pileup.txt.gz | Raw pileup calls |
| sample/pileup_barcode_matrices.h5 | Raw pileup barcode count matrices |
| sample/cells.h5 | Serialized python flammkuchen file with high quality cell and expression data |
| sample/summary.h5 | Serialized HDF5 file with raw cell barcode statistics and expression data data|
| sample/snv_map.txt.gz | SNV co-expression map |
| sample/snv_base_counts.txt.gz | histograms of collapsed base lengths based on the number of reads and pcr duplicates in the molecules |
| sample/snv_reads.txt.gz | SNV co-expression molecule data, the SNV_idx is defined in the snv_map file or the passed_snvs file | 
| sample/barcode_bases.txt.gz | Maximum collapsed read lengths per barcode |


The pileup file can be annotated with various databases using the scsnvmisc annotate command.  The serialized flammkuchen files can be viewed using ddls after flammkuchen is installed. (pip install flammkuchen)

#### Merge -- Merge the annotated output from multiple alignment tools to calculate accuracy using a matched whole genome or exome sequence
##### Command:
```bash
scsnvmisc merge -n scsnv -n cellranger --vcf <snvs>.vcf -o sample/merged scsnv/pileup cellranger/pileup
```

A -n argument is required for each alignment method. 

#### accuracy -- Calculate barcode counts and coverage information for the SNVs in the file outputted by the merge command
##### Command:
```bash
scsnv accuracy -n genome -a G -b sample/passed_barcodes.txt.gz -r genome_fasta -i index_prefix -o sample/accuracy_genome.txt.gz -s sample/merged.txt.gz -l V2 WGS.bam
scsnv accuracy -n scsnv -a S -b sample/passed_barcodes.txt.gz -r genome_fasta -i index_prefix -o sample/accuracy_scsnv.txt.gz -s sample/merged.txt.gz -l V2 collapsed.bam
scsnv accuracy -n scsnv_merged -a D -b sample/passed_barcodes.txt.gz -r genome_fasta -i index_prefix -o sample/accuracy_scsnv_merged.txt.gz -s sample/merged.txt.gz -l V2 collapsed.bam

scsnv accuracy -n cellranger -a C -b sample/passed_barcodes.txt.gz -r genome_fasta -i index_prefix -o sample/accuracy_cellranger.txt.gz -s sample/merged.txt.gz -l V2 cellranger/outs/possorted_genome_bam.bam
```
#### results -- Merge the results from multiple scsnv accuracy runs

scsnvmisc results -o results.txt.gz sample/accuracy

The output file will contain merged counts from all of the files as well as some useful annotations from the merged_snvs.txt.gz file that was generated from the merge command. This command will read all files with the prefix ```sample/accuracy_*.txt.gz``` and merge them together



More documentation on the output formats and optional commands forthcoming 

