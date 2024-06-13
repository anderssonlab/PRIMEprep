# CAGE

Scripts for preprocesses, QC and mapping of CAGE sequencing data. Produced data compatible with [PRIME](https://github.com/anderssonlab/PRIME).

## Overview

1. Quality check before filtering using FastQC
2. Trim and filter reads using `fastp`
   - Number of trimmed bases should match barcode length
3. rDNA filtering (optional)
   - Filter reads using rRNA blacklist with `rRNAdust`
4. Quality check after filtering using FastQC
5. Map reads using STAR
6. Optional: Using VCF for variant-aware mapping**
7. Index resulting BAM file using `samtools`
8. Calculate alignment complexity using `preseq`
9. Calculate alignment statistics using `samtools`
10. Optional but recommended: Remove unmatched G additions and create bed files
    - Identify and process reads with G additions on the 5' end using `samtools` and `bedtools`


## Output
- `${OUT}/QC` contains all QC reports. These can be consolidated into a single report using MultiQC.
- `${OUT}/bam_files` contains the mapped reads. The G correction step does not alter these files.
- `${OUT}/bed_files` contains corresponding bed files. If G correction is performed, unmatched G's at the 5' end will be removed from these files.
- `${OUT}/bw_files` contains BigWig files compatible with [PRIME](https://github.com/anderssonlab/PRIME).




## Dependencies

In script directory:
- `rRNAdust`: A tool for rRNA filtering, found in the script directory's `bin` folder. (version 1.02)
- `bedGraphToBigWig`: A tool to convert bedGraph files to BigWig files, found in the script directory's `bin` folder. (version 4.0)

As executables:
- `STAR`: RNA-seq aligner for mapping reads to the genome. (version 2.7.3a)
- `Preseq`: For estimating library complexity. (version 2.0)
- `fastp`: A tool for quality control and filtering of sequencing data. (version 0.23.4)
- `perl`: Required for running various scripts. (version 5.38.0)
- `openjdk`: Java runtime environment. (version 20.0.0)
- `fastqc`: Quality control tool for high throughput sequence data. (version 0.12.1)
- `samtools`: Tools for manipulating next-generation sequencing data in BAM format. (version 1.20.0)
- `bedtools`: A powerful toolset for genome arithmetic. (version 2.31.0)

Files:
STAR index in ${GENOME_PATH}/${GENOME}/STAR
chrom size file in ${GENOME_PATH}/${GENOME}/${GENOME}.chrom.sizes


## Parameters

- **f [STRING]**: Fastq.gz file (required).
- **g [STRING]**: Reference genome (default=${GENOME}).
- **b [INTEGER]**: Number of trimmed bases (default=${FIRST_BASE}).
- **t [INTEGER]**: Number of threads used (default=${THREADS}).
- **d [STRING]**: rRNA blacklist (default=${DUSTFILE}).
- **o [STRING]**: Output directory (default=${OUT}).
- **i [BOOL]**: rDNA filtering (default=${FILTER}).
- **a [BOOL]**: Correct G additions and create bed files (default=${G_CORRECT}).
- **v [STRING]**: VCF path (enables VCF usage).
- **p [STRING]**: Genome path (required).
- **s [STRING]**: Script directory path (required).
- **h**: This help message.

## Example
```bash
./CAGE_STAR_mapping.sh -f [FASTQ_FILE] -g [GENOME] -b [TRIM_BASES] -t [THREADS] -d [DUSTFILE] -o [OUTPUT_DIR] -i [FILTER] -a [G_CORRECT] -v [VCF_PATH] -p [GENOME_PATH] -s [SCRIPTDIR]


## Citation

Update with PRIME citation
