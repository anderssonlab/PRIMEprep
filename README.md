# CAGE

Preprocessing, QC and mapping of CAGE sequencing data. Produces data compatible with [PRIME](https://github.com/anderssonlab/PRIME).

## Overview

1. Quality check before filtering using FastQC
2. Trim and filter reads using `fastp`
   - Number of trimmed bases should match barcode length
3. rDNA filtering (optional)
   - Filter reads using rRNA blacklist with `rRNAdust`
4. Quality check after filtering using `FastQC`
5. Map reads using `STAR`
6. Optional: Using VCF for variant-aware mapping**
7. Index resulting BAM file using `samtools`
8. Calculate alignment complexity using `preseq`
9. Calculate alignment statistics using `samtools`
10. Optional but recommended: Remove unmatched G additions and create bed files
    - Identify and process reads with G additions on the 5' end using `samtools` and `bedtools`


## Output directories & content
- **QC**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contains all QC reports. These can be consolidated into a single report using [`MultiQC`](https://multiqc.info).
- **bam_files**&nbsp;&nbsp;&nbsp;&nbsp;Contains the mapped reads. The G correction step does not alter these files.
- **bed_files**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contains bed files. If G correction is performed, unmatched G's at the 5' end will be removed from these files.
- **bw_files**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contains bigWig files compatible with [`PRIME`](https://github.com/anderssonlab/PRIME).


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
STAR index in `${GENOME_PATH}/${GENOME}/STAR`
chrom size file in `${GENOME_PATH}/${GENOME}/${GENOME}.chrom.sizes`


## Parameters
```
-h               This help message.
-f [STRING]      Fastq.gz file(s).           [required]
-g [STRING]      Specify reference genome.   [required]
-b [INTEGER]     Number of trimmed bases.    [default = 3]
-t [INTEGER]     Number of threads used.     [default = 6]
-p [STRING]      Genome path.                [required]
-s [STRING]      Script directory path.      [required]
-o [STRING]      Output directory.           [default = "./"]
-d [STRING]      rRNA blacklist.             [optional, default = false, enables rRNA filtering]
-a [BOOL]        Correct G additions.        [optional, default = true]
-v [STRING]      VCF path.                   [optional, default = false, enables VCF usage]
```

## Example
```bash
./CAGE_STAR_mapping.sh -f [FASTQ_FILE] -g [GENOME] -b [TRIM_BASES] -t [THREADS] -d [DUSTFILE] -o [OUTPUT_DIR] -i [FILTER] -a [G_CORRECT] -v [VCF_PATH] -p [GENOME_PATH] -s [SCRIPTDIR]

```

## Citation

To be added
