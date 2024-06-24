# CAGE

Preprocessing, QC and mapping of CAGE sequencing data. Produces data compatible with [PRIME](https://github.com/anderssonlab/PRIME).

## Overview

**1.**&nbsp;&nbsp;&nbsp;&nbsp;Quality check before filtering using [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). \
**2.**&nbsp;&nbsp;&nbsp;&nbsp;Trimming and filtering reads using [`fastp`](https://github.com/OpenGene/fastp).
   - Number of trimmed bases should match barcode length.
**3.**&nbsp;&nbsp;&nbsp;&nbsp;rDNA filtering (optional, e.g., human rDNA [U13369.1](https://www.ncbi.nlm.nih.gov/nuccore/U13369.1)).
   - Filtering reads using rRNA blacklist with [`rRNAdust`](https://fantom.gsc.riken.jp/5/suppl/rRNAdust/).
**4.**&nbsp;&nbsp;&nbsp;&nbsp;Quality check after filtering using [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). \
**5.**&nbsp;&nbsp;&nbsp;&nbsp;Mapping reads using [`STAR`](https://github.com/alexdobin/STAR). \
**6.**&nbsp;&nbsp;&nbsp;&nbsp;Using VCF for variant-aware mapping (optional). \
**7.**&nbsp;&nbsp;&nbsp;&nbsp;Indexing the generated BAM file using [`samtools`](http://www.htslib.org). \
**8.**&nbsp;&nbsp;&nbsp;&nbsp;Calculating the alignment complexity using [`preseq`](https://preseq.readthedocs.io/en/latest/). \
**9.**&nbsp;&nbsp;&nbsp;&nbsp;Calculating alignment statistics using [`samtools`](http://www.htslib.org). \
**10.**&nbsp;&nbsp;&nbsp;&nbsp;Removing unmatched G additions (optional, recommended).
    - Identify and process reads with G additions on the 5' end using [`samtools`](http://www.htslib.org) and [`bedtools`](https://bedtools.readthedocs.io/en/latest/).


## Output directories & content
- **QC**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contains all QC reports. These can be consolidated into a single report using [`MultiQC`](https://multiqc.info).
- **bam_files**&nbsp;&nbsp;&nbsp;&nbsp;Contains the mapped reads. The G correction step does not alter these files.
- **bed_files**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contains bed files. If G correction is performed, unmatched G's at the 5' end will be removed from these files.
- **bw_files**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contains bigWig files compatible with [`PRIME`](https://github.com/anderssonlab/PRIME).


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

## Example for single-end data
```bash
./CAGE_STAR_mapping.sh \
   -f *_001.fastq.gz \
   -g hg38 \
   -t 3 \
   -p  /path/to/genome/hg38/STAR/ \
   -s /path/to/script/directory/ \
   -d /path/to/human/rDNA/dustfile/U13369.1.fa \
   -o . \
   -a true
```

## Example for paired-end data
```bash
./CAGE_STAR_mapping.sh \
   -f *R1_001.fastq.gz \
   -f *R2_001.fastq.gz \
   -g hg38 \
   -t 6 \
   -p  /path/to/genome/hg38/STAR/ \
   -s /path/to/script/directory/ \
   -d /path/to/human/rDNA/dustfile/U13369.1.fa
```

## Citation

To be added
