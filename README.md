# PRIMEprep

Shell sript to preprocess single- or paired-end CAGE-sequencing data as generated by the SLIC-CAGE protocol ([Cvetesic et al., 2018](https://genome.cshlp.org/content/28/12/1943.long)). The script requires fastq file(s) as input, trims off barcode sequences from the reads before mapping them to the corresponding genome with [`STAR`](https://github.com/alexdobin/STAR) and finally generates bigWig files that serve as input for downstream analysis with [`PRIME`](https://github.com/anderssonlab/PRIME) or [`CAGEfightR`](https://github.com/MalteThodberg/CAGEfightR). Optional pre-processing steps include 1) filtering of reads derived from abundant, ribosomal RNAs via [`rRNAdust`](https://fantom.gsc.riken.jp/5/suppl/rRNAdust/) as [recommended](https://fantom.gsc.riken.jp/5/sstar/Protocols:rRNAdust) by the FANTOM5 consortium, 2) STAR mapping supported by personal variants that overlap the aigned reads for which vcf files would be required and 3) removal of Gs at the read's 5' end that don't match the genome sequence and are due to the 3' overhang produced by the Superscript reverse transcriptase during cDNA first-strand synthesis. Please see **Overview of data processing steps** for further details on all pre-processing steps.

## Overview of data processing steps

**1.**&nbsp;&nbsp;&nbsp;&nbsp;Quality check before filtering using [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). \
**2.**&nbsp;&nbsp;&nbsp;&nbsp;Trimming and filtering reads using [`fastp`](https://github.com/OpenGene/fastp). (Number of trimmed bases should match barcode length.) \
**3.**&nbsp;&nbsp;&nbsp;&nbsp;rRNA filtering with [`rRNAdust`](https://fantom.gsc.riken.jp/5/suppl/rRNAdust/) & an rRNA blacklist (optional, e.g., human rDNA [U13369.1](https://www.ncbi.nlm.nih.gov/nuccore/U13369.1)). \
**4.**&nbsp;&nbsp;&nbsp;&nbsp;Quality check after filtering using [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). \
**5.**&nbsp;&nbsp;&nbsp;&nbsp;Mapping reads using [`STAR`](https://github.com/alexdobin/STAR). \
**6.**&nbsp;&nbsp;&nbsp;&nbsp;Using VCF for variant-aware mapping (optional). \
**7.**&nbsp;&nbsp;&nbsp;&nbsp;Indexing the generated BAM file using [`samtools`](http://www.htslib.org). \
**8.**&nbsp;&nbsp;&nbsp;&nbsp;Calculating the alignment complexity using [`preseq`](https://preseq.readthedocs.io/en/latest/). \
**9.**&nbsp;&nbsp;&nbsp;&nbsp;Calculating alignment statistics using [`samtools`](http://www.htslib.org). \
**10.**&nbsp;&nbsp;&nbsp;Removing unmatched G additions (optional, recommended). \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Identify and process reads with G additions on the 5' end using [`samtools`](http://www.htslib.org) and [`bedtools`](https://bedtools.readthedocs.io/en/latest/).

## Output directories & content

- **QC**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contains all QC reports. These can be consolidated into a single report using [`MultiQC`](https://multiqc.info).
- **bam_files**&nbsp;&nbsp;&nbsp;&nbsp;Contains the mapped reads. The G correction step does not alter these files.
- **bed_files**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contains bed files. If G correction is performed, unmatched G's at the 5' end will be removed.
- **bw_files**&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Contains bigWig files compatible with [`PRIME`](https://github.com/anderssonlab/PRIME).

## Dependencies

**[1.] The following software must be installed and permissions granted to execute them:**
- **[`rRNAdust`](https://fantom.gsc.riken.jp/5/suppl/rRNAdust/)**(v1.02)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Tool to filter abundant rRNA species.
- **[`bedGraphToBigWig`](https://www.encodeproject.org/software/bedgraphtobigwig/)**(v4.0)&nbsp;&nbsp;Allows converting bedGraph to bigWig files.
- **[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**&nbsp;(v0.12.1)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Quality control tool for high throughput sequence data. 
- **[`fastp`](https://github.com/OpenGene/fastp)**&nbsp;(v0.23.4)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A tool for quality control and filtering of sequencing data. 
- **[`STAR`](https://github.com/alexdobin/STAR)**&nbsp;(v2.7.3a)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;RNA-seq aligner for mapping reads to the genome. 
- **[`samtools`](http://www.htslib.org)**&nbsp;(v1.20.0)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Tools for manipulating next-generation sequencing data in BAM format.
- **[`preseq`](https://preseq.readthedocs.io/en/latest/)**&nbsp;(v2.0)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For estimating library complexity.
- **[`Perl`](https://www.perl.org/get.html)**&nbsp;(v5.38.0)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Required for running various scripts. 
- **[`openjdk`](https://openjdk.org)**&nbsp;(v20.0.0)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Java runtime environment.
- **[`bedtools`](https://bedtools.readthedocs.io/en/latest/)**&nbsp;(v2.31.0))&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Set of tools for genome arithmetic.

**[2.] Paths for the following, pre-computed files need to be provided:**
- The path to the STAR genome index corresponding to the species the datasets are derived from (Provided by `-g`, see **Parameters**). For instance, for the hg38 genome assembly, the corresponding human STAR genome index is prepared from the corresponding [fasta](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz) and [gtf](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz) files:
```
STAR   \
   --runMode genomeGenerate \
   --runThreadN 3 \
   --genomeDir hg38/STAR \
   --genomeFastaFiles fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
   --sjdbGTFfile ./gencode.v43.chr_patch_hapl_scaff.annotation.gtf
```
- The path to the fasta file including the rDNA sequence corresponding to the reference genome, e.g., [U13369.1](https://www.ncbi.nlm.nih.gov/nuccore/U13369.1) for the human genme assembly hg38.

## Parameters
```
-h              Help message.
-f   [STRING]   Fastq.gz file(s).                         [required]
-g   [STRING]   Path to STAR index of reference genome.   [required]
-b   [INTEGER]  Number of trimmed bases.                  [default=3]
-t   [INTEGER]  Number of threads used.                   [default=6]
-o   [STRING]   Output directory.                         [default=.]
-d   [STRING]   rDNA blacklist fasta file.                [default=""]
-a   [BOOL]     Correct G additions.                      [default=true]
-v   [BOOL]     VCF path.                                 [default=""]
```

## Example for single-end data with G-correction
```bash
./CAGE_STAR_mapping.sh \
   -f *_001.fastq.gz \
   -g /path/to/genome/STAR/ \
   -t 3 \
   -d /path/to/rDNA/blacklist/<blacklist>.fa \
   -o . \
   -a true
```

## Example for paired-end data without G-correction
```bash
./CAGE_STAR_mapping.sh \
   -f *R1_001.fastq.gz \
   -f *R2_001.fastq.gz \
   -g /path/to/genome/STAR/ \
   -t 6 \
   -d /path/to/rDNA/blacklist/<blacklist>.fa
   -a false
```

## Citation

To be added.
