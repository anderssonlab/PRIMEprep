#!/bin/bash

# Set -e to exit immediately if a command exits with a non-zero status
set -e

DATE=$(date +%Y%m%d)

# Default parameters
FASTQ=""
GENOME="hg38"
FIRST_BASE=3 
THREADS=6
GENOME_PATH=""
SCRIPTDIR=""
DUSTFILE="U13369.1"
OUT=$(pwd)
G_CORRECT=true
FILTER=false
USE_VCF=false
VCF_PATH=""

# Help message
help_message="\
Options:
    -h                  This help message.
    -f [STRING]         Fastq.gz file (required).
    -g [STRING]         Reference genome (default=${GENOME}).
    -b [INTEGER]        Number of trimmed bases (default=${FIRST_BASE}).
    -t [INTEGER]        Number of threads used (default=${THREADS}).
    -d [STRING]         rRNA blacklist (default=${DUSTFILE}).
    -o [STRING]         Output directory (default=${OUT}).
    -i [BOOL]           rDNA filtering (default=${FILTER}).
    -a [BOOL]           Correct G additions and create bed files (default=true).
    -v [STRING]         VCF path (enables VCF usage).
    -p [STRING]         Genome path (required).
    -s [STRING]         Script directory path (required).
"

# Parse command-line options
while getopts 'f:g:b:t:d:o:i:a:v:p:s:h' opt; do
    case ${opt} in
        f) FASTQ=${OPTARG} ;;
        g) GENOME=${OPTARG} ;;
        b) FIRST_BASE=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        d) DUSTFILE=${OPTARG} ;;
        o) OUT=${OPTARG} ;;
        i) FILTER=${OPTARG} ;;
        v) USE_VCF=true; VCF_PATH=${OPTARG} ;;
        a) G_CORRECT=true ;;
        p) GENOME_PATH=${OPTARG} ;;
        s) SCRIPTDIR=${OPTARG} ;;
        h) echo -e ${help_message}; exit 0 ;;
        \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    esac
done

# Ensure required arguments are provided
if [[ -z ${FASTQ} ]]; then
    echo "Fastq file is required."
    echo -e ${help_message}
    exit 1
fi

if [[ ! -f ${FASTQ} ]]; then
    echo "No such file or directory: ${FASTQ}"
    exit 1
fi

if [[ -z ${GENOME_PATH} ]]; then
    echo "Genome path is required."
    echo -e ${help_message}
    exit 1
fi

if [[ -z ${SCRIPTDIR} ]]; then
    echo "Script directory path is required."
    echo -e ${help_message}
    exit 1
fi

# Create output subdirectories
mkdir -p ${OUT}/{bam_files,bed_files,bw_files,QC}
chmod -R g+ws ${OUT}/{bam_files,bed_files,bw_files,QC}

# Log parameters
PREFIXRAW=$(basename ${FASTQ} | sed -e "s/.fastq.*//")
PREFIX="${PREFIXRAW}_STAR"

log_parameters() {
    local log_file=${OUT}/QC/${PREFIX}_parameter.log
    local params=("DATE__${DATE}" "WorkingDir__${OUT}" "Threads__${THREADS}" "File__${FASTQ}" "Prefix__${PREFIX}" "Genome__${GENOME}" "TrimmedBp__${FIRST_BASE}" "GenomeFasta__${GENOME_PATH}" "Blacklist__${DUSTFILE}")

    printf "Parameters\n" > ${log_file}
    printf "\055%.s" {1..10} >> ${log_file}
    printf "\n%s" "${params[@]}" >> ${log_file}
}

log_parameters

# Create temporary directory within the output directory
TEMPDIR=$(mktemp -d -p ${OUT})
trap "rm -rf ${TEMPDIR}" EXIT
echo "Temporary directory: ${TEMPDIR}"

# Quality check before filtering
printf "> Quality check before filtering: %s\n" ${PREFIXRAW}
fastqc -o ${OUT}/QC/ ${FASTQ} >> ${OUT}/QC/${PREFIX}_parameter.log 2>&1

# Trim and filter reads
printf "> Trim and filter reads: %s\n" ${PREFIXRAW}
fastp -i ${FASTQ} -o ${TEMPDIR}/${PREFIX}_trimmed.fastq.gz --trim_front1 ${FIRST_BASE} --length_required 30 --qualified_quality_phred 20 --unqualified_percent_limit 50 --thread ${THREADS} --json ${OUT}/QC/${PREFIX}.fastp.json --html ${OUT}/QC/${PREFIX}.fastp.html &> ${OUT}/QC/${PREFIX}_fastp.log

# rDNA filtering
if ${FILTER}; then
    printf "> Filter reads: %s\n" ${PREFIXRAW}
    gunzip -c ${TEMPDIR}/${PREFIX}_trimmed.fastq.gz | ${SCRIPTDIR}/bin/rRNAdust ${GENOME_PATH}/${DUSTFILE}/fasta/${DUSTFILE}.fa -t ${THREADS} 2>> ${OUT}/QC/${PREFIX}_parameter.log | gzip > ${TEMPDIR}/${PREFIX}_filtered.fastq.gz
else
    printf "> Skip filtering reads.\n"
    mv ${TEMPDIR}/${PREFIX}_trimmed.fastq.gz ${TEMPDIR}/${PREFIX}_filtered.fastq.gz
fi

# Quality check after filtering
printf "> Quality check after filtering: %s\n" ${PREFIXRAW}
fastqc -o ${OUT}/QC/ ${TEMPDIR}/${PREFIX}_filtered.fastq.gz >> ${OUT}/QC/${PREFIX}_parameter.log 2>&1

# Run mapping and sorting
printf "> Map reads: %s\n" ${PREFIXRAW}
STAR_CMD=("STAR" "--readFilesType Fastx" "--alignEndsType Extend5pOfRead1" "--readFilesIn ${TEMPDIR}/${PREFIX}_filtered.fastq.gz" "--readFilesCommand gunzip -c" "--readQualityScoreBase 33" "--outFilterMultimapNmax 1" "--outFilterMultimapScoreRange 1" "--genomeDir ${GENOME_PATH}/${GENOME}/STAR" "--runThreadN ${THREADS}" "--outFileNamePrefix ${OUT}/bam_files/${PREFIX}_filtered." "--outSAMtype BAM SortedByCoordinate" "--outBAMcompression 1" "--outBAMsortingThreadN ${THREADS}")

if ${USE_VCF} && [[ -n ${VCF_PATH} ]]; then
    STAR_CMD+=("--varVCFfile ${VCF_PATH}/${PREFIXRAW}.vcf" "--waspOutputMode SAMtag" "--outSAMattributes vA vG vW")
else
    STAR_CMD+=("--outSAMattributes NH HI AS MD nM")
fi

"${STAR_CMD[@]}" >> ${OUT}/QC/${PREFIX}_parameter.log 2>&1

# Remove temporary mapping files
rm -r ${OUT}/bam_files/${PREFIX}_filtered._STARtmp

# Index bam file
samtools index ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam

# Calculate alignment complexity
printf "> Calculate alignment complexity: %s\n" ${PREFIXRAW}
preseq c_curve -output ${OUT}/QC/${PREFIX}_filtered_preseq.txt -bam ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam

# Calculate alignment stats
printf "> Calculate alignment stats: %s\n" ${PREFIXRAW}
samtools stats ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam > ${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.stats
samtools flagstat ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam > ${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.fstats

#Remove unmatched G addition and create bed file (if G_CORRECT is true)
if ${G_CORRECT}; then
  printf ">Removing unmatched G on 5' end and creating bedfiles:\t\t%s\n" ${PREFIXRAW}


  #Remove unmatched G addition on + strand
  samtools view -F 16 ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam | \
    awk 'BEGIN {FS="\t"} { if ($10 ~ /^GG.+/ && $0 ~ /MD:Z:0[A-Z]0[A-Z]/) print 2; else if ($10 ~ /^G.+/ && $0 ~ /MD:Z:0/) print 1; else print 0 }' \
    > ${TEMPDIR}/${PREFIX}.plus.counts

  samtools view -b -F 16 ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam | \
    bamToBed > ${TEMPDIR}/${PREFIX}.plus.bed

  paste ${TEMPDIR}/${PREFIX}.plus.bed ${TEMPDIR}/${PREFIX}.plus.counts \
    > ${TEMPDIR}/${PREFIX}.combined.plus.bed

  awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $2 + $7, $3, $4, $5, $6, $7}' ${TEMPDIR}/${PREFIX}.combined.plus.bed \
    > ${TEMPDIR}/${PREFIX}.final.plus.bed


  #Remove unmatched G addition on - strand
  samtools view -f 16 ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam | \
    awk 'BEGIN {FS="\t"} { if ($10 ~ /.+CC$/ && $0 ~ /MD:Z:.+[A-Z]0[A-Z]0/) print 2; else if ($10 ~ /.+C$/ && $0 ~ /MD:Z:.+[A-Z]0/) print 1; else print 0 }' \
    > ${TEMPDIR}/${PREFIX}.minus.counts

  samtools view -b -f 16 ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam | \
    bamToBed > ${TEMPDIR}/${PREFIX}.minus.bed

  paste ${TEMPDIR}/${PREFIX}.minus.bed ${TEMPDIR}/${PREFIX}.minus.counts \
    > ${TEMPDIR}/${PREFIX}.combined.minus.bed

  awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $2 , $3 - $7, $4, $5, $6, $7}' ${TEMPDIR}/${PREFIX}.combined.minus.bed \
    > ${TEMPDIR}/${PREFIX}.final.minus.bed


  #Merge +/- and sort final bed file
  cat ${TEMPDIR}/${PREFIX}.final.plus.bed ${TEMPDIR}/${PREFIX}.final.minus.bed > ${TEMPDIR}/${PREFIX}.final.bed

  sort -k 1,1 -k2,2n ${TEMPDIR}/${PREFIX}.final.bed > ${OUT}/bed_files/${PREFIX}.bed

  printf ">Finished removing unmatched G and creating bed file:\t%s\n" ${PREFIXRAW}
  
else

  bamToBed -i ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam > ${OUT}/bed_files/${PREFIX}.final.bed
  printf ">Skipped removing unmatched G:\t%s\n" ${PREFIXRAW}
  
fi


#Format CTSS bed and bigwig files
awk 'BEGIN{OFS="\t"}{if ($6 == "+") {print $1, $2, $2+1, $6} else {print $1, $3-1, $3, $6}}' ${OUT}/bed_files/${PREFIX}.bed | \
  sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2, $3, $4, $2 ":" $3 "-" $4 "," $5, $1, $5}' | \
  sort -k 1,1 -k 2,2n > ${OUT}/bed_files/${PREFIX}.ctss.bed


grep -F ",-" ${OUT}/bed_files/${PREFIX}.ctss.bed | cut -f 1,2,3,5 | sort -k1,1 -k 2,2n > ${TEMPDIR}/${PREFIX}.minus.bedgraph
grep -F ",+" ${OUT}/bed_files/${PREFIX}.ctss.bed | cut -f 1,2,3,5 | sort -k1,1 -k 2,2n > ${TEMPDIR}/${PREFIX}.plus.bedgraph


${SCRIPTDIR}/bin/bedGraphToBigWig ${TEMPDIR}/${PREFIX}.minus.bedgraph ${GENOME_PATH}/${GENOME}/chrom_size/${GENOME}.chrom.sizes ${OUT}/bw_files/${PREFIX}.minus.bw
${SCRIPTDIR}/bin/bedGraphToBigWig ${TEMPDIR}/${PREFIX}.plus.bedgraph ${GENOME_PATH}/${GENOME}/chrom_size/${GENOME}.chrom.sizes ${OUT}/bw_files/${PREFIX}.plus.bw

# Remove_tmp_files
rm -rf ${TEMPDIR}

# Completion message
echo "Script execution completed successfully."

exit 0
