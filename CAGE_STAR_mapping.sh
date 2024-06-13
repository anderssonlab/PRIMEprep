#!/bin/bash

DATE=$(date +%Y%m%d)

#CAGE mapping with STAR 

#Optional steps
  #Correction of mapping bias using WASP mapping: Requires samples specific VCF file with prefix matching fastq file in WASP directory
  
  #Filtering of rDNA reads using rDNAdust
  
  #Correction of unmatched G additions on 5' end of mapped reads


#The ouputs are bam, bed, ctss.bed, bigwwig and QC related files


#Load required modules
module load STAR
module load Preseq
module load fastp
module load perl
module load openjdk
module load fastqc
module load samtools
module load bedtools


##Default_parameter
FASTQ=""
GENOME="hg38"
FIRST_BASE=3 
THREADS=6
GENOME_PATH="/projects/ralab/data/genome"
SCRIPTDIR="/projects/ralab/apps"
DUSTFILE="U13369.1"
OUT=$(pwd)
G_CORRECT=true  #Default: remove added Gs and create bed file  
FILTER=false    #Default: do not filter rDNA
USE_VCF=false   #Default is false, enabled when VCF_PATH is given
VCF_PATH="" 



#Help_message
help_message="\
Options:\n\
\t-h\tThis help message.\n\
\t-f\t[STRING]\tFastq.gz file\t\t\t[required]\n\
\t-g\t[STRING]\tReference genome.\t\t[default=${GENOME}]\n\
\t-b\t[INTEGER]\tNumber of trimmed bases.\t[default=${FIRST_BASE}]\n\
\t-t\t[INTEGER]\tNumber of threads used.\t\t[default=${THREADS}]\n\
\t-d\t[STRING]\trRNA blacklist.\t\t\t[default=${DUSTFILE}]\n\
\t-o\t[STRING]\tOutput directory.\t\t[default=${OUT}]\n\
\t-i\t[BOOL]\t\trDNA filtering.\t\t[default=${FILTER}]\n\
\t-a\t[BOOL]\t\tAdd G and create bed file.\t[default=${G_CORRECT}]\n"


##Read command-line options
while getopts 'f:g:b:q:s:t:d:o:i:a:v:h' opt; do
  case ${opt} in
    f)
      FASTQ=${OPTARG};;
    g)
      GENOME=${OPTARG};;
    b)
      FIRST_BASE=${OPTARG};;
    t)
      THREADS=${OPTARG};;
    d)
      DUSTFILE=${OPTARG};;
    o)
      OUT=${OPTARG};;
    i)
      FILTER=${OPTARG};;
    v)
      USE_VCF=true
      VCF_PATH=${OPTARG};;
    a)
      G_CORRECT=${OPTARG};;
    h)
      echo -e ${help_message}
      exit 1;;
    \?)
      echo -e "Invalid option: -${OPTARG}" >&2
      exit 1;;
  esac
done
printf ">Running:\t\t\t\t%s\n" ${FASTQ}



#Create output subdirectories
if [ ! -d ${OUT} ]; then mkdir -p ${OUT}; chmod -R g+ws ${OUT}; fi
if [ ! -d ${OUT}/bam_files ]; then mkdir -p ${OUT}/bam_files; chmod -R g+ws ${OUT}/bam_files; fi
if [ ! -d ${OUT}/bed_files ]; then mkdir -p ${OUT}/bed_files; chmod -R g+ws ${OUT}/bed_files; fi
if [ ! -d ${OUT}/bw_files ]; then mkdir -p ${OUT}/bw_files; chmod -R g+ws ${OUT}/bw_files; fi
if [ ! -d ${OUT}}/QC ]; then mkdir -p ${OUT}/QC; chmod -R g+ws ${OUT}/QC; fi


##check_param
if [ '${FASTQ}' = "" ]; then \
		echo "Unable to find fastq file: ${FASTQ}"; \
		echo -e ${help_message}; \
		exit 1; \
	fi
if [ '${FASTQ}' != "" ] && [ ! -f $FASTQ ]; then \
		echo "No such file or directory: ${FASTQ}"; \
		echo -e ${help_message}; \
		exit 1; \
	fi
if [ ${FIRST_BASE} = "" ]; then \
		echo "Badly set first base: ${FIRST_BASE}"; \
		echo -e ${help_message}; \
		exit 1; \
	fi
if [ ${GENOME} = "" ]; then \
		echo "Badly set reference genome: ${GENOME}"; \
		echo -e ${help_message}; \
		exit 1; \
	fi
if [ ${GENOME} != "" ] && [ ! -f "${GENOME_PATH}/${GENOME}/STAR/SA" ]; then \
		echo "No reference genome files found: ${GENOME}"; \
		echo -e ${help_message}; \
		exit 1; \
	fi
if [ ! -f "${GENOME_PATH}/${DUSTFILE}/fasta/${DUSTFILE}.fa" ]; then \
		echo "No blacklist file found: ${DUSTFILE}"; \
		echo -e ${help_message}; \
		exit 1; \
	fi
if [ ${GENOME} != "hg38" ] && [ ${DUSTFILE} = "${GENOME_PATH}/${DUSTFILE}/fasta/${DUSTFILE}.fa" ]; then \
		echo "Please specify the blacklist file path for ${GENOME}"; \
		echo -e ${help_message}; \
		exit 1; \
	fi


##get_prefix_of_fastq
PREFIXRAW=$( basename $FASTQ | sed -e "s/.fastq.*//" )
PREFIX="${PREFIXRAW}_STAR"


##log_parameters
param=(\
       "DATE__${DATE}" \
       "WorkingDir__${OUT}" \
       "Threads__${THREADS}" \
       "File__${FASTQ}" \
       "Prefix__${PREFIX}" \
       "Genome__${GENOME}" \
       "TrimmedBp__${FIRST_BASE}"  \
       "GenomeFasta__${GENOME_PATH}" \
       "Blacklist__${DUSTFILE}" 
)
printf "Parameters\n" \
>> ${OUT}/QC/${PREFIX}_parameter.log
printf "\055%.s" {1..10} \
>> ${OUT}/QC/${PREFIX}_parameter.log
printf "\n%s" ${param[@]} \
>> ${OUT}/QC/${PREFIX}_parameter.log


##create_tmp_directory
TEMPDIR=$( mktemp -d -p /projects/ralab/scratch/ ) #Remove -p flag?
echo "tmp directory = ${TEMPDIR}"

printf "\nTempDir__${TEMPDIR}\n" \
>> ${OUT}/QC/${PREFIX}_parameter.log


##first_QC
printf ">Quality check before filtering:\t%s\n" ${PREFIXRAW}
fastqc \
  -o ${OUT}/QC/ \
  ${FASTQ} \
  >> ${OUT}/QC/${PREFIX}_parameter.log 2>&1

##trim_reads and filter reads
printf ">Trim and filter  reads:\t\t\t\t%s\n" ${PREFIXRAW}

fastp \
  -i ${FASTQ} \
  -o ${TEMPDIR}/${PREFIX}_trimmed.fastq.gz \
  --trim_front1 ${FIRST_BASE} \
  --length_required 30 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 50 \
  --thread ${THREADS}  \
  --json ${OUT}/QC/${PREFIX}.fastp.json \
  --html ${OUT}/QC/${PREFIX}.fastp.html \
  &> ${OUT}/QC/${PREFIX}_fastp.log


printf "\nFilteredFile__${TEMPDIR}/${PREFIX}_trimmed.fastq.gz\n" \
  >> ${OUT}/QC/${PREFIX}_parameter.log


##rDNA filtering
if ${FILTER}; then \
		printf ">Filter reads:\t\t\t\t%s\n" ${PREFIXRAW}; \
		gunzip -c ${TEMPDIR}/${PREFIX}_trimmed.fastq.gz \
				| ${SCRIPTDIR}/bin/rRNAdust \
					${GENOME_PATH}/${DUSTFILE}/fasta/${DUSTFILE}.fa \
					-t ${THREADS} \
					2>> ${OUT}/QC/${PREFIX}_parameter.log \
				|	gzip > ${TEMPDIR}/${PREFIX}_filtered.fastq.gz

		printf "\nFilteredFile__${TEMPDIR}/${PREFIX}_filtered.fastq.gz\n" \
			>> ${OUT}/QC/${PREFIX}_parameter.log
	else \
		printf ">Skip filtering reads.\n"
		mv ${TEMPDIR}/${PREFIX}_trimmed.fastq.gz ${TEMPDIR}/${PREFIX}_filtered.fastq.gz 
		printf "\nUnfilteredFile__${TEMPDIR}/${PREFIX}_filtered.fastq.gz\n" \
		 >> ${OUT}/QC/${PREFIX}_parameter.log 
	fi


##second_QC
printf ">Quality check after filtering:\t\t%s\n" ${PREFIXRAW}
fastqc \
  -o ${OUT}/QC/ \
  ${TEMPDIR}/${PREFIX}_filtered.fastq.gz \
  >> ${OUT}/QC/${PREFIX}_parameter.log 2>&1
  

##Run mapping and sorting
printf ">Map reads:\t\t\t\t%s\n" ${PREFIXRAW}
STAR_CMD=(
  "STAR"
  "--readFilesType Fastx"
  "--alignEndsType Extend5pOfRead1"
  "--readFilesIn ${TEMPDIR}/${PREFIX}_filtered.fastq.gz"
  "--readFilesCommand gunzip -c"
  "--readQualityScoreBase 33"
  "--outFilterMultimapNmax 1"
  "--outFilterMultimapScoreRange 1"
  "--genomeDir ${GENOME_PATH}/${GENOME}/STAR"
  "--runThreadN ${THREADS}"
  "--outFileNamePrefix ${OUT}/bam_files/${PREFIX}_filtered."
  "--outSAMtype BAM SortedByCoordinate"
  "--outBAMcompression 1"
  "--outBAMsortingThreadN ${THREADS}"
)

if ${USE_VCF} && [ -n "${VCF_PATH}" ]; then
  STAR_CMD+=("--varVCFfile ${VCF_PATH}/${PREFIXRAW}.vcf")
  STAR_CMD+=("--waspOutputMode SAMtag")
  STAR_CMD+=("--outSAMattributes vA vG vW")
else
  STAR_CMD+=("--outSAMattributes NH HI AS MD nM")
fi

"${STAR_CMD[@]}" >> ${OUT}/QC/${PREFIX}_parameter.log 2>&1


#Remove temporary mapping files
rm -r ${OUT}/bam_files/${PREFIX}_filtered._STARtmp

printf "\nBamFile__${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam" \
		>> ${OUT}/QC/${PREFIX}_parameter.log


##Inde bam file
samtools index \
	${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam

printf "\nBaiFile__${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.bai" \
		>> ${OUT}/QC/${PREFIX}_parameter.log


##alignment_complexity
printf ">Calculate alignment complexity:\t%s\n" ${PREFIXRAW}
preseq c_curve \
	-output ${OUT}/QC/${PREFIX}_filtered_preseq.txt \
	-bam ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam

printf "\nPreseqFile__${OUT}/QC/${PREFIX}_filtered_preseq.txt" \
		>> ${OUT}/QC/${PREFIX}_parameter.log


##alignment_stats
printf ">Calculate alignment stats:\t\t%s\n" ${PREFIXRAW}
samtools stats \
	${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
	> ${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.stats
samtools flagstat \
	${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
	> ${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.fstats

printf "\nBamStats__${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.stats" \
		>> ${OUT}/QC/${PREFIX}_parameter.log
printf "\nFlagStats__${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.fstats\n" \
		>> ${OUT}/QC/${PREFIX}_parameter.log


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


#Clean up (can cleaned up in script)
mv ${OUT}/bam_files/${PREFIX}*.out ${OUT}/QC/
mv ${OUT}/bam_files/${PREFIX}*.tab ${OUT}/QC/


##remove_tmp_files
rm -rf ${TEMPDIR}

printf ">Finished:\t\t\t\t%s\n" ${PREFIXRAW}
