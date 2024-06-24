#!/bin/bash

#CAGE mapping with STAR

#Optional steps
  #Correction of mapping bias using WASP mapping: Requires samples specific VCF file with prefix matching fastq file in WASP directory

  #Filtering of rDNA reads using rDNAdust

  #Correction of unmatched G additions on 5' end of mapped reads

#The ouput consists of bam, bed, ctss.bed, bigwwig and QC-related files.

# Set -e to exit immediately if a command exits with a non-zero status
set -e

DATE=$(date +%Y%m%d)

##Default_parameters
##------------------
GENOME=""
FIRST_BASE=3
THREADS=6
GENOME_PATH=""
SCRIPTDIR=""
OUT=$(pwd)

FILTER=false    #Default is false: do not filter rDNA
DUSTFILE=""     #e.g., U13369.1 for human rDNA complete repeating unit

G_CORRECT=true  #Default is true: remove added Gs and create bed file

USE_VCF=false   #Default is false: enabled when VCF_PATH is given
VCF_PATH=""

##Help_message
##------------
help_message="\
Options:\n\
\t-h\t\t\tThis help message.\n\
\t-f\t[STRING]\tFastq.gz file(s)\t\t\t[required]\n\
\t-g\t[STRING]\tReference genome.\t\t\t[default=${GENOME}]\n\
\t-b\t[INTEGER]\tNumber of trimmed bases.\t\t[default=${FIRST_BASE}]\n\
\t-t\t[INTEGER]\tNumber of threads used.\t\t\t[default=${THREADS}]\n\
\t-p\t[STRING]\tGenome path.\t\t\t\t[required]\n\
\t-s\t[STRING]\tScript directory path.\t\t\t[required]\n\
\t-o\t[STRING]\tOutput directory.\t\t\t[default=${OUT}]\n\
\t-d\t[STRING]\trRNA blacklist.\t\t\t\t[default=${DUSTFILE}]\n\
\t-a\t[BOOL]\t\tCorrect G additions.\t[default=${G_CORRECT}]\n
\t-v\t[BOOL]\t\tVCF path.\t\t\t\t[default=${VCF_PATH}]\n\
"

##Parse_utility_options
##---------------------
while getopts 'f:g:b:t:p:s:o:d:a:v:h' opt; do
  case ${opt} in
    f)
      FASTQ+=(${OPTARG})
      ;;
    g)
      GENOME=${OPTARG}
      ;;
    b)
      FIRST_BASE=${OPTARG}
      ;;
    t)
      THREADS=${OPTARG}
      ;;
    p)
      GENOME_PATH=${OPTARG}
      ;;
    s)
      SCRIPTDIR=${OPTARG}
      ;;
    o)
      OUT=${OPTARG}
      ;;
    d)
      FILTER=true
      DUSTFILE=${OPTARG}
      ;;
    a)
      G_CORRECT=true
      ;;
    v)
      USE_VCF=true
      VCF_PATH=${OPTARG}
      ;;
    h)
      echo -e ${help_message}
      exit 0
      ;;
    \?)
      echo -e "Invalid option: -${OPTARG}" >&2
      exit 1
      ;;
  esac
done

##Create_output_subdirectories
##----------------------------
mkdir -p ${OUT}/{bam_files,bed_files,bw_files,QC}
chmod -R g+ws ${OUT}/{bam_files,bed_files,bw_files,QC}

##Ensure_required_arguments_are_provided
##--------------------------------------
if [ ${#FASTQ[@]} -eq 0 ]; then
		echo "Unable to find fastq file(s): ${FASTQ[@]}"
		echo -e ${help_message}
		exit 1
	fi
if [ ${#FASTQ[@]} -gt 2 ]; then
		echo "More than two *.fastq files provided: ${#FASTQ[@]}"
		echo -e ${help_message}
		exit 1
	fi
for f in ${FASTQ[@]}; do
		if [ ! -f ${f} ]; then
			echo "No such file or directory: ${f}"
			echo -e ${help_message}
			exit 1
		fi
	done
if [ ${GENOME} = "" ]; then
		echo "Badly set reference genome: ${GENOME}"
		echo -e ${help_message}
		exit 1
	fi
if [ ${GENOME} != "" ] && [ ! -f "${GENOME_PATH}/${GENOME}/STAR/SA" ]; then
		echo "No reference genome files found: ${GENOME_PATH}"
		echo -e ${help_message}
		exit 1
	fi
if [ ${FIRST_BASE} = "" ]; then
		echo "Badly set first base: ${FIRST_BASE}"
		echo -e ${help_message}
		exit 1
	fi
if [[ -z ${SCRIPTDIR} ]]; then
    echo "Script directory path is required."
    echo -e ${help_message}
    exit 1
  fi
if [ ! -f ${DUSTFILE} ]; then
		echo "No blacklist file found: ${DUSTFILE}"
		echo -e ${help_message}
		exit 1
	fi

##Get_prefix_of_fastqs
##--------------------
if [ ${#FASTQ[@]} -gt 1 ]; then \
		# echo -e "Number of FASTQ-files:\t${#FASTQ[@]}"; \
		PREFIXS=(); \
		for f in ${FASTQ[@]}; do \
				PREFIXS+=($(basename ${f} | sed -e "s/.fastq.*//")); \
			done
		# echo -e "Prefixes:\t${PREFIXS[0]}\t${PREFIXS[1]}"; \
		if [ ${#PREFIXS[0]} -gt ${#PREFIXS[1]} ]; then \
			   long=${PREFIXS[0]}; \
				 short=${PREFIXS[1]}; \
			else \
				long=${PREFIXS[1]}; \
				short=${PREFIXS[0]}; \
			fi
		# echo -e "Short:\t${short}\tLong:\t${long}"; \
			lshort=${#short}
			score=0
			for (( i=0; i<lshort-score; i++ )); do \
			   for (( l=score+1; l<=lshort-i; l++ )); do \
						sub=${short:i:l}; \
						# echo "${sub}";
			      [[ ${long} != *${sub}* ]] && break; \
			      PREFIX=${sub}; \
						score=$l; \
			   done; \
			done

			PREFIX=$(echo ${PREFIX} | sed 's/_*$//')
			# echo -e "Score:\t${score}\tPrefix:\t${PREFIX}"

		if [ ${score} -eq 0 ]; then \
				echo "Paired end files lack similar names: ${FASTQ[@]}"; \
				echo -e ${help_message}; \
				exit 1; \
			fi
	else
		PREFIX=$(basename ${FASTQ} | sed -e "s/.fastq.*//"); \
		PREFIXS=${PREFIX}
	fi

printf "Running:\t\t\t\t%s\n" ${PREFIX}

##Log_parameters
##--------------
declare -A param=( \
		[Date]=${DATE} \
    [Files]=$(echo -e ${FASTQ[@]}) \
    [Genome]=${GENOME} \
    [TrimmedBp]=${FIRST_BASE} \
    [Threads]=${THREADS} \
    [GenomeFasta]=${GENOME_PATH} \
    [ScriptDir]=${SCRIPTDIR} \
		[WorkingDir]=${OUT} \
    [Filter]=${FILTER} \
		[Blacklist]=${DUSTFILE} \
    [G-correct]=${G_CORRECT} \
    [UseVCF]=${USE_VCF}\
    [PathVcf]=${VCF_PATH} \
		[Prefixes]=$(echo -e ${PREFIXS[@]}) \
		[Prefix]=${PREFIX}
	)
printf "Parameters\n" \
	>> ${OUT}/${PREFIX}_parameter.log
printf "\055%.s" {1..10} \
	>> ${OUT}/${PREFIX}_parameter.log
for key in ${!param[@]}; do \
	printf "\n%s\t\t\t%s" ${key} ${param[${key}]} \
		>> ${OUT}/${PREFIX}_parameter.log
done

##Create_tmp_directory
##--------------------
TEMPDIR=$( mktemp -d -p ${OUT} ) #Remove -p flag?
trap "rm -rf ${TEMPDIR}" EXIT
printf "Temporary directory:\t\t\t\t%s\n" ${TEMPDIR}

printf "\n%s\t\t\t%s" "TempDir" ${TEMPDIR} \
  >> ${OUT}/${PREFIX}_parameter.log

##Quality_check_before_filtering
##------------------------------
printf "Quality check before filtering:\t%s\n" ${PREFIX}
for f in ${FASTQ[@]}; do \
		printf "\nQuality check: %s\n" ${f} \
			>> ${OUT}/${PREFIX}_parameter.log; \
		printf "\055%.s" {1..$(( ${#f}+16 ))} \
			>> ${OUT}/${PREFIX}_parameter.log; \
		printf "\n" \
			>> ${OUT}/${PREFIX}_parameter.log; \
		printf "Quality check:\t\t\t\t%s\n" ${f}; \
		fastqc \
			-o ${OUT}/QC/ \
			${f} \
			>> ${OUT}/QC/${PREFIX}_parameter.log 2>&1
	done

##Trim_reads
##----------
printf "Trim and filter reads:\t\t\t%s\n" ${PREFIX}
FASTP_CMD=(
  "fastp"
    "-i ${FASTQ[0]}"
    "-o ${TEMPDIR}/${FASTQ[0]/.fastq.gz/_trimmed.fastq.gz}"
    "--trim_front1 ${FIRST_BASE}"
    "--length_required 30"
    "--qualified_quality_phred 20"
    "--unqualified_percent_limit 50"
    "--thread ${THREADS}"
    "--json ${OUT}/QC/${PREFIX}.fastp.json"
    "--html ${OUT}/QC/${PREFIX}.fastp.html"
  )
if [ ${#FASTQ[@]} -gt 1 ]; then \
    FASTP_CMD+=("-I ${FASTQ[1]}")
    FASTP_CMD+=("-O ${TEMPDIR}/${FASTQ[1]/.fastq.gz/_trimmed.fastq.gz}")
    FASTP_CMD+=("--trim_front2 0")
  fi

eval "${FASTP_CMD[@]}" \
	>> ${OUT}/QC/${PREFIX}_parameter.log 2>&1

printf "Trimming of reads done:\t\t%s\n" ${PREFIX}

for f in ${FASTQ[@]}; do \
    printf "\n%s\t\t\t%s" "TrimmedFile" ${TEMPDIR}/${f/.fastq.gz/_trimmed.fastq.gz} \
      >> ${OUT}/${PREFIX}_parameter.log
  done

##rDNA_filtering
##--------------
if ${FILTER}; then \

		printf "Filter reads:\t\t\t\t%s\n" ${PREFIX}; \
    for f in ${FASTQ[@]}; do \
    		gunzip -c ${TEMPDIR}/${f/.fastq.gz/_trimmed.fastq.gz} \
    				| ${SCRIPTDIR}/bin/rRNAdust \
    					${GENOME_PATH}/${DUSTFILE}/fasta/${DUSTFILE}.fa \
    					-t ${THREADS} \
    					2>> ${OUT}/QC/${PREFIX}_parameter.log \
    				|	gzip > ${TEMPDIR}/${f/.fastq.gz/_filtered.fastq.gz}

    		printf "\n%s\t\t\t%s" "FilteredFile" ${TEMPDIR}/${f/.fastq.gz/_filtered.fastq.gz} \
    			>> ${OUT}/${PREFIX}_parameter.log
      done

	else

		printf "Skip filtering reads.\n"
    for f in ${FASTQ[@]}; do \
  		mv ${TEMPDIR}/${f/.fastq.gz/_trimmed.fastq.gz} ${TEMPDIR}/${f/.fastq.gz/_filtered.fastq.gz}
  		printf "\n%s\t\t\t%s" "UnfilteredFile" ${TEMPDIR}/${f/.fastq.gz/_filtered.fastq.gz} \
  		 >> ${OUT}/${PREFIX}_parameter.log
     done

  fi

##Quality_check_after_filtering
##------------------------------
printf "Quality check after filtering:\t\t%s\n" ${PREFIX}
for f in ${FASTQ[@]}; do
		printf "\n>Quality check: %s\n" ${f} \
			>> ${OUT}/${PREFIX}_parameter.log
		printf "\055%.s" {1..$(( ${#f}+16 ))} \
			>> ${OUT}/${PREFIX}_parameter.log
		printf "\n" \
			>> ${OUT}/${PREFIX}_parameter.log
		printf "Quality check:\t\t\t\t%s\n" ${f}
		fastqc \
			-o ${OUT}/QC/ \
			${TEMPDIR}/${f/.fastq.gz/_filtered.fastq.gz} \
			>> ${OUT}/QC/${PREFIX}_parameter.log 2>&1
	done


##Mapping_and_sorting
##-------------------
printf "Map reads:\t\t\t\t%s\n" ${PREFIX}

FASTQN=() # New vector for adding file extensions
for f in ${FASTQ[@]}; do
    n="${TEMPDIR}/${f/.fastq.gz/_filtered.fastq.gz}";
    FASTQN+=(${n})
  done

STAR_CMD=(
  "STAR"
  "--readFilesType Fastx"
  "--alignEndsType Extend5pOfRead1"
  "--readFilesIn ${FASTQN[@]}"
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
    STAR_CMD+=("--varVCFfile ${VCF_PATH}/${PREFIX}.vcf")
    STAR_CMD+=("--waspOutputMode SAMtag")
    STAR_CMD+=("--outSAMattributes vA vG vW")
  else
    STAR_CMD+=("--outSAMattributes NH HI AS MD nM")
  fi

cmd="${STAR_CMD[@]}"
printf "\n%s\t\t\t%s" "STAR command" ${cmd} \
  >> ${OUT}/${PREFIX}_parameter.log
eval "${STAR_CMD[@]}" \
  >> ${OUT}/${PREFIX}_parameter.log 2>&1

##Remove_temporary_mapping_files
##------------------------------
rm -r ${OUT}/bam_files/${PREFIX}_filtered._STARtmp

printf "\n%s\t\t\t%s" "BamFile" ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
		>> ${OUT}/${PREFIX}_parameter.log

##Index_bam_file
##--------------
samtools index \
	${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam

printf "\n%s\t\t\t%s" "BaiFile" ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.bai \
		>> ${OUT}/${PREFIX}_parameter.log

##Calculate_alignment_complexity
##------------------------------
printf "Calculate alignment complexity:\t%s\n" ${PREFIX}
preseq c_curve \
	-output ${OUT}/QC/${PREFIX}_filtered_preseq.txt \
	-bam ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam

printf "\n%s\t\t\t%s" "PreseqFile" ${OUT}/bam_files/${PREFIX}_filtered_preseq.txt \
		>> ${OUT}/${PREFIX}_parameter.log

##Calculate_alignment_stats
##-------------------------
printf "Calculate alignment stats:\t\t%s\n" ${PREFIX}
samtools stats \
	${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
	> ${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.stats
samtools flagstat \
	${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
	> ${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.fstats

printf "\n%s\t\t\t%s" "BamStats" ${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.stats \
		>> ${OUT}/${PREFIX}_parameter.log
printf "\n%s\t\t\t%s" "FlagStats" ${OUT}/QC/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam.fstats \
		>> ${OUT}/${PREFIX}_parameter.log

##Remove_unmatched_G_addition
##---------------------------
##create_final_bed_file_if_G_CORRECT_is_true
if ${G_CORRECT}; then
  printf "Removing unmatched G on 5' end and creating bedfiles:\t\t%s\n" ${PREFIX}

  #Remove_unmatched_G_addition_on_+_strand
  #---------------------------------------
  samtools view -F 16 \
    ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
      | awk 'BEGIN {FS="\t"} {
        if ($10 ~ /^GG.+/ && $0 ~ /MD:Z:0[A-Z]0[A-Z]/)
          print 2;
        else if ($10 ~ /^G.+/ && $0 ~ /MD:Z:0/)
          print 1;
        else
          print 0
        }' \
    > ${TEMPDIR}/${PREFIX}.plus.counts

  samtools view -b -F 16 \
    ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
    | bamToBed \
    > ${TEMPDIR}/${PREFIX}.plus.bed

  paste ${TEMPDIR}/${PREFIX}.plus.bed ${TEMPDIR}/${PREFIX}.plus.counts \
    > ${TEMPDIR}/${PREFIX}.combined.plus.bed

  cat  ${TEMPDIR}/${PREFIX}.combined.plus.bed \
    | awk 'BEGIN{FS="\t";OFS="\t"} {
        print $1, $2 + $7, $3, $4, $5, $6, $7
      }' \
    > ${TEMPDIR}/${PREFIX}.final.plus.bed

  #Remove_unmatched_G_addition_on_-_strand
  #---------------------------------------
  samtools view -f 16 \
    ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
      | awk 'BEGIN {FS="\t"} {
          if ($10 ~ /.+CC$/ && $0 ~ /MD:Z:.+[A-Z]0[A-Z]0/)
            print 2;
          else if ($10 ~ /.+C$/ && $0 ~ /MD:Z:.+[A-Z]0/)
            print 1;
          else
            print 0
        }' \
    > ${TEMPDIR}/${PREFIX}.minus.counts

  samtools view -b -f 16 \
    ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
      | bamToBed \
      > ${TEMPDIR}/${PREFIX}.minus.bed

  paste ${TEMPDIR}/${PREFIX}.minus.bed ${TEMPDIR}/${PREFIX}.minus.counts \
    > ${TEMPDIR}/${PREFIX}.combined.minus.bed

  cat ${TEMPDIR}/${PREFIX}.combined.minus.bed \
    | awk 'BEGIN{FS="\t";OFS="\t"} {
        print $1, $2 , $3 - $7, $4, $5, $6, $7
      }'  \
    > ${TEMPDIR}/${PREFIX}.final.minus.bed

  #Merge_+/-_and_sort_final_bed_file
  #---------------------------------
  cat ${TEMPDIR}/${PREFIX}.final.plus.bed ${TEMPDIR}/${PREFIX}.final.minus.bed \
    > ${TEMPDIR}/${PREFIX}.final.bed

  cat ${TEMPDIR}/${PREFIX}.final.bed \
    | sort -k 1,1 -k2,2n  \
    > ${OUT}/bed_files/${PREFIX}.bed

  printf "Finished removing unmatched G and creating bed file:\t%s\n" ${PREFIX}

else

  bamToBed \
    -i ${OUT}/bam_files/${PREFIX}_filtered.Aligned.sortedByCoord.out.bam \
    > ${OUT}/bed_files/${PREFIX}.bed
  printf "Skipped removing unmatched G:\t\t%s\n" ${PREFIX}

fi


##Format_CTSS_bed
##---------------
cat ${OUT}/bed_files/${PREFIX}.bed \
  | awk 'BEGIN{OFS="\t"}{
      if ($6 == "+") {
          print $1, $2, $2+1, $6
        } else {
          print $1, $3-1, $3, $6
        }
      }' \
      | sort \
      | uniq -c \
      | awk 'BEGIN{OFS="\t"}{
          print $2, $3, $4, $2 ":" $3 "-" $4 "," $5, $1, $5
        }' \
      | sort -k 1,1 -k 2,2n \
      > ${OUT}/bed_files/${PREFIX}.ctss.bed

cat ${OUT}/bed_files/${PREFIX}.ctss.bed \
  | grep -F ",-" \
  | cut -f 1,2,3,5 \
  | sort -k 1,1 -k 2,2n \
  > ${TEMPDIR}/${PREFIX}.minus.bedgraph
cat ${OUT}/bed_files/${PREFIX}.ctss.bed \
  | grep -F ",+" \
  | cut -f 1,2,3,5 \
  | sort -k 1,1 -k 2,2n \
  > ${TEMPDIR}/${PREFIX}.plus.bedgraph

##Convert_to_bigwig
##-----------------
${SCRIPTDIR}/bin/bedGraphToBigWig \
  ${TEMPDIR}/${PREFIX}.minus.bedgraph \
  ${GENOME_PATH}/${GENOME}/chrom_size/${GENOME}.chrom.sizes \
  ${OUT}/bw_files/${PREFIX}.minus.bw
${SCRIPTDIR}/bin/bedGraphToBigWig \
  ${TEMPDIR}/${PREFIX}.plus.bedgraph \
  ${GENOME_PATH}/${GENOME}/chrom_size/${GENOME}.chrom.sizes \
  ${OUT}/bw_files/${PREFIX}.plus.bw

##Remove_tmp_files
##----------------
rm -rf ${TEMPDIR}

printf "Finished:\t\t\t\t%s\n" ${PREFIX}

exit 0
