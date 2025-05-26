#!/usr/bin/env bash
#PBS -N star_mapping
#PBS -q omp
#PBS -l mem=64G
#PBS -l ncpus=32
#PBS -l walltime=20:00:00
#PBS -W depend=afterany:2919084

cd ${PBS_O_WORKDIR}

# define path, in and out
path="/home/datawork-lpba/Blueremediomics_Lou/Metatranscriptomics/data_analysis"
idx="/home/datawork-lpba/Blueremediomics_Lou/Metatranscriptomics/Génomes/Prymnesium_parvum/PrymneIndexStar"


SAMPLE=PS4_bact_J14
R1=${path}/2_rrna_removed/${SAMPLE}/${SAMPLE}_otherRNA.R1.fastq.gz
R2=${path}/2_rrna_removed/${SAMPLE}/${SAMPLE}_otherRNA.R2.fastq.gz
out=${path}/3_mapping_pp/${SAMPLE}/ARNm

## Remove empty reads ##
#zcat ${R1} | awk 'BEGIN {RS="@"; ORS=""; FS="\n"; OFS="\n";} {if ($2 != "" && $4 != "") print "@"$0}' | gzip > ${R1%.R1.fastq.gz}_clean.R1.fastq.gz
#zcat ${R2} | awk 'BEGIN {RS="@"; ORS=""; FS="\n"; OFS="\n";} {if ($2 != "" && $4 != "") print "@"$0}' | gzip > ${R2%.R2.fastq.gz}_clean.R2.fastq.gz

## Sorted sequences ##
source /appli/bioinfo/seqkit/2.9.0/env.sh
if [[ -f "${R1%.fastq.gz}.paired.fastq.gz" && -f "${R2%.fastq.gz}.paired.fastq.gz" ]];
then

     echo "Les lectures sont déjà triées"

else

     seqkit pair -1 ${R1} -2 ${R2} -o ${path}/2_rrna_removed/${SAMPLE} >& ${path}/report/${PBS_JOBNAME}.e${PBS_JOBID%.*}_seqkit 2>&1

fi


source /appli/bioinfo/seqkit/2.9.0/delenv.sh

## Mapping ##
source /appli/bioinfo/star/2.7.10b/env.sh

     STAR --runThreadN ${NCPUS} \
          --genomeDir ${idx} \
          --readFilesIn ${R1%.fastq.gz}.paired.fastq.gz ${R2%.fastq.gz}.paired.fastq.gz \
          --readFilesCommand gunzip -c \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMmapqUnique 10 \
          --alignIntronMin 10 \
          --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
          --outFileNamePrefix ${out}/${SAMPLE}_mRNA. >& ${path}/report/${PBS_JOBNAME}.e${PBS_JOBID%.*}_star 2>&1

source /appli/bioinfo/star/2.7.10b/delenv.sh

