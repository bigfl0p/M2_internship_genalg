#!/usr/bin/env bash
#PBS -N remove_rRNA
#PBS -q omp
#PBS -l mem=32G
#PBS -l ncpus=32
#PBS -l walltime=20:00:00
### PBS -W depend=afterany:

# Load bowtie2 environment
source /appli/bioinfo/bowtie2/2.5.4/env.sh

cd ${PBS_O_WORKDIR}

# Define the path to the data
path="/home/datawork-lpba/Blueremediomics_Lou/Metatranscriptomics/data_analysis"
d_in="1_trimmed_reads"
d_out="2_rrna_removed"
db="/home1/datahome/fpetrill/bt2index_rRNA/rRNA_of_references"

SAMPLE=B3_bact_J14
R1=${path}/${d_in}/${SAMPLE}/${SAMPLE}_R1_no_adapter.fastq.gz
R2=${path}/${d_in}/${SAMPLE}/${SAMPLE}_R2_no_adapter.fastq.gz
f_out="${path}/${d_out}/${SAMPLE}/${SAMPLE}"

bowtie2 -q -p ${NCPUS} \
        -1 ${R1} \
        -2 ${R2} \
        -x ${db} \
        --sensitive-local \
        --al-conc-gz ${f_out}_rRNA.R%.fastq.gz \
        --un-conc-gz ${f_out}_otherRNA.R%.fastq.gz >& ${path}/report/${PBS_JOBNAME}.e${PBS_JOBID%.*} 2>&1
        
        #--un-conc-gz ${f_out}_otherRNA.R%.fastq.gz



