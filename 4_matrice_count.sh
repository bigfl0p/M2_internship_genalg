#!/usr/bin/env bash
#PBS -N matrice_count
#PBS -q omp
#PBS -l ncpus=64
#PBS -l mem=32gb
#PBS -l walltime=20:00:00


source /appli/bioinfo/htseq-count/0.13.5/env.sh

in="/home/datawork-lpba/Blueremediomics_Lou/Metatranscriptomics/data_analysis/3_mapping_pp"
path="/home/datawork-lpba/Blueremediomics_Lou/Metatranscriptomics/data_analysis/4_counts"
gff="/home/datawork-lpba/Blueremediomics_Lou/Metatranscriptomics/Génomes/SynCom/species_annotations/all_genomes.gff"


bam="/home/datawork-lpba/Blueremediomics_Lou/Metatranscriptomics/data_analysis/3_mapping_bact/BS2_bact_J14/SC/ARNm/BS2_bactAligned.sortedByCoord.out.bam"

htseq-count -m union \
        --nonunique random \
        -f bam \
        -r name \
        -a 2 \
        -s reverse \
        -t CDS \
        -i ID  \
        $bam \
        $gff > ${path}/test3.csv

    
