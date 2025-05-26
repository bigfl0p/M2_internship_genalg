#!/usr/bin/env bash
#PBS -N matrice_count
#PBS -q omp
#PBS -l ncpus=64
#PBS -l mem=32gb
#PBS -l walltime=20:00:00


source /appli/bioinfo/htseq-count/0.13.5/env.sh

in="/path/in"
path="/path/out"
gff="/path/genomes.gff"
bam="/path/aligned.bam"

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

    
