#!/usr/bin/env bash
#PBS -N cutadapt_fpetrilli
#PBS -q omp
#PBS -l ncpus=32
#PBS -l mem=64gb
#PBS -l walltime=100:00:00

source /appli/bioinfo/cutadapt/4.1/env.sh

# Chemins des données
PATH_DATA="path/data"
OUTPATH="/path/out"

# Adaptateurs Illumina standard
ADAPTER_FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Boucle sur chaque sous-dossier
for dir in "$PATH_DATA"/*/; do
    DIR_PATH="${dir%/}"
    SAMPLE=$(basename "${DIR_PATH}")

    # Chemins des fichiers d'entrée et de sortie
    R1="${DIR_PATH}/${SAMPLE}_R1.fastq.gz"
    R2="${DIR_PATH}/${SAMPLE}_R2.fastq.gz"
    OUTPUT_R1="${OUTPATH}/${SAMPLE}_R1_no_adapter.fastq.gz"
    OUTPUT_R2="${OUTPATH}/${SAMPLE}_R2_no_adapter.fastq.gz"

    # Vérifier si les fichiers de sortie existent déjà
    if [[ -f "$OUTPUT_R1" && -f "$OUTPUT_R2" ]]; then
        echo "Les fichiers $OUTPUT_R1 et $OUTPUT_R2 existent déjà."
        # Vérifier si les fichiers d'entrée existent
        if [[ -f "$R1" && -f "$R2" ]]; then
            cutadapt \
                -a $ADAPTER_FWD \
                -A $ADAPTER_REV \
                -o $OUTPUT_R1 \
                -p $OUTPUT_R2 \
                -m 50 \
                $R1 $R2
        else
            echo "Les fichiers $R1 et/ou $R2 n'existent pas pour l'échantillon $SAMPLE."
        fi
    fi
done

echo "Tous les échantillons ont été traités."
