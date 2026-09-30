#!/bin/bash
#PBS -N bakta_all_genomes
#PBS -l nodes=1:ppn=32
#PBS -l mem=64gb
#PBS -o /home/ottena/mags_project/bakta_analysis/bakta_all_genomes.out
#PBS -e /home/ottena/mags_project/bakta_analysis/bakta_all_genomes.err

cd "$PBS_O_WORKDIR"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bakta_env

INPUT_DIR="/home/ottena/mags_project/pankmer_analysis/all_genomes"
OUTPUT_DIR="/home/ottena/mags_project/bakta_analysis/output"
DB="/data/minichj_shared/software/databases/bakta/bakta_db/db"

mkdir -p "$OUTPUT_DIR"

for GENOME in "$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fna; do
    [ -e "$GENOME" ] || continue   # skip if one of the two patterns matches nothing

    NAME=$(basename "$GENOME")
    NAME="${NAME%.*}"   # strip extension

    mkdir -p "$OUTPUT_DIR/$NAME"

    bakta \
        --db "$DB" \
        --output "$OUTPUT_DIR/$NAME" \
        --prefix "$NAME" \
        --force \
        --threads 32 \
        "$GENOME"
done
