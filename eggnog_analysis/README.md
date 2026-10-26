# EggNOG-mapper Analysis

## Overview
EggNOG-mapper (v2.1.15) was used to perform functional annotation of all 61 genomes. EggNOG-mapper searches protein sequences against the EggNOG database of orthologous groups to assign COG functional categories, KEGG pathways, GO terms, and gene names/descriptions to each protein.

## Input
- Protein FASTA files (`.faa`) produced by Bakta for each genome
- Located at: `/home/ottena/mags_project/bakta_analysis/output/<genome_name>/`

## Database
- EggNOG database version: 5.0.2
- Located at: `/data/minichj_shared/databases/eggnog/v5/`

## Method
- Run as a PBS array job (1 job per genome, 61 total)
- Used `-m diamond` search mode with `--dmnd_algo ctg`
- 8 threads per job, 32gb memory, 4 hour walltime
- Conda environment: `eggnog_env`

## Script
- PBS array script: `eggnog_array.pbs`

## Output
- Individual `.emapper.annotations` files per genome in `output/` folder
- Each row represents one protein with columns including: query, seed ortholog, evalue, score, eggNOG OGs, max annotation level, COG category, description, preferred name, GO terms, EC number, KEGG ko, KEGG pathway, KEGG module, KEGG reaction, KEGG rclass, BRITE, KEGG TC, CAZy, BiGG reaction, PFAMs
- Summary files for downstream analysis:
  - `cog_per_genome.tsv` — COG category per protein per genome
  - `kegg_per_genome.tsv` — KEGG pathway per protein per genome

## Key Results
- 136,140 total protein annotations across 61 genomes (~2,200 proteins per genome on average)
- 46,690 KEGG pathway annotations
- All genomes dominated by metabolism and genetic information processing categories
- MAG genomes show similar functional profiles to each other regardless of village
- NCBI reference genomes have more total annotations reflecting complete genome assemblies vs MAG fragments

## Figures
- `COG_all_genomes_v3.pdf/png` — COG functional categories across all genomes
- `COG_MAG_vs_Reference_v3.pdf/png` — COG category comparison between MAG and reference genomes
- `KEGG_stacked_bar.pdf/png` — KEGG pathway categories per genome stacked bar chart
