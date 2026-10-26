# AMRFinderPlus Analysis

## Overview
AMRFinderPlus (v4.2.7) was used to identify antimicrobial resistance (AMR) genes across all 61 genomes. AMRFinderPlus searches protein sequences against the NCBI Antimicrobial Resistance Gene Database to identify genes conferring resistance to antibiotics, biocides, and heavy metals.

## Input
- Protein FASTA files (`.faa`) and GFF3 annotation files (`.gff3`) produced by Bakta for each genome
- Located at: `/home/ottena/mags_project/bakta_analysis/output/<genome_name>/`

## Database
- AMRFinderPlus database version: 2026-05-15.1
- Located at: `/data/minichj_shared/software/databases/bakta/bakta_db/db/amrfinderplus-db/2026-05-15.1`

## Method
- Run as a PBS array job (1 job per genome, 61 total)
- Used `--annotation_format bakta` to integrate protein and genomic coordinate information
- 8 threads per job, 16gb memory, 2 hour walltime

## Script
- PBS array script: `amrfinder_array.pbs`

## Output
- Individual TSV files per genome in `output/` folder
- Combined master table: `amrfinder_all_genomes.tsv`
- Each row represents one AMR gene hit with columns: genome name, gene symbol, sequence name, genomic coordinates, strand, element type, class, subclass, method, target length, reference length, coverage, percent identity, alignment length, accession, gene name, product name, scope

## Key Results
- 71 total AMR gene hits across 61 genomes
- 27 out of 61 genomes had at least one hit
- Dominant resistance classes: Tetracycline (33 hits) and Beta-lactam (20 hits)
- All Masinjere MAGs consistently carry tet(M) + blaACI-1
- Limera MAGs had zero AMR hits
- Pig genomes carried the most diverse AMR profiles including optrA (last-resort antibiotic resistance)
