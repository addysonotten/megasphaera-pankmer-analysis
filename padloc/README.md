# PADLOC Analysis

Detection of anti-phage defense systems across 61 bacterial genomes using [PADLOC](https://github.com/padlocbio/padloc).

## Genomes
- 25 NCBI reference genomes (cow, pig, chicken, and human gut)
- 36 Malawian MAGs (20 Masinjere village, 16 Limera village)

## Software
- PADLOC v2.0.0 (via Apptainer/Singularity container)
- Database: PADLOC-DB v2.0.0
- Container: quay.io/biocontainers/padloc:2.0.0--hdfd78af_1

## Installation Notes
PADLOC conda installation failed due to R dependency conflicts (r-yaml, r-tidyverse).
Solution: pulled via Apptainer from quay.io biocontainers and bound the database directory at runtime.

## Input
Both .faa and .gff3 files from Bakta annotation (required by PADLOC):
/home/ottena/mags_project/bakta_analysis/output/<genome>/

## Output
Per genome, three output files in output/<genome>/:
- *_padloc.csv — detected defense systems
- *_padloc.gff — GFF format output
- *.domtblout — raw HMMER domain table

## Results Summary
- 2,734 total defense system hits across 61 genomes (~45 per genome)
- Most prevalent: CRISPR-Cas (749), Restriction-Modification (468), Innate-like/Other (1040)
- Masinjere MAGs had the highest total hits (987), followed by NCBI Reference (922) and Limera MAGs (825)

## Usage
qsub run_padloc.pbs
