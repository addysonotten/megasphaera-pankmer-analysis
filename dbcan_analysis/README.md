# dbCAN CAZyme Analysis

## Overview
dbCAN3 (v5.2.9) was used to identify Carbohydrate-Active enZymes (CAZymes) across all 61 genomes. CAZymes are enzymes involved in the synthesis, degradation, and modification of carbohydrates. dbCAN uses three complementary methods (DIAMOND, HMM, and dbCAN-sub HMM) to annotate CAZyme families, with hits called by at least 2 out of 3 methods considered confident.

## Input
- Protein FASTA files (`.faa`) produced by Bakta for each genome
- Located at: `/home/ottena/mags_project/bakta_analysis/output/<genome_name>/`

## Database
- dbCAN database version: downloaded August 2026
- Located at: `/data/minichj_shared/software/databases/dbcan/`
- Contains: `CAZy.dmnd`, `dbCAN.hmm`, `dbCAN-sub.hmm`, `TF.hmm`, `STP.hmm`, `TCDB.dmnd`, `PUL.dmnd`, `sulfatlas_db.dmnd`, `peptidase_db.dmnd`

## Method
- Run as a PBS array job (1 job per genome, 61 total)
- Used `--mode protein` since input is protein sequences from Bakta
- All three methods used: `--methods diamond,hmm,dbCANsub`
- 8 threads per job, 32gb memory, 4 hour walltime
- Conda environment: `dbcan-env`

## Script
- PBS array script: `dbcan_array.pbs`

## Output
- Individual output folder per genome in `output/` containing:
  - `overview.tsv` — main results combining all three methods, confident hits (≥2 tools agree)
  - `dbCAN_hmm_results.tsv` — raw HMM search results
  - `dbCANsub_hmm_results.tsv` — raw dbCAN-sub HMM results
  - `dbCANsub_hmm_raw.tsv` — unfiltered dbCAN-sub results
  - `diamond.out` — raw DIAMOND search results
  - `uniInput.faa` — reformatted input protein sequences

## CAZyme Classes
- **GH** — Glycoside Hydrolases: break down carbohydrates
- **GT** — Glycosyltransferases: build carbohydrates and cell wall components
- **PL** — Polysaccharide Lyases: degrade uronic acid-containing polysaccharides
- **CE** — Carbohydrate Esterases: remove ester modifications from carbohydrates
- **AA** — Auxiliary Activities: oxidoreductases that modify carbohydrates
- **CBM** — Carbohydrate Binding Modules: non-catalytic domains that target carbohydrate substrates

## Key Results
- All genomes dominated by GH and GT classes
- Cow/rumen genomes have the highest CAZyme counts consistent with specialized plant cell wall degradation
- MAG genomes from both Masinjere and Limera show simple CAZyme profiles (primarily GH and GT only)
- No significant AA or PL enzymes detected in MAG genomes
- One Limera MAG clusters with human gut reference genomes based on CAZyme profile

## Figures
- `CAZyme_stacked_per_genome.pdf/png` — CAZyme class distribution per genome
- `CAZyme_grouped_by_host.pdf/png` — CAZyme classes grouped by host/source type
- `CAZyme_family_heatmap.pdf/png` — top 40 CAZyme families presence/absence heatmap
