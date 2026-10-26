# MetaVF Virulence Factor Analysis

## Overview
MetaVF was used to identify virulence factor (VF) genes across all 61 genomes. Virulence factors are genes that help bacteria cause disease by facilitating colonization, immune evasion, toxin production, and other pathogenic mechanisms. MetaVF uses nucleotide BLAST against the pVFGC and eVFGC databases derived from VFDB (Virulence Factor Database).

## Input
- Genome assembly FASTA files (`.fna`) produced by Bakta for each genome
- Symlinked into: `input/` folder
- Original location: `/home/ottena/mags_project/bakta_analysis/output/<genome_name>/`

## Database
- pVFGC database: 17,306 sequences (prokaryotic virulence gene clusters)
- eVFGC database: 62,332 sequences (extended virulence gene clusters)
- Located at: `/data/minichj_shared/software/MetaVF_toolkit/databases/`

## Method
- Run as a single PBS job processing all 61 genomes simultaneously via Snakemake
- Used `-m draft` mode for draft genome assemblies
- Identity threshold: 90% (`-ti 90`)
- Coverage threshold: 80% (`-tc 80`)
- 16 cores, 64gb memory, 8 hour walltime
- Conda environment: `MetaVF_toolkit`
- MetaVF toolkit path: `/data/minichj_shared/software/MetaVF_toolkit/`

## Script
- PBS script: `metavf_run.pbs`

## Output
- Individual output folder per genome in `output/mags_project/` containing:
  - `<genome>.summary` — VF gene profiling per gene with TPM abundance, species specificity, mobile element association, and VF category
  - `<genome>.VF_info.summary` — VF profiling per VF with host species, mobile element association, category, and completeness
  - `<genome>.VF.fasta` — FASTA sequences of identified VF genes
  - `<genome>.blastn.filter.best` — filtered best BLAST hits
  - `intermediate_files/` — raw and filtered BLAST output

## Key Results
- No confident virulence factor hits detected across any of the 61 genomes at 90% identity / 80% coverage thresholds
- One weak raw BLAST hit found in GCF_001304715.1 (rtxB from Aeromonas hydrophila, 76.8 score) but filtered by thresholds
- The MetaVF database is heavily biased toward known pathogens (E. coli, Aeromonas etc.)
- The near-zero VF hits is biologically meaningful — these gut commensals lack classical pathogen virulence factor repertoires
- Follow-up VFDB analysis planned with relaxed thresholds and protein-level search (DIAMOND) to detect distant homologs

## Notes
- VFDB direct analysis planned as follow-up using protein sequences and relaxed identity thresholds (60-70%)
- MetaVF conda environment installed manually due to channel priority conflicts with original yaml file
