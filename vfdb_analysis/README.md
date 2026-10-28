# VFDB Virulence Factor Analysis

## Overview
Virulence factor (VF) genes were identified across all 61 genomes using DIAMOND protein-level search against the VFDB (Virulence Factor Database) Set B database. Two searches were run with different stringency thresholds to distinguish true functional homologs from distant evolutionary relatives.

## Input
- Protein FASTA files (`.faa`) produced by Bakta for each genome
- Located at: `/home/ottena/mags_project/bakta_analysis/output/<genome_name>/`

## Database
- VFDB Set B (full dataset including experimentally verified and predicted homologs)
- Downloaded: October 2026 from http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
- Size: ~5.6 MB compressed
- Located at: `/home/ottena/mags_project/vfdb_analysis/database/VFDB_setB_pro.dmnd`
- DIAMOND database built with: `diamond makedb`

## Methods

### Run 1 — VFDB Relaxed
- Tool: DIAMOND blastp v2.2.8
- Identity threshold: 30%
- Query coverage threshold: 50%
- E-value: 1e-5
- Mode: sensitive
- Max target sequences: 1 (best hit only)
- Threads: 8 per job, 16gb memory, 2 hour walltime
- PBS array job: 1 job per genome, 61 total
- Conda environment: bakta_env
- Script: `vfdb_array.pbs`
- Result: 248-331 hits per genome, all 61 genomes had hits
- Interpretation: At 30% identity most hits represent very distant evolutionary homologs — genes sharing ancient ancestry with VFDB entries but functioning as normal metabolic/housekeeping genes in these gut commensals. Dominated by nutritional/metabolic factor homologs.

### Run 2 — VFDB Strict (PRIMARY ANALYSIS)
- Tool: DIAMOND blastp v2.2.8
- Identity threshold: 50%
- Query coverage threshold: 70%
- E-value: 1e-10
- Mode: sensitive
- Max target sequences: 1 (best hit only)
- Threads: 8 per job, 16gb memory, 2 hour walltime
- PBS array job: 1 job per genome, 61 total
- Conda environment: bakta_env
- Script: `vfdb_strict_array.pbs`
- Result: 23-44 hits per genome, all 61 genomes had hits
- Interpretation: At 50% identity these represent genuine functional homologs of known virulence factors. Categories dominated by immune modulation (polysaccharide capsule, LPS, LOS, PSA), nutritional/metabolic factors (heme biosynthesis, GuaA/B), and adherence (EF-Tu).

## Output Files
- `output/` — relaxed threshold results, one TSV per genome (`<genome>_vfdb.tsv`)
- `output_strict/` — strict threshold results, one TSV per genome (`<genome>_vfdb_strict.tsv`)
- Each TSV contains columns: query protein ID, VFDB subject ID, percent identity, alignment length, mismatches, gap opens, query start/end, subject start/end, e-value, bitscore, subject title (includes VF name, VF ID, VF category, and source organism)

## Key Results
- Relaxed (30% identity): 248-331 hits per genome — mostly distant metabolic homologs
- Strict (50% identity): 23-44 hits per genome — genuine functional homologs
- Dominant VF categories in strict results: immune modulation (polysaccharide capsule, LPS, LOS, PSA), nutritional/metabolic factors (heme biosynthesis, GuaA/B), adherence (EF-Tu), effector delivery (T6SS-II)
- Immune modulation dominance is biologically meaningful — gut commensals use surface polysaccharides to persist in the host without triggering excessive immune responses
- No classical toxin or invasion genes detected — consistent with commensal rather than pathogenic lifestyle

## Important Note
Hits in gut commensals against VFDB represent functional homologs of virulence-associated genes, not evidence of pathogenicity. The same gene families used by pathogens for virulence (capsule biosynthesis, iron acquisition) are used by gut commensals for normal host interaction and colonization. These results reflect host interaction capability rather than pathogenic potential.

## Software
- DIAMOND v2.2.8 (available in bakta_env)
