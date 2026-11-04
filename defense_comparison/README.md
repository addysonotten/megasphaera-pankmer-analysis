# Defense System Comparison: DefenseFinder vs PADLOC

Comparative analysis of anti-phage defense systems detected by DefenseFinder and PADLOC across 61 bacterial genomes.

## Overview
Both tools detect anti-phage defense systems but use different approaches:
- DefenseFinder uses MacSyFinder with co-localization rules (multi-gene systems)
- PADLOC uses HMM profiles and is more sensitive for single-gene systems

## Results Summary
- DefenseFinder: 682 total systems (~11 per genome)
- PADLOC: 2,734 total systems (~45 per genome)
- PADLOC detects ~4x more systems, reflecting its broader HMM database
- Both tools agree on major functional categories: RM systems and CRISPR-Cas dominate

## Figures
- defense_comparison_fig1_bar.png — Normalized systems per genome by group
- defense_comparison_fig2_scatter.png — Per-genome correlation between tools
- defense_comparison_fig3_functional.png — Functional category proportions
- defense_comparison_fig4_venn.png — Overlap of system types detected
- defense_comparison_fig5_group_func.png — Functional breakdown by genome group

## Notebook
defense_comparison.ipynb — full analysis and figure generation
