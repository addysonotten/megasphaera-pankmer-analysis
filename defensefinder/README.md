# DefenseFinder Analysis

Systematic detection of anti-phage defense systems across 61 bacterial genomes using [DefenseFinder](https://github.com/mdmparis/defense-finder).

## Genomes
- 25 NCBI reference genomes (cow, pig, chicken, and human gut)
- 36 Malawian MAGs (20 Masinjere village, 16 Limera village)

## Software
- DefenseFinder v2.0.1
- MacSyFinder v2.1.4
- defense-finder-models v3.1.0
- CasFinder v3.1.1

## Installation Notes
CasFinder 3.1.1 ships with XML model definitions using `vers="2.1"`, which is incompatible with MacSyFinder 2.1.4 (supports `vers="2.0"` only). Fix applied:
```bash
find db/CasFinder -name "*.xml" -exec sed -i 's/vers="2\.1"/vers="2.0"/g' {} \;
```

## Input
Protein FASTA files (`.faa`) from Bakta annotation:
`/home/ottena/mags_project/bakta_analysis/all_faa/<genome>.faa`

## Output
Per genome, three output files in `output/<genome>/`:
- `*_defense_finder_systems.tsv` — detected defense systems
- `*_defense_finder_genes.tsv` — genes associated with defense systems
- `*_defense_finder_hmmer.tsv` — raw HMMER hits

## Usage
```bash
qsub run_defensefinder.pbs
```
