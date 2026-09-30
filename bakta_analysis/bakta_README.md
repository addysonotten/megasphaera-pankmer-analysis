# Bakta Annotation Pipeline — Megasphaera MAGs + NCBI Reference Genomes

This documents the workflow used to annotate all 61 genomes (36 MAGs + 25 NCBI
references) with Bakta, a rapid bacterial genome annotation tool.

## 1. Environment setup

This HPC has no `conda`/`mamba`/`anaconda` module and no pre-installed Bakta
module, so Bakta's install used a private Miniconda installation. Bakta
depends on several external compiled tools (tRNAscan-SE, Aragorn, Infernal,
Diamond, HMMER, PILER-CR, AMRFinderPlus) that plain `pip` cannot install —
this is exactly what bioconda packages exist to solve.

```bash
# (Miniconda was already present at ~/miniconda3 on this system)
source ~/miniconda3/etc/profile.d/conda.sh

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# One-time ToS acceptance required by newer conda versions:
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

conda create -n bakta_env -c conda-forge -c bioconda bakta -y
```

Reactivate every new session with:
```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bakta_env
```

## 2. Locating the reference database

Rather than downloading a fresh ~84GB copy, an existing shared lab copy was
found and reused:

```bash
find /data -maxdepth 4 -iname "*bakta*" 2>/dev/null
```

Database path used for all runs:
```
/data/minichj_shared/software/databases/bakta/bakta_db/db
```

(A labmate's existing PBS script at
`/data/minichj_shared/tmcgarry/Scripts/run_bakta_08-16-26.sh` confirmed this
path and the correct command-line flag pattern.)

## 3. PBS batch job

This cluster uses PBS (not SLURM). Bakta annotates one genome per run, so the
job loops over every genome file and writes each genome's output to its own
subfolder.

`run_bakta.sh`:
```bash
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
    [ -e "$GENOME" ] || continue

    NAME=$(basename "$GENOME")
    NAME="${NAME%.*}"

    mkdir -p "$OUTPUT_DIR/$NAME"

    bakta \
        --db "$DB" \
        --output "$OUTPUT_DIR/$NAME" \
        --prefix "$NAME" \
        --force \
        --threads 32 \
        "$GENOME"
done
```

Submit and monitor:
```bash
qsub run_bakta.sh
qstat -u ottena
```

Result: one output subfolder per genome under `output/`, each containing
`.gff3`, `.gbff`, `.tsv`, `.faa`, `.ffn`, `.txt`, `.png`/`.svg`, and `.json`.

## 4. Consolidating protein sequences

Bakta produces two FASTA protein files per genome — the full set and a
subset of "hypothetical" (unknown-function) proteins. Keep these separate,
since a broad `*.faa` glob will grab both (61 genomes × 2 files = 122, not 61):

```bash
mkdir -p all_faa
for dir in output/*/; do
    name=$(basename "$dir")
    cp "$dir/${name}.faa" all_faa/
done
# -> 61 files: full protein sets, input for eggNOG/dbCAN3/VFDB/etc.

mkdir -p all_hypotheticals
for dir in output/*/; do
    name=$(basename "$dir")
    cp "$dir/${name}.hypotheticals.faa" all_hypotheticals/
done
# -> 61 files: hypothetical-protein subsets only
```

## 5. Summarizing annotation stats across all genomes

Each genome's `.txt` file reports basic assembly + annotation stats (CDS
count, tRNAs, rRNAs, CRISPR arrays, coding density, etc.). `summarize_bakta.py`
parses every genome's `.txt` file into one combined CSV, and relabels each row
using the same `label_lookup_v3.csv` built during the PanKmer clustering step
(see `../pankmer_analysis/README.md`), so genome names are human-readable
rather than raw NCBI accessions or contig headers.

```bash
python3 summarize_bakta.py
```
Output: `summary/all_summaries.csv`

## Output folder reference

| Folder | Contents |
|---|---|
| `output/` | Full Bakta results, one subfolder per genome |
| `all_faa/` | Consolidated full protein FASTA, one file per genome |
| `all_hypotheticals/` | Consolidated hypothetical-protein FASTA |
| `summary/` | Combined, relabeled annotation-stats CSV |

## Notes / gotchas encountered
- No conda module on this HPC — used a private Miniconda install instead;
  channels must be `conda-forge` + `bioconda` with strict priority for
  bioconda packages to resolve correctly.
- Newer conda versions require explicitly accepting Anaconda's Terms of
  Service for the `defaults`/`main`/`r` channels before any install will run.
- The PBS scheduler (server `bcm11`) had a temporary outage where `qsub`/
  `qstat` timed out entirely — this was a cluster-side issue, not a script
  problem, and resolved on its own; worth escalating to HPC support if it
  persists.
- `--metric`/database flags aside, always re-verify shared paths for
  permission errors before assuming a folder is usable — an initial
  candidate path (`/data/minichj_shared/databases/bakta`) returned
  "Permission denied," while a second candidate
  (`/data/minichj_shared/software/databases/bakta/bakta_db/db`) worked and
  matched a labmate's already-working script.
- Bakta's `.faa` output includes a `.hypotheticals.faa` variant alongside the
  full `.faa` — a naive `*.faa` glob will double-count files.

## Still to do (per project checklist)
- [ ] Annotate the same genomes with PGAP (NCBI's own pipeline) for
      cross-validation against Bakta
- [ ] AMRFinderPlus (Bakta can run this inline via `--amrfinderplus`, or as a
      standalone step)
- [ ] eggNOG-mapper functional annotation
- [ ] CAZyme analysis (dbCAN3)
- [ ] Virulence factor analysis (VFDB/MetaVF)
- [ ] Defense systems (DefenseFinder, PADLOC)
- [ ] Plasmid/virus detection (geNomad)
- [ ] Biosynthetic gene clusters (antiSMASH, Gecco)
- [ ] DRAM & METABOLIC
- [ ] Secretion systems (MacSyFinder/TXSScan)
- [ ] Mobile genetic elements (MobileOG-db, MOB-suite)

## ⚠️ Sensitive data note
`summary/all_summaries.csv` and any `renamed/` outputs include `genome_label`
values built from participant village, age, and ID — treat these the same as
`rename_filled.csv` from the PanKmer pipeline: do not commit to a public
repository.
