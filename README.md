# PanKmer Pipeline — Megasphaera MAGs vs. NCBI Reference Genomes

This documents the full workflow used to build a k-mer-based clustermap comparing
36 in-house MAGs against 25+ NCBI reference genomes, to test whether genome
similarity reflects host species, geographic origin, or neither.

## 1. Environment setup

PanKmer requires Python 3.8–3.12; Python 3.14 causes a segfault (incompatible
compiled extensions). Use a specific, older module version.

```bash
module load python/3.12.8
cd ~/mags_project/pankmer_analysis
python3 -m venv pankmer_env
source pankmer_env/bin/activate
pip install --upgrade pip setuptools wheel
pip install pankmer
```

Reactivate every new session with:
```bash
module load python/3.12.8
source ~/mags_project/pankmer_analysis/pankmer_env/bin/activate
```

## 2. Splitting the original multi-FASTA into individual MAG files

```bash
cd ~/mags_project/pankmer_analysis/input_fastas
mkdir -p split_output
awk '/^>/{
  if (out) close(out)
  header = substr($0, 2)
  gsub(/[ \/:,]/, "_", header)
  out = "split_output/" header ".fasta"
}
{ print > out }' "Megasphaera_MAGs (1).fasta"
```

Result: 36 individual `.fasta` files, one per MAG.

## 3. Downloading NCBI reference genomes

Downloaded as a batch zip (`ncbi_dataset.zip`) from the NCBI Genome browser
(selected genomes → Download → Genome sequences FASTA), uploaded to
`input_fastas/`, then unzipped:

```bash
cd ~/mags_project/pankmer_analysis/input_fastas
unzip ncbi_dataset.zip -d ncbi_genomes
```

Actual `.fna` files are nested under:
`ncbi_genomes/ncbi_dataset/data/<accession>/<accession>_genomic.fna`

## 4. Combining all genomes into one folder

```bash
mkdir -p ~/mags_project/pankmer_analysis/all_genomes
cp ~/mags_project/pankmer_analysis/input_fastas/split_output/*.fasta \
   ~/mags_project/pankmer_analysis/all_genomes/
find ~/mags_project/pankmer_analysis/input_fastas/ncbi_genomes/ -name "*.fna" \
   -exec cp {} ~/mags_project/pankmer_analysis/all_genomes/ \;
```

**Duplicate removal:** NCBI sometimes returns both a GenBank (`GCA_`) and
RefSeq (`GCF_`) accession for the same underlying assembly. Kept RefSeq,
removed GenBank duplicates:

```bash
cd ~/mags_project/pankmer_analysis/all_genomes
for f in GCA_*; do
  suffix=$(echo "$f" | sed -E 's/^GCA_[0-9]+\.[0-9]+_//')
  match=$(ls GCF_*_"$suffix" 2>/dev/null)
  if [ -n "$match" ]; then
    rm "$f"
  fi
done
```

Final count: 61 genome files (36 MAGs + 25 NCBI, after de-duplication).

## 5. Building the k-mer index

```bash
pankmer index -g ~/mags_project/pankmer_analysis/all_genomes/ \
              -o ~/mags_project/pankmer_analysis/pankmer_output_v2/
```

## 6. Generating the adjacency matrix

```bash
pankmer adj-matrix -i ~/mags_project/pankmer_analysis/pankmer_output_v2/ \
                    -o ~/mags_project/pankmer_analysis/adj_matrix.csv
```
`adj-matrix` has no `--metric` option — it just computes the raw comparison
data. The similarity metric is chosen later, at plotting time.

## 7. Relabeling the matrix with readable sample names

Raw filenames (NCBI accessions, contig headers) are unreadable on a 61-genome
plot, so genomes were relabeled based on collected metadata
(`rename_filled.csv`: host/source for NCBI genomes; village/age/ID for MAGs).

`relabel_matrix.py`:
```python
import pandas as pd
import re

df = pd.read_csv("adj_matrix.csv", index_col=0)
meta = pd.read_csv("rename_filled.csv")

def strip_ext(name):
    return re.sub(r"\.(fna|fasta|fa)$", "", str(name), flags=re.IGNORECASE)

def strip_word(val, word):
    return re.sub(rf"\b{word}\b", "", str(val), flags=re.IGNORECASE).strip()

def clean(val):
    val = str(val).strip()
    return re.sub(r"[^A-Za-z0-9]+", "", val) if val and val.lower() != "nan" else ""

def clean_number(val):
    val = str(val).strip()
    if not val or val.lower() == "nan":
        return ""
    return re.sub(r"[^0-9.]", "", val)  # keeps literal decimal point

meta["filename_noext"] = meta["filename"].apply(strip_ext)
meta = meta.set_index("filename_noext")

def build_label(fname):
    key = strip_ext(fname)
    if key not in meta.index:
        return fname[:20]
    row = meta.loc[key]
    gtype = row.get("type", "")
    if gtype == "NCBI":
        host = clean(strip_word(row.get("host", ""), "NCBI"))
        source = clean(strip_word(row.get("source", ""), "NCBI"))
        parts = [p for p in [host, source] if p]
        label = "NCBI_" + "_".join(parts) if parts else "NCBI_" + fname[:15]
    else:
        village = clean(strip_word(row.get("village", ""), "MAG"))
        age = clean_number(row.get("age", ""))
        pid = clean(row.get("id", ""))
        parts = [village]
        if pid:
            parts.append(f"ID{pid}")
        if age:
            parts.append(f"Age{age}")
        parts = [p for p in parts if p]
        label = "MAG_" + "_".join(parts) if parts else "MAG_" + fname[:15]
    return label

raw_labels = [build_label(f) for f in df.index]

seen = {}
final_labels = []
for lbl in raw_labels:
    if lbl in seen:
        seen[lbl] += 1
        final_labels.append(f"{lbl}_{seen[lbl]}")
    else:
        seen[lbl] = 0
        final_labels.append(lbl)

pd.DataFrame({
    "plot_label": final_labels,
    "original_name": df.index
}).to_csv("label_lookup_v3.csv", index=False)

matched = sum(1 for f in df.index if strip_ext(f) in meta.index)
print(f"Matched {matched} of {len(df.index)} genomes to metadata")

df.index = final_labels
df.columns = final_labels
df.to_csv("adj_matrix_relabeled.csv")
```

Run with:
```bash
python3 relabel_matrix.py
```

## 8. Building the clustermap

```bash
pankmer clustermap -i ~/mags_project/pankmer_analysis/adj_matrix_relabeled.csv \
                    -o ~/mags_project/pankmer_analysis/clustermap_v3.pdf \
                    --width 24 --height 24
```

Default `--metric` is **ANI** (Average Nucleotide Identity) — this run did not
override it, so ANI is the similarity metric shown in `clustermap_v3.pdf`.

To generate a Jaccard-based version for comparison:
```bash
pankmer clustermap -i ~/mags_project/pankmer_analysis/adj_matrix_relabeled.csv \
                    -o ~/mags_project/pankmer_analysis/clustermap_jaccard.pdf \
                    --metric jaccard --width 24 --height 24
```

## 9. Interactive metadata entry tool

`interactive_rename.py` — prompts one-by-one for host/source (NCBI genomes) or
village/age/ID (MAGs), auto-detecting type by filename pattern, and saves
progress after every entry so it can be safely stopped/resumed:

```bash
python3 interactive_rename.py
```

## Key finding

Clustering (ANI, complete linkage) shows MAGs group primarily by **village of
origin** (Masinjere vs. Limera) rather than by host species, and both village
clusters are distinct from all NCBI reference genomes, including other
human-associated *Megasphaera*. This suggests population-level/geographic
structure as the dominant signal — worth confirming isn't a batch/technical
artifact before drawing biological conclusions.

## Notes / gotchas encountered
- `python3 -m venv` bakes in whatever Python version was active at creation
  time — reactivating an existing venv does NOT change its underlying Python
  version. Must delete and recreate the venv after loading a different module.
- NCBI zip downloads can include both GCA and GCF accessions for the same
  assembly — check for and remove duplicates before indexing.
- `rename_filled.csv` contains participant-level metadata (village, age, ID)
  — treat as sensitive; do not commit to a public repository.
