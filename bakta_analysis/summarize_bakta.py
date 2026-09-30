import pandas as pd
import os
import re

lookup = pd.read_csv("/home/ottena/mags_project/pankmer_analysis/label_lookup_v3.csv")

def strip_ext(name):
    return re.sub(r"\.(fna|fasta|fa)$", "", str(name), flags=re.IGNORECASE)

lookup["key"] = lookup["original_name"].apply(strip_ext)
lookup = lookup.set_index("key")

# Only keep these fields — skips section headers and static software/DB metadata
KEEP_FIELDS = {
    "Length", "Count", "GC", "N50", "N90", "N ratio", "coding density",
    "tRNAs", "tmRNAs", "rRNAs", "ncRNAs", "ncRNA regions", "CRISPR arrays",
    "CDSs", "pseudogenes", "hypotheticals", "sORFs", "gaps",
    "oriCs", "oriVs", "oriTs"
}

bakta_output = "/home/ottena/mags_project/bakta_analysis/output"
rows = []

for genome_dir in sorted(os.listdir(bakta_output)):
    full_dir = os.path.join(bakta_output, genome_dir)
    txt_path = os.path.join(full_dir, f"{genome_dir}.txt")
    if not os.path.isfile(txt_path):
        continue

    label = lookup.loc[genome_dir, "plot_label"] if genome_dir in lookup.index else genome_dir

    stats = {"genome_original": genome_dir, "genome_label": label}
    with open(txt_path) as f:
        for line in f:
            line = line.strip()
            if ":" in line:
                key, val = line.split(":", 1)
                key = key.strip()
                val = val.strip()
                if key in KEEP_FIELDS and val:
                    stats[key] = val

    rows.append(stats)

df = pd.DataFrame(rows)
df.to_csv("/home/ottena/mags_project/bakta_analysis/summary/all_summaries.csv", index=False)
print(f"Wrote summary for {len(df)} genomes")
print(df.columns.tolist())
