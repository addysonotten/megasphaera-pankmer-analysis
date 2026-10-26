import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import warnings
warnings.filterwarnings('ignore')

BASE = "/home/ottena/mags_project/eggnog_analysis"
COG_FILE = f"{BASE}/cog_per_genome.tsv"
LABEL_FILE = "/home/ottena/mags_project/pankmer_analysis/label_lookup_v3.csv"

cog = pd.read_csv(COG_FILE, sep="\t", header=None, names=["genome", "COG"])
labels = pd.read_csv(LABEL_FILE)

def strip_ext(name):
    return re.sub(r"\.(fna|fasta|fa)$", "", str(name), flags=re.IGNORECASE)

labels["original_stripped"] = labels["original_name"].apply(strip_ext)
name_map = dict(zip(labels["original_stripped"], labels["plot_label"]))

def rename_genome(name):
    key = strip_ext(str(name))
    return name_map.get(key, name)

cog["genome"] = cog["genome"].apply(rename_genome)
cog["type"] = cog["genome"].apply(lambda x: "MAG" if x.startswith("MAG") else "Reference")

cog_labels = {
    "A": "RNA processing", "B": "Chromatin structure",
    "C": "Energy metabolism", "D": "Cell cycle control",
    "E": "Amino acid metabolism", "F": "Nucleotide metabolism",
    "G": "Carbohydrate metabolism", "H": "Coenzyme metabolism",
    "I": "Lipid metabolism", "J": "Translation",
    "K": "Transcription", "L": "Replication/repair",
    "M": "Cell wall/membrane", "N": "Cell motility",
    "O": "Post-translational modification", "P": "Inorganic ion transport",
    "Q": "Secondary metabolite biosynthesis", "R": "General function",
    "S": "Unknown function", "T": "Signal transduction",
    "U": "Intracellular trafficking", "V": "Defense mechanisms",
    "W": "Extracellular structures", "X": "Mobilome",
    "Z": "Cytoskeleton"
}

cog_clean = cog[~cog["COG"].isin(["-", ""])].dropna(subset=["COG"])
cog_expanded = cog_clean.copy()
cog_expanded["COG"] = cog_expanded["COG"].astype(str)
cog_expanded = cog_expanded.assign(COG=cog_expanded["COG"].str.split("")).explode("COG")
cog_expanded = cog_expanded[cog_expanded["COG"].str.match(r"[A-Z]")]

# ── GRAPH 1 ────────────────────────────────────────────────────────────────────
print("Making Graph 1...")
cog_counts = cog_expanded.groupby("COG").size().reset_index(name="count")
cog_counts["description"] = cog_counts["COG"].map(cog_labels)
cog_counts = cog_counts.sort_values("count", ascending=False)

fig, ax = plt.subplots(figsize=(18, 9))
bars = ax.bar(range(len(cog_counts)), cog_counts["count"],
              color=sns.color_palette("tab20", len(cog_counts)))

ax.set_xticks(range(len(cog_counts)))
ax.set_xticklabels(
    [f"{row.COG} — {row.description}" for _, row in cog_counts.iterrows()],
    rotation=35, ha="right", fontsize=9
)
ax.set_ylabel("Gene Count", fontsize=12)
ax.set_title("COG Functional Categories Across All Genomes", fontsize=14, fontweight="bold")
ax.set_ylim(0, cog_counts["count"].max() * 1.1)
plt.tight_layout()
plt.savefig(f"{BASE}/COG_all_genomes_v3.pdf", bbox_inches="tight")
plt.savefig(f"{BASE}/COG_all_genomes_v3.png", dpi=150, bbox_inches="tight")
plt.close()
print("Graph 1 saved!")

# ── GRAPH 2 ────────────────────────────────────────────────────────────────────
print("Making Graph 2...")
cog_type = cog_expanded.groupby(["type", "COG"]).size().reset_index(name="count")
cog_type["total"] = cog_type.groupby("type")["count"].transform("sum")
cog_type["proportion"] = cog_type["count"] / cog_type["total"]

cog_order = cog_type.groupby("COG")["proportion"].mean().sort_values(ascending=False).index

fig, ax = plt.subplots(figsize=(18, 9))
colors = {"MAG": "#E87D72", "Reference": "#53B4E5"}
x = range(len(cog_order))
width = 0.35

for i, genome_type in enumerate(["MAG", "Reference"]):
    data = cog_type[cog_type["type"] == genome_type].set_index("COG")
    values = [data.loc[c, "proportion"] if c in data.index else 0 for c in cog_order]
    ax.bar([xi + i * width for xi in x], values, width,
           label=genome_type, color=colors[genome_type], alpha=0.85)

ax.set_xticks([xi + width / 2 for xi in x])
ax.set_xticklabels(
    [f"{c} — {cog_labels.get(c, c)}" for c in cog_order],
    rotation=35, ha="right", fontsize=9
)
ax.set_ylabel("Proportion of Genes", fontsize=12)
ax.set_title("COG Categories: MAG vs Reference Genomes", fontsize=14, fontweight="bold")
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(f"{BASE}/COG_MAG_vs_Reference_v3.pdf", bbox_inches="tight")
plt.savefig(f"{BASE}/COG_MAG_vs_Reference_v3.png", dpi=150, bbox_inches="tight")
plt.close()
print("Graph 2 saved!")

print("\nDone!")
