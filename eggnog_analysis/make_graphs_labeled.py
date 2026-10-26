import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ── PATHS ──────────────────────────────────────────────────────────────────────
BASE = "/home/ottena/mags_project/eggnog_analysis"
COG_FILE = f"{BASE}/cog_per_genome.tsv"
KEGG_FILE = f"{BASE}/kegg_per_genome.tsv"
LABEL_FILE = "/home/ottena/mags_project/pankmer_analysis/label_lookup_v3.csv"

# ── LOAD DATA ──────────────────────────────────────────────────────────────────
cog = pd.read_csv(COG_FILE, sep="\t", header=None, names=["genome", "COG"])
kegg = pd.read_csv(KEGG_FILE, sep="\t", header=None, names=["genome", "pathway"])
labels = pd.read_csv(LABEL_FILE)

# ── BUILD NAME MAPPING ─────────────────────────────────────────────────────────
import re

def strip_ext(name):
    return re.sub(r"\.(fna|fasta|fa)$", "", str(name), flags=re.IGNORECASE)

labels["original_stripped"] = labels["original_name"].apply(strip_ext)
name_map = dict(zip(labels["original_stripped"], labels["plot_label"]))

def rename_genome(name):
    key = strip_ext(str(name))
    return name_map.get(key, name)

cog["genome"] = cog["genome"].apply(rename_genome)
kegg["genome"] = kegg["genome"].apply(rename_genome)

# Label genome type
cog["type"] = cog["genome"].apply(lambda x: "MAG" if x.startswith("MAG") else "Reference")
kegg["type"] = kegg["genome"].apply(lambda x: "MAG" if x.startswith("MAG") else "Reference")

# COG descriptions
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

# ── CLEAN COG DATA ─────────────────────────────────────────────────────────────
cog_clean = cog[~cog["COG"].isin(["-", ""])].dropna(subset=["COG"])
cog_expanded = cog_clean.copy()
cog_expanded["COG"] = cog_expanded["COG"].astype(str)
cog_expanded = cog_expanded.assign(
    COG=cog_expanded["COG"].str.split("")
).explode("COG")
cog_expanded = cog_expanded[cog_expanded["COG"].str.match(r"[A-Z]")]

# ── GRAPH 1: COG bar chart — all genomes combined ──────────────────────────────
print("Making Graph 1...")
cog_counts = cog_expanded.groupby("COG").size().reset_index(name="count")
cog_counts["description"] = cog_counts["COG"].map(cog_labels)
cog_counts = cog_counts.sort_values("count", ascending=False)

fig, ax = plt.subplots(figsize=(16, 7))
bars = ax.bar(cog_counts["COG"], cog_counts["count"],
              color=sns.color_palette("tab20", len(cog_counts)))
ax.set_xlabel("COG Category", fontsize=12)
ax.set_ylabel("Gene Count", fontsize=12)
ax.set_title("COG Functional Categories Across All Genomes", fontsize=14, fontweight="bold")

for bar, (_, row) in zip(bars, cog_counts.iterrows()):
    ax.text(bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 50,
            row["description"], ha="center", va="bottom",
            fontsize=7, rotation=45)

ax.set_ylim(0, cog_counts["count"].max() * 1.4)
plt.tight_layout()
plt.savefig(f"{BASE}/COG_all_genomes_labeled.pdf", bbox_inches="tight")
plt.savefig(f"{BASE}/COG_all_genomes_labeled.png", dpi=150, bbox_inches="tight")
plt.close()
print("Graph 1 saved!")

# ── GRAPH 2: COG comparison — MAG vs Reference ────────────────────────────────
print("Making Graph 2...")
cog_type = cog_expanded.groupby(["type", "COG"]).size().reset_index(name="count")
cog_type["total"] = cog_type.groupby("type")["count"].transform("sum")
cog_type["proportion"] = cog_type["count"] / cog_type["total"]
cog_type["description"] = cog_type["COG"].map(cog_labels)

cog_order = cog_type.groupby("COG")["proportion"].mean().sort_values(ascending=False).index

fig, ax = plt.subplots(figsize=(16, 7))
colors = {"MAG": "#E87D72", "Reference": "#53B4E5"}
x = range(len(cog_order))
width = 0.35

for i, genome_type in enumerate(["MAG", "Reference"]):
    data = cog_type[cog_type["type"] == genome_type].set_index("COG")
    values = [data.loc[c, "proportion"] if c in data.index else 0 for c in cog_order]
    ax.bar([xi + i * width for xi in x], values, width,
           label=genome_type, color=colors[genome_type], alpha=0.85)

ax.set_xticks([xi + width / 2 for xi in x])
ax.set_xticklabels(cog_order, fontsize=10)
ax.set_xlabel("COG Category", fontsize=12)
ax.set_ylabel("Proportion of Genes", fontsize=12)
ax.set_title("COG Categories: MAG vs Reference Genomes", fontsize=14, fontweight="bold")
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig(f"{BASE}/COG_MAG_vs_Reference_labeled.pdf", bbox_inches="tight")
plt.savefig(f"{BASE}/COG_MAG_vs_Reference_labeled.png", dpi=150, bbox_inches="tight")
plt.close()
print("Graph 2 saved!")

# ── GRAPH 3: KEGG pathway heatmap ─────────────────────────────────────────────
print("Making Graph 3...")
kegg_clean = kegg.copy()
kegg_clean["pathway"] = kegg_clean["pathway"].astype(str)
kegg_clean["pathway"] = kegg_clean["pathway"].str.split(",")
kegg_clean = kegg_clean.explode("pathway")
kegg_clean = kegg_clean[kegg_clean["pathway"].str.startswith("map")]

top_pathways = (kegg_clean.groupby("pathway").size()
                .sort_values(ascending=False)
                .head(40).index.tolist())

kegg_filtered = kegg_clean[kegg_clean["pathway"].isin(top_pathways)]
kegg_matrix = kegg_filtered.drop_duplicates(["genome", "pathway"])
kegg_matrix = kegg_matrix.copy()
kegg_matrix["present"] = 1
kegg_pivot = kegg_matrix.pivot(index="genome", columns="pathway", values="present").fillna(0)

kegg_pivot["type"] = kegg_pivot.index.map(lambda x: "MAG" if x.startswith("MAG") else "Reference")
kegg_pivot = kegg_pivot.sort_values("type")
kegg_pivot = kegg_pivot.drop(columns="type")

row_colors = pd.Series(
    ["#E87D72" if g.startswith("MAG") else "#53B4E5" for g in kegg_pivot.index],
    index=kegg_pivot.index
)

fig, ax = plt.subplots(figsize=(20, 22))
sns.heatmap(kegg_pivot, cmap=["#f0f0f0", "#2166ac"],
            linewidths=0.3, linecolor="white",
            cbar_kws={"label": "Present/Absent", "shrink": 0.3},
            ax=ax, vmin=0, vmax=1)

for i, (genome, color) in enumerate(zip(kegg_pivot.index, row_colors)):
    ax.add_patch(mpatches.Rectangle((-0.8, i), 0.5, 1,
                                     color=color, clip_on=False))

ax.set_title("Top 40 KEGG Pathways Across All Genomes", fontsize=14, fontweight="bold", pad=20)
ax.set_xlabel("KEGG Pathway", fontsize=11)
ax.set_ylabel("Genome", fontsize=11)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)

mag_patch = mpatches.Patch(color="#E87D72", label="MAG")
ref_patch = mpatches.Patch(color="#53B4E5", label="Reference")
ax.legend(handles=[mag_patch, ref_patch], loc="upper right",
          bbox_to_anchor=(1.12, 1.05), fontsize=10, title="Genome Type")

plt.tight_layout()
plt.savefig(f"{BASE}/KEGG_heatmap_labeled.pdf", bbox_inches="tight")
plt.savefig(f"{BASE}/KEGG_heatmap_labeled.png", dpi=150, bbox_inches="tight")
plt.close()
print("Graph 3 saved!")

print("\nAll done! Files saved to:", BASE)
