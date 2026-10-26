import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import re
import os
import warnings
warnings.filterwarnings('ignore')

# ── PATHS ──────────────────────────────────────────────────────────────────────
BASE = "/home/ottena/mags_project/dbcan_analysis"
DBCAN_OUT = f"{BASE}/output"
LABEL_FILE = "/home/ottena/mags_project/pankmer_analysis/label_lookup_v3.csv"

# ── LOAD LABELS ────────────────────────────────────────────────────────────────
labels = pd.read_csv(LABEL_FILE)

def strip_ext(name):
    return re.sub(r"\.(fna|fasta|fa)$", "", str(name), flags=re.IGNORECASE)

labels["original_stripped"] = labels["original_name"].apply(strip_ext)
name_map = dict(zip(labels["original_stripped"], labels["plot_label"]))

def rename_genome(name):
    key = strip_ext(str(name))
    return name_map.get(key, name)

# ── LOAD AND COMBINE ALL OVERVIEW FILES ────────────────────────────────────────
print("Loading all overview.tsv files...")
all_data = []
for genome_dir in os.listdir(DBCAN_OUT):
    overview_path = os.path.join(DBCAN_OUT, genome_dir, "overview.tsv")
    if os.path.exists(overview_path):
        df = pd.read_csv(overview_path, sep="\t")
        df["genome"] = genome_dir
        all_data.append(df)

combined = pd.concat(all_data, ignore_index=True)
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import re
import os
import warnings
warnings.filterwarnings('ignore')

# ── PATHS ──────────────────────────────────────────────────────────────────────
BASE = "/home/ottena/mags_project/dbcan_analysis"
DBCAN_OUT = f"{BASE}/output"
LABEL_FILE = "/home/ottena/mags_project/pankmer_analysis/label_lookup_v3.csv"

# ── LOAD LABELS ────────────────────────────────────────────────────────────────
labels = pd.read_csv(LABEL_FILE)

def strip_ext(name):
    return re.sub(r"\.(fna|fasta|fa)$", "", str(name), flags=re.IGNORECASE)

labels["original_stripped"] = labels["original_name"].apply(strip_ext)
name_map = dict(zip(labels["original_stripped"], labels["plot_label"]))

def rename_genome(name):
    key = strip_ext(str(name))
    return name_map.get(key, name)

# ── LOAD AND COMBINE ALL OVERVIEW FILES ────────────────────────────────────────
print("Loading all overview.tsv files...")
all_data = []
for genome_dir in os.listdir(DBCAN_OUT):
    overview_path = os.path.join(DBCAN_OUT, genome_dir, "overview.tsv")
    if os.path.exists(overview_path):
        df = pd.read_csv(overview_path, sep="\t")
        df["genome"] = genome_dir
        all_data.append(df)

combined = pd.concat(all_data, ignore_index=True)

# Keep only confident calls (>=2 tools agree)
confident = combined[combined["Recommend Results"] != "-"].copy()
confident = confident[confident["Recommend Results"].notna()]

# Extract CAZyme family (e.g. GH13_1 -> GH13, GT2 -> GT2)
def get_family(result):
    result = str(result)
    # Remove subfamilies and domain coordinates like (1-200)
    result = re.sub(r"\(.*?\)", "", result)
    result = re.sub(r"_\d+$", "", result.strip())
    return result.strip()

def get_class(family):
    for prefix in ["GH", "GT", "PL", "CE", "AA", "CBM"]:
        if family.startswith(prefix):
            return prefix
    return "Other"

confident["family"] = confident["Recommend Results"].apply(get_family)
confident["class"] = confident["family"].apply(get_class)

# Apply genome labels
confident["label"] = confident["genome"].apply(rename_genome)

# Assign host type
def get_host(label):
    label_lower = label.lower()
    if "pig" in label_lower or "swine" in label_lower:
        return "Pig"
    elif "cow" in label_lower or "rumen" in label_lower or "bostaurus" in label_lower or "holstein" in label_lower:
        return "Cow"
    elif "gallus" in label_lower or "chicken" in label_lower or "caecal" in label_lower:
        return "Chicken"
    elif "human" in label_lower or "dental" in label_lower or "feces" in label_lower:
        return "Human"
    elif "spoiled" in label_lower or "beer" in label_lower:
        return "Other"
    elif "masinjere" in label_lower or "masenjere" in label_lower:
        return "MAG_Masinjere"
    elif "limera" in label_lower:
        return "MAG_Limera"
    else:
        return "Other"

confident["host"] = confident["label"].apply(get_host)

# CAZyme class colors
class_colors = {
    "GH":  "#2166ac",
    "GT":  "#4dac26",
    "PL":  "#d01c8b",
    "CE":  "#f1a340",
    "AA":  "#998ec3",
    "CBM": "#d73027",
    "Other": "#999999"
}

caz_classes = ["GH", "GT", "PL", "CE", "AA", "CBM", "Other"]

# ── GRAPH 1: Stacked bar chart — CAZyme classes per genome ─────────────────────
print("Making Graph 1: Stacked bar per genome...")

genome_class = confident.groupby(["label", "class"]).size().reset_index(name="count")

# Sort genomes by host type then total count
genome_totals = genome_class.groupby("label")["count"].sum().reset_index()
genome_totals["host"] = genome_totals["label"].apply(get_host)
genome_totals = genome_totals.sort_values(["host", "count"], ascending=[True, False])
genome_order = genome_totals["label"].tolist()

fig, ax = plt.subplots(figsize=(24, 9))
bottom = pd.Series([0.0] * len(genome_order), index=genome_order)

for caz_class in caz_classes:
    values = []
    for genome in genome_order:
        subset = genome_class[(genome_class["label"] == genome) & (genome_class["class"] == caz_class)]
        values.append(subset["count"].sum() if len(subset) > 0 else 0)
    ax.bar(range(len(genome_order)), values, bottom=bottom.values,
           label=caz_class, color=class_colors[caz_class], alpha=0.9)
    bottom += pd.Series(values, index=genome_order)

# Add dividing lines between host groups
host_groups = genome_totals.groupby("host").size()
current = 0
host_order_plot = ["Chicken", "Cow", "Human", "Other", "Pig", "MAG_Limera", "MAG_Masinjere"]
prev_host = genome_totals.iloc[0]["host"]
for i, row in genome_totals.iterrows():
    if row["host"] != prev_host:
        ax.axvline(x=current - 0.5, color="black", linestyle="--", linewidth=1.2)
        prev_host = row["host"]
    current += 1

ax.set_xticks(range(len(genome_order)))
ax.set_xticklabels(genome_order, rotation=90, fontsize=5.5)
ax.set_ylabel("Number of CAZyme Annotations", fontsize=12)
ax.set_title("CAZyme Class Distribution per Genome", fontsize=14, fontweight="bold")
ax.legend(title="CAZyme Class", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
plt.tight_layout()
plt.savefig(f"{BASE}/CAZyme_stacked_per_genome.pdf", bbox_inches="tight")
plt.savefig(f"{BASE}/CAZyme_stacked_per_genome.png", dpi=150, bbox_inches="tight")
plt.close()
print("Graph 1 saved!")

# ── GRAPH 2: Grouped bar chart — CAZyme classes by host type ──────────────────
print("Making Graph 2: Grouped bar by host type...")

host_order = ["Cow", "Pig", "Chicken", "Human", "MAG_Masinjere", "MAG_Limera"]
host_class = confident.groupby(["host", "class"]).size().reset_index(name="count")
host_class = host_class[host_class["host"].isin(host_order)]

n_classes = len(caz_classes) - 1  # exclude Other for cleanliness
plot_classes = ["GH", "GT", "PL", "CE", "AA", "CBM"]
width = 0.6 / len(plot_classes)

fig, ax = plt.subplots(figsize=(16, 8))

for i, caz_class in enumerate(plot_classes):
    data = host_class[host_class["class"] == caz_class].set_index("host")
    values = [data.loc[h, "count"] if h in data.index else 0 for h in host_order]
    x_positions = [j + (i - len(plot_classes) / 2) * width + width / 2 for j in range(len(host_order))]
    ax.bar(x_positions, values, width=width,
           label=caz_class, color=class_colors[caz_class], alpha=0.9)

ax.axvline(x=3.5, color="black", linestyle="--", linewidth=1.5)
ax.text(1.5, ax.get_ylim()[1] * 0.97, "NCBI Reference Genomes",
        ha="center", fontsize=10, fontstyle="italic")
ax.text(4.5, ax.get_ylim()[1] * 0.97, "MAG Genomes",
        ha="center", fontsize=10, fontstyle="italic")

ax.set_xticks(range(len(host_order)))
ax.set_xticklabels(["Cow", "Pig", "Chicken", "Human", "Masinjere MAGs", "Limera MAGs"],
                   fontsize=11)
ax.set_ylabel("Total CAZyme Annotations", fontsize=12)
ax.set_title("CAZyme Classes by Host/Source Type", fontsize=14, fontweight="bold")
ax.legend(title="CAZyme Class", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
plt.tight_layout()
plt.savefig(f"{BASE}/CAZyme_grouped_by_host.pdf", bbox_inches="tight")
plt.savefig(f"{BASE}/CAZyme_grouped_by_host.png", dpi=150, bbox_inches="tight")
plt.close()
print("Graph 2 saved!")

# ── GRAPH 3: Heatmap — top CAZyme families across genomes ─────────────────────
print("Making Graph 3: CAZyme family heatmap...")

# Get top 40 most common families
top_families = (confident.groupby("family").size()
                .sort_values(ascending=False)
                .head(40).index.tolist())

fam_genome = confident[confident["family"].isin(top_families)]
fam_counts = fam_genome.groupby(["label", "family"]).size().reset_index(name="count")
fam_pivot = fam_counts.pivot(index="label", columns="family", values="count").fillna(0)

# Sort by host type
fam_pivot["host"] = fam_pivot.index.map(get_host)
fam_pivot["host_order"] = fam_pivot["host"].map({
    "Cow": 0, "Pig": 1, "Chicken": 2, "Human": 3,
    "Other": 4, "MAG_Limera": 5, "MAG_Masinjere": 6
})
fam_pivot = fam_pivot.sort_values("host_order")
fam_pivot = fam_pivot.drop(columns=["host", "host_order"])

# Row colors
def host_color(label):
    host = get_host(label)
    colors = {
        "Cow": "#4dac26", "Pig": "#d01c8b", "Chicken": "#f1a340",
        "Human": "#2166ac", "MAG_Masinjere": "#E87D72",
        "MAG_Limera": "#998ec3", "Other": "#999999"
    }
    return colors.get(host, "#999999")

row_colors = pd.Series([host_color(g) for g in fam_pivot.index], index=fam_pivot.index)

fig, ax = plt.subplots(figsize=(22, 18))
sns.heatmap(fam_pivot, cmap="YlOrRd",
            linewidths=0.2, linecolor="white",
            cbar_kws={"label": "Gene Count", "shrink": 0.4},
            ax=ax)

for i, (genome, color) in enumerate(zip(fam_pivot.index, row_colors)):
    ax.add_patch(mpatches.Rectangle((-1.2, i), 0.8, 1,
                                     color=color, clip_on=False))

ax.set_title("Top 40 CAZyme Families Across All Genomes", fontsize=14, fontweight="bold", pad=20)
ax.set_xlabel("CAZyme Family", fontsize=11)
ax.set_ylabel("Genome", fontsize=11)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)

# Legend
legend_patches = [
    mpatches.Patch(color="#4dac26", label="Cow"),
    mpatches.Patch(color="#d01c8b", label="Pig"),
    mpatches.Patch(color="#f1a340", label="Chicken"),
    mpatches.Patch(color="#2166ac", label="Human"),
    mpatches.Patch(color="#E87D72", label="MAG Masinjere"),
    mpatches.Patch(color="#998ec3", label="MAG Limera"),
]
ax.legend(handles=legend_patches, loc="upper right",
          bbox_to_anchor=(1.15, 1.05), fontsize=9, title="Host/Source")

plt.tight_layout()
plt.savefig(f"{BASE}/CAZyme_family_heatmap.pdf", bbox_inches="tight")
plt.savefig(f"{BASE}/CAZyme_family_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("Graph 3 saved!")

print(f"\nAll done! Files saved to: {BASE}")
