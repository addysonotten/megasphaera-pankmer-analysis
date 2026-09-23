import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re

# Load adjacency matrix
df = pd.read_csv("adj_matrix.csv", index_col=0)

# Load your filled-in metadata
meta = pd.read_csv("rename_filled.csv")
meta = meta.set_index("filename")

def clean(val):
    val = str(val).strip()
    return re.sub(r"[^A-Za-z0-9]+", "", val) if val and val.lower() != "nan" else ""

def build_label(fname):
    if fname not in meta.index:
        return fname[:20]  # fallback if somehow missing from metadata
    row = meta.loc[fname]
    gtype = row.get("type", "")
    if gtype == "NCBI":
        host = clean(row.get("host", ""))
        source = clean(row.get("source", ""))
        parts = [p for p in [host, source] if p]
        label = "NCBI_" + "_".join(parts) if parts else "NCBI_" + fname[:15]
    else:
        village = clean(row.get("village", ""))
        pid = clean(row.get("id", ""))
        age = clean(row.get("age", ""))
        parts = [p for p in [village, pid, age] if p]
        label = "MAG_" + "_".join(parts) if parts else "MAG_" + fname[:15]
    return label

# Build initial labels
raw_labels = [build_label(f) for f in df.index]

# De-duplicate labels (e.g. multiple Human_gut entries) by appending _1, _2, ...
seen = {}
final_labels = []
for lbl in raw_labels:
    if lbl in seen:
        seen[lbl] += 1
        final_labels.append(f"{lbl}_{seen[lbl]}")
    else:
        seen[lbl] = 0
        final_labels.append(lbl)

# Build type list (for color strip) in the SAME order as df.index
types = []
for f in df.index:
    if f in meta.index:
        types.append(meta.loc[f].get("type", "Unknown"))
    else:
        types.append("Unknown")

# Save the lookup table
lookup_df = pd.DataFrame({
    "plot_label": final_labels,
    "original_name": df.index,
    "type": types
})
lookup_df.to_csv("label_lookup_v3.csv", index=False)
print("Saved label_lookup_v3.csv")

# Apply new labels to the matrix
df.index = final_labels
df.columns = final_labels

# Build color strip: NCBI = steelblue, MAG = darkorange
color_map = {"NCBI": "steelblue", "MAG": "darkorange", "Unknown": "gray"}
row_colors = pd.Series(types, index=final_labels).map(color_map)

g = sns.clustermap(
    df,
    figsize=(24, 24),
    cmap="viridis",
    xticklabels=True,
    yticklabels=True,
    row_colors=row_colors,
    col_colors=row_colors
)
g.ax_heatmap.tick_params(labelsize=7)

# Add a legend for the color strip
handles = [mpatches.Patch(color=c, label=t) for t, c in color_map.items() if t in types]
g.ax_heatmap.legend(handles=handles, title="Type", bbox_to_anchor=(1.15, 1), loc="upper left")

g.savefig("clustermap_v3.pdf", dpi=150)
print("Saved clustermap_v3.pdf")
