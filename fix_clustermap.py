import pandas as pd
import seaborn as sns
import re

df = pd.read_csv("adj_matrix.csv", index_col=0)

def make_label(name, counter):
    match = re.match(r"^(GC[AF]_\d+\.\d+)", str(name))
    if match:
        return match.group(1), counter
    return f"MAG_{counter:02d}", counter + 1

new_labels = []
lookup = []
counter = 1
for name in df.index:
    label, counter = make_label(name, counter)
    new_labels.append(label)
    lookup.append({"plot_label": label, "original_name": name})

df.index = new_labels
df.columns = new_labels

pd.DataFrame(lookup).to_csv("label_lookup.csv", index=False)

g = sns.clustermap(df, figsize=(22, 22), cmap="viridis", xticklabels=True, yticklabels=True)
g.ax_heatmap.tick_params(labelsize=8)
g.savefig("clustermap_v2_fixed.pdf", dpi=150)
print("Done!")
