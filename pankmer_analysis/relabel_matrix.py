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
    # keeps digits and a single decimal point (e.g. "2.5" stays "2.5", turned into "2p5" for label safety)
    val = str(val).strip()
    if not val or val.lower() == "nan":
        return ""
    val = re.sub(r"[^0-9.]", "", val)
    return val  # dots can be awkward in labels, so use "p" instead (2.5 -> 2p5)

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
print("Saved label_lookup_v3.csv")

matched = sum(1 for f in df.index if strip_ext(f) in meta.index)
print(f"Matched {matched} of {len(df.index)} genomes to metadata")

df.index = final_labels
df.columns = final_labels
df.to_csv("adj_matrix_relabeled.csv")
print("Saved adj_matrix_relabeled.csv")
