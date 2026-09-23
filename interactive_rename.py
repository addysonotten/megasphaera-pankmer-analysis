import os
import csv
import re

genome_dir = "all_genomes"
output_file = "rename_filled.csv"

files = sorted(os.listdir(genome_dir))
fieldnames = ["filename", "type", "host", "source", "village", "age", "id"]

# Resume support: skip files already done if you stop partway through
done = {}
if os.path.exists(output_file):
    with open(output_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            done[row["filename"]] = row

rows = []
total = len(files)

def is_ncbi(fname):
    return re.match(r"^GC[AF]_\d+\.\d+", fname) is not None

for i, fname in enumerate(files, 1):
    if fname in done:
        rows.append(done[fname])
        continue

    print(f"\n[{i}/{total}] {fname}")

    if is_ncbi(fname):
        host = input("  Host (e.g. Human, Pig, Cow) [blank to skip]: ").strip()
        source = input("  Source (e.g. gut, feces, rumen) [blank to skip]: ").strip()
        row = {"filename": fname, "type": "NCBI", "host": host, "source": source,
               "village": "", "age": "", "id": ""}
    else:
        village = input("  Village [blank to skip]: ").strip()
        age = input("  Age [blank to skip]: ").strip()
        pid = input("  ID [blank to skip]: ").strip()
        row = {"filename": fname, "type": "MAG", "host": "", "source": "",
               "village": village, "age": age, "id": pid}

    rows.append(row)

    # Save progress after every entry, so nothing is lost if you stop early
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

print(f"\nDone! Saved to {output_file}")
