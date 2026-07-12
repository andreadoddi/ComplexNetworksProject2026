
import os
import csv
from collections import Counter

TSV = "adhd200_preprocessed_phenotypics.tsv"
SITE_DIR = "NYU"


disk_ids = set()
for name in os.listdir(SITE_DIR):
    if os.path.isdir(os.path.join(SITE_DIR, name)) and name.isdigit():
        disk_ids.add(int(name))

mappa ID -> info dal TSV, solo Site==5 e QC==1
tsv_info = {}
id_col = None
with open(TSV, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for cand in ("ScanDir ID", "ScanDirID", "ScanDir_ID", "ID"):
        if cand in reader.fieldnames:
            id_col = cand
            break
    if id_col is None:
        print("Colonne:", reader.fieldnames)
        raise SystemExit("Colonna ID non trovata.")
    for row in reader:
        raw = row[id_col].strip()
        if not raw.isdigit():
            continue
        if row["Site"].strip() != "5":
            continue
        if row["QC_Athena"].strip() != "1":
            continue
        tsv_info[int(raw)] = row


usable = sorted(disk_ids & set(tsv_info))
print(f"Cartelle su disco:                {len(disk_ids)}")
print(f"Validi nel TSV (Site==5, QC==1):  {len(tsv_info)}")
print(f"UTILIZZABILI (intersezione):      {len(usable)}")

 
dx = Counter(tsv_info[i]["DX"].strip() for i in usable)
ctrl = sum(1 for i in usable if tsv_info[i]["DX"].strip() == "0")
adhd = sum(1 for i in usable if tsv_info[i]["DX"].strip() in ("1", "2", "3"))
print("\nDistribuzione DX tra gli utilizzabili:")
for k, v in sorted(dx.items()):
    print(f"  DX = {k:<8} -> {v}")
print(f"\nControlli: {ctrl}    ADHD: {adhd}")


final = [i for i in usable if tsv_info[i]["DX"].strip() in ("0", "1", "2", "3")]


with open("nyu_usable.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "group", "dx", "age", "gender"])
    for i in final:
        r = tsv_info[i]
        grp = "control" if r["DX"].strip() == "0" else "ADHD"
        gender = "M" if r["Gender"].strip() == "1" else "F"
        w.writerow([f"{i:07d}", grp, r["DX"].strip(), r["Age"].strip(), gender])
print(f"\nSalvato nyu_usable.csv  ({len(final)} soggetti con gruppo definito)")
