

import os
import csv
import glob
import numpy as np

MANIFEST = "nyu_matched.csv"
SITE_DIR = "NYU"
RUN = "rest_1"
PREFIX = "sfnwmrda"         


def load_1d(path):
    # Carica un .1D Athena -> array (T timepoint x 116 regioni).
    rows = []
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            rows.append([float(x) for x in parts[2:]])
    return np.asarray(rows, dtype=float)


def find_file(subject_id):
    folder = os.path.join(SITE_DIR, subject_id)
    if not os.path.isdir(folder):
        return None
    pattern = os.path.join(folder, f"{PREFIX}{subject_id}_session_1_{RUN}_aal_TCs.1D")
    for h in glob.glob(pattern):
        if "*" not in os.path.basename(h):
            return h
    return None


def main():
    subjects = []
    with open(MANIFEST, newline="") as f:
        for r in csv.DictReader(f):
            subjects.append((r["id"].strip(), r["group"].strip(), r["pair_id"].strip()))

    ids, groups, pairs, mats, skipped = [], [], [], [], []

    for sid, grp, pid in subjects:
        path = find_file(sid)
        if path is None:
            skipped.append((sid, "file mancante o placeholder"))
            continue
        ts = load_1d(path)
        if ts.shape[1] != 116:
            skipped.append((sid, f"colonne={ts.shape[1]}"))
            continue
        if ts.shape[0] < 30:
            skipped.append((sid, f"timepoint={ts.shape[0]}"))
            continue
        C = np.corrcoef(ts.T)
        C = np.nan_to_num(C, nan=0.0)
        np.fill_diagonal(C, 0.0)
        ids.append(sid); groups.append(grp); pairs.append(pid); mats.append(C)

    C_all = np.stack(mats) if mats else np.empty((0, 116, 116))
    np.savez_compressed("nyu_matrices.npz",
                        ids=np.array(ids),
                        groups=np.array(groups),
                        pairs=np.array(pairs),
                        C=C_all)

    n_adhd = sum(1 for g in groups if g == "ADHD")
    n_ctrl = sum(1 for g in groups if g == "control")
    print(f"Caricati: {len(ids)} soggetti  ({n_adhd} ADHD, {n_ctrl} controlli)")
    print(f"Matrici salvate in nyu_matrices.npz  ->  shape {C_all.shape}")
    if skipped:
        print(f"\nScartati: {len(skipped)}")
        for sid, reason in skipped:
            print(f"  {sid}: {reason}")



if __name__ == "__main__":
    main()
