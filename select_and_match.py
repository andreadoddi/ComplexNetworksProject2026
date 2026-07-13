

import csv
import numpy as np
from scipy.optimize import linear_sum_assignment

USABLE = "nyu_usable.csv"
AGE_TOL = 1.0   


def match_within(subset):

    adhd = [r for r in subset if r["group"] == "ADHD"]
    ctrl = [r for r in subset if r["group"] == "control"]
    if not adhd or not ctrl:
        return []
    A = np.array([r["age"] for r in adhd])
    C = np.array([r["age"] for r in ctrl])
    cost = np.abs(A[:, None] - C[None, :])       # |eta_i - eta_j|
    row_idx, col_idx = linear_sum_assignment(cost)  # min somma dei costi
    return [(adhd[i], ctrl[j], cost[i, j]) for i, j in zip(row_idx, col_idx)]


def main():
    rows = list(csv.DictReader(open(USABLE, newline="")))
    for r in rows:
        r["age"] = float(r["age"])

    pairs = []
    for sex in ("M", "F"):
        sub = [r for r in rows if r["gender"].strip() == sex]
        pairs += match_within(sub)

    # tieni solo coppie entro la tolleranza d'eta
    kept = [(a, c, d) for a, c, d in pairs if d <= AGE_TOL]
    kept.sort(key=lambda t: (t[0]["gender"], t[0]["age"]))

    # report
    dif = [d for *_, d in kept]
    nM = sum(1 for a, _, _ in kept if a["gender"] == "M")
    nF = len(kept) - nM
    print(f"Coppie formate: {len(kept)}  ({nM} maschili, {nF} femminili)")
    print(f"Soggetti totali: {2*len(kept)}  ({len(kept)} ADHD + {len(kept)} controlli)")
    print(f"Scarto d'eta: medio={np.mean(dif):.3f}  max={np.max(dif):.3f} anni")

    
    with open("nyu_matched.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pair_id", "id", "group", "dx", "age", "gender"])
        for k, (a, c, d) in enumerate(kept):
            for r in (a, c):
                w.writerow([k, r["id"], r["group"], r["dx"],
                            f"{r['age']:.2f}", r["gender"]])
    print(f"\nSalvato nyu_matched.csv  ({2*len(kept)} righe, {len(kept)} coppie)")


if __name__ == "__main__":
    main()
