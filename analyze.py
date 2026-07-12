

import numpy as np
import networkx as nx
from scipy import stats
from collections import defaultdict
import csv

NPZ = "nyu_matrices.npz"
DENSITIES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
METRIC_NAMES = ["eff_glob", "clustering", "eff_loc", "modularity", "path_len"]
FDR_Q = 0.05   # livello FDR (proporzione attesa di falsi positivi tollerata)


def threshold_binary(M, density):
    M = M.copy()
    M[M < 0] = 0
    iu = np.triu_indices_from(M, k=1)
    w = M[iu]
    n_keep = int(round(density * len(w)))
    if n_keep < 1:
        return None
    thr = np.sort(w)[::-1][n_keep - 1]
    A = (M >= thr).astype(int)
    np.fill_diagonal(A, 0)
    return A


def _char_path_len(G):
    if G.number_of_edges() == 0:
        return np.nan
    comp = max(nx.connected_components(G), key=len)
    sub = G.subgraph(comp)
    if sub.number_of_nodes() < 2:
        return np.nan
    return nx.average_shortest_path_length(sub)


def compute_metrics(A):
    G = nx.from_numpy_array(A)
    comms = nx.community.louvain_communities(G, seed=0)
    return {
        "eff_glob": nx.global_efficiency(G),
        "clustering": nx.average_clustering(G),
        "eff_loc": nx.local_efficiency(G),
        "modularity": nx.community.modularity(G, comms),
        "path_len": _char_path_len(G),
    }


def overall_fc(M):
    iu = np.triu_indices_from(M, k=1)
    w = M[iu]
    return w[w > 0].mean()


def paired_pvalues(av, cv, d_ofc):
    
    ok = ~(np.isnan(av) | np.isnan(cv))
    dm = av[ok] - cv[ok]
    try:
        _, p_raw = stats.wilcoxon(dm)
    except ValueError:
        p_raw = np.nan
    do = d_ofc[ok]
    if len(dm) > 2 and np.std(do) > 0:
        beta = np.polyfit(do, dm, 1)
        resid = dm - np.polyval(beta, do)
        try:
            _, p_corr = stats.wilcoxon(resid)
        except ValueError:
            p_corr = np.nan
    else:
        p_corr = np.nan
    return p_raw, p_corr


def fdr(pvals):

    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    mask = ~np.isnan(p)
    if mask.sum() == 0:
        return out
    try:
        out[mask] = stats.false_discovery_control(p[mask], method="bh")
    except AttributeError:
        # fallback manuale BH
        idx = np.where(mask)[0]
        order = idx[np.argsort(p[idx])]
        m = len(order)
        prev = 1.0
        for rank in range(m - 1, -1, -1):
            i = order[rank]
            q = p[i] * m / (rank + 1)
            prev = min(prev, q)
            out[i] = prev
    return out


def main():
    d = np.load(NPZ, allow_pickle=True)
    C, groups, pairs, ids = d["C"], d["groups"], d["pairs"], d["ids"]
    n = len(ids)
    print(f"Caricati {n} soggetti  ({np.sum(groups=='ADHD')} ADHD, "
          f"{np.sum(groups=='control')} controlli)\n")

    ofc = np.array([overall_fc(C[k]) for k in range(n)])
    # indici di coppia (servono sia per il test overall FC sia per le metriche)
    pair_idx = defaultdict(dict)
    for k in range(n):
        pair_idx[pairs[k]][groups[k]] = k
    complete = [p for p, dd in pair_idx.items()
                if "ADHD" in dd and "control" in dd]
    a_k = np.array([pair_idx[p]["ADHD"] for p in complete])
    c_k = np.array([pair_idx[p]["control"] for p in complete])
   
    d_ofc = ofc[a_k] - ofc[c_k]
    _, p_ofc = stats.wilcoxon(d_ofc)
    print("=== Overall FC (confondente) ===")
    print(f"ADHD={ofc[a_k].mean():.4f}  controlli={ofc[c_k].mean():.4f}  "
          f"Wilcoxon appaiato p={p_ofc:.4f}")
    print("  -> " + ("differisce: correzione FC NECESSARIA" if p_ofc < 0.05
                     else "non differisce: artefatto FC limitato") + "\n")

    results = {dd: {m: np.full(n, np.nan) for m in METRIC_NAMES} for dd in DENSITIES}
    rows_out = []
    for k in range(n):
        for dd in DENSITIES:
            A = threshold_binary(C[k], dd)
            if A is None:
                continue
            mt = compute_metrics(A)
            for m in METRIC_NAMES:
                results[dd][m][k] = mt[m]
            rows_out.append([ids[k], groups[k], pairs[k], dd, ofc[k]] +
                            [mt[m] for m in METRIC_NAMES])
    with open("metrics_by_density.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "group", "pair_id", "density", "overall_fc"] + METRIC_NAMES)
        w.writerows(rows_out)
    print("Salvato metrics_by_density.csv\n")

    d_ofc = ofc[a_k] - ofc[c_k]
    print(f"Coppie complete: {len(complete)}\n")

    # raccogli p grezzi e corretti-FC per ogni (densita', metrica)
    raw = {dd: {} for dd in DENSITIES}
    corr = {dd: {} for dd in DENSITIES}
    for dd in DENSITIES:
        for m in METRIC_NAMES:
            vals = results[dd][m]
            p_r, p_c = paired_pvalues(vals[a_k], vals[c_k], d_ofc)
            raw[dd][m] = p_r
            corr[dd][m] = p_c

    #  FDR applicato sull'insieme dei test (tutte densita' x metriche) 
    flat_keys = [(dd, m) for dd in DENSITIES for m in METRIC_NAMES]
    raw_flat = np.array([raw[dd][m] for dd, m in flat_keys])
    corr_flat = np.array([corr[dd][m] for dd, m in flat_keys])
    q_raw = fdr(raw_flat)
    q_corr = fdr(corr_flat)
    q_raw_d = {k: v for k, v in zip(flat_keys, q_raw)}
    q_corr_d = {k: v for k, v in zip(flat_keys, q_corr)}

    def fmt(x):
        return "  nan " if np.isnan(x) else f"{x:.3f}"

    print(f"=== Test APPAIATO (Wilcoxon) — p grezzi ===")
    print("densita'" + "".join(f"{m:>12}" for m in METRIC_NAMES))
    for dd in DENSITIES:
        print(f"{dd:>6.0%}  " + "".join(f"{fmt(raw[dd][m]):>12}" for m in METRIC_NAMES))

    print(f"\n=== q-value FDR (Benjamini-Hochberg, q={FDR_Q}) sui p grezzi ===")
    print("(significativo se q < %.2f; marcato con *)" % FDR_Q)
    print("densita'" + "".join(f"{m:>12}" for m in METRIC_NAMES))
    for dd in DENSITIES:
        cells = []
        for m in METRIC_NAMES:
            q = q_raw_d[(dd, m)]
            star = "*" if (not np.isnan(q) and q < FDR_Q) else " "
            cells.append(f"{fmt(q)}{star}".rjust(12))
        print(f"{dd:>6.0%}  " + "".join(cells))

    print(f"\n=== q-value FDR sui p corretti per overall FC ===")
    print("densita'" + "".join(f"{m:>12}" for m in METRIC_NAMES))
    for dd in DENSITIES:
        cells = []
        for m in METRIC_NAMES:
            q = q_corr_d[(dd, m)]
            star = "*" if (not np.isnan(q) and q < FDR_Q) else " "
            cells.append(f"{fmt(q)}{star}".rjust(12))
        print(f"{dd:>6.0%}  " + "".join(cells))

    n_sig_raw = np.nansum(q_raw < FDR_Q)
    n_sig_corr = np.nansum(q_corr < FDR_Q)
    print(f"\nTest significativi dopo FDR: {int(n_sig_raw)} (grezzi), "
          f"{int(n_sig_corr)} (corretti per FC), su {len(flat_keys)} test.")



if __name__ == "__main__":
    main()
