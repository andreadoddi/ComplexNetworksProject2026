# ComplexNetworksProject2026
# Functional brain network organization in ADHD: how robust are case-control differences?

A graph-theoretic case–control analysis of resting-state functional brain
networks in children with ADHD, built on the **ADHD-200** dataset. The project
asks not only whether ADHD and control networks differ topologically, but
whether any such difference is **robust** to the methodological choices involved
in constructing and comparing the networks (thresholding, group matching,
multiple-comparison correction, and control for overall connectivity).

**Main result:** no statistically robust topological difference between the two
groups. Weak trends in global efficiency and characteristic path length, in the
direction reported by earlier studies, do not survive FDR correction and are
fully accounted for by individual variation in overall functional connectivity.

---

## Overview

For each subject, a functional connectivity matrix is obtained from the Pearson
correlations between the mean BOLD time series of the 116 AAL regions. Each
matrix is thresholded into a binary graph (proportional thresholding, positive
correlations only) across a range of densities, and five standard graph metrics
are computed. Groups are compared with a paired, sex- and age-matched design,
using non-parametric tests, FDR correction, and a control for overall
functional connectivity following van den Heuvel et al. (2017).

Pipeline at a glance:

```
phenotypic TSV + NYU/ time-series (.1D)
        │
        ▼
  usable_nyu.py        → nyu_usable.csv      (190 usable subjects)
        │
        ▼
  select_and_match.py  → nyu_matched.csv     (64 matched pairs, Hungarian)
        │
        ▼
  build_matrices.py    → nyu_matrices.npz    (128 × 116×116 correlation matrices)
        │
        ▼
  analyze.py           → metrics_by_density.csv + p/q tables
```

---

## Data

The analysis uses the **ADHD-200 Preprocessed** repository (Athena pipeline),
which provides, for each subject, regional mean BOLD time series already
extracted with the AAL atlas (116 regions).

- Repository (NITRC): https://www.nitrc.org/frs/?group_id=383
- Phenotypic key (site / diagnosis / QC codes):
  https://fcon_1000.projects.nitrc.org/indi/adhd200/general/ADHD-200_PhenotypicKey.pdf

Files needed locally (not included in this repo — download from NITRC):

- `adhd200_preprocessed_phenotypics.tsv` — phenotypic table (labels, QC, demographics)
- `ADHD200_AAL_TCs_filtfix` — band-pass filtered AAL time courses; after
  extraction this gives one folder per site. Only the `NYU/` folder is used.

### Cohort selection

- **Single site:** NYU (`Site == 5`), to remove between-site batch effects.
- **Quality control:** keep only `QC_Athena == 1` (pass).
- **Diagnosis:** ADHD = `DX ∈ {1, 2, 3}` (combined / hyperactive-impulsive /
  inattentive); control = `DX == 0` (typically developing); `pending` discarded.
- **Time series:** for the filtered variant, use `sfnwmrda*` files, first run
  (`rest_1`) only.

This yields **190 usable subjects** (87 controls, 103 ADHD). Sex-stratified
one-to-one matching (see below) then produces **64 pairs = 128 subjects**.

---

## Requirements

- Python 3.10+
- `numpy`, `scipy`, `networkx`, `matplotlib`

Optional, for the anatomical brain figures:

- `nilearn`, `nibabel`, `plotly`

```bash
pip install numpy scipy networkx matplotlib
pip install nilearn plotly          # optional, for glass-brain figures
```

---

## Usage

Run from the folder containing the phenotypic TSV and the extracted `NYU/`
folder. Each step writes the input for the next.

```bash
python usable_nyu.py         # → nyu_usable.csv
python select_and_match.py   # → nyu_matched.csv
python build_matrices.py     # → nyu_matrices.npz
python analyze.py            # → metrics_by_density.csv + prints p/q tables
```

---

## Scripts

### Core pipeline

| Script | Input | Output | What it does |
|---|---|---|---|
| `usable_nyu.py` | phenotypic TSV, `NYU/` | `nyu_usable.csv` | Intersects on-disk subjects with `Site==5` and `QC==1`; assigns group. |
| `select_and_match.py` | `nyu_usable.csv` | `nyu_matched.csv` | Optimal 1:1 matching (Hungarian algorithm), within sex, minimizing total age difference; keeps pairs with age gap ≤ 1 year. |
| `build_matrices.py` | `nyu_matched.csv`, `.1D` files | `nyu_matrices.npz` | Loads each subject's regional time series, computes the 116×116 Pearson correlation matrix. |
| `analyze.py` | `nyu_matrices.npz` | `metrics_by_density.csv` | Proportional thresholding across densities; computes five graph metrics; paired Wilcoxon tests (raw and FC-corrected) with FDR correction. |



---

## Method notes

- **Thresholding.** Proportional (density-based) thresholding keeps the same
  number of edges per subject, enabling fair comparison. Because the choice of
  density is arbitrary, all metrics are reported across a range (5–30%).
- **Overall FC confound.** Proportional thresholding interacts with overall
  connectivity: subjects with lower overall FC accumulate more spurious edges,
  inflating measures like global efficiency. Following van den Heuvel et al.
  (2017), overall FC is tested between groups and regressed out of each paired
  comparison as a confound.
- **Matching.** Sex-stratified optimal matching guarantees identical sex
  composition and minimal residual age difference. The strong sex imbalance of
  the ADHD group limits the number of female pairs; NYU was chosen over the
  larger Peking site because it allows a more balanced number of female pairs.
- **Statistics.** Paired Wilcoxon signed-rank tests throughout (no normality
  assumption), with Benjamini–Hochberg FDR correction over the 5 metrics × 6
  densities.

---

## Key references

-ADHD-200 Consortium. The ADHD-200 Consortium: A Model to Advance the Translational Potential of Neuroimaging in Clinical Neuroscience. Frontiers in Systems   Neuroscience, 6:62, 2012. doi:10.3389/fnsys.2012.00062.


-ADHD-200 Consortium. ADHD-200 Preprocessed Repository.Neuroimaging Informatics Tools and Resources Clearinghouse (NITRC), Preprocessed Connectomes Project. Available at:https://www.nitrc.org/frs/?group_id=383. Accessed July 2026.


-ADHD-200 Consortium. ADHD-200 Phenotypic Key.Supplementary documentation for the ADHD-200 Sample. Neuroimaging Informatics Tools and Resources Clearinghouse (NITRC). Available at: https://fcon_1000.projects.nitrc.org/indi/adhd200/general/ADHD-200_PhenotypicKey.pdf. Accessed July 2026.

-Fornito A., Zalesky A., Bullmore E.T. Fundamentals of Brain Network Analysis. Academic Press, 2016.

-Van den Heuvel M.P., de Lange S.C., Zalesky A., Seguin C., Yeo B.T.T., Schmidt R. Proportional thresholding in resting-state fMRI functional connectivity networks and consequences for patient-control connectome studies: Issues and recommendations. NeuroImage 152:437--449, 2017.

-Wang L., Zhu C., He Y., Zang Y., Cao Q., Zhang H., Zhong Q., Wang Y. Altered small-world brain functional networks in children with attention-deficit/hyperactivity disorder. Human Brain Mapping 30(2):638--649, 2009.

-Zalesky A., Fornito A., Bullmore E.T. Network-based statistic: identifying differences in brain networks. NeuroImage 53(4):1197--1207, 2010.



---

## Notes

- The raw ADHD-200 data are not redistributed here; download them from NITRC.
- Results are reproducible: community detection and network layouts use fixed
  random seeds.
