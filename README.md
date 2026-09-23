# Interactive UMAP + HDBSCAN for genome embeddings

> **Status: early prototype.** Interfaces and defaults may change.

Turn one or more sets of genome embeddings (e.g. from DNABERT-S) into a single interactive map. Genomes are projected to 2-D with [UMAP](https://umap-learn.readthedocs.io/), grouped with [HDBSCAN](https://hdbscan.readthedocs.io/), and written to a standalone HTML file you can open in any browser. You can hover over a point to see the genome ID, and click legend entries to show or hide clusters and datasets.

Features:

- Accepts any number of datasets, each with its own label and marker shape.
- **Smaller datasets are drawn on top.** Any dataset making up less than 40% of all genomes is overlaid on the larger ones with a thin black outline, so it doesn't get hidden.
- **Stability check.** The analysis is repeated with several random seeds, and a plain-language report says how consistent the clusters were.
- **Coarse clusters by default** (usually fewer than 10) for a first-pass overview. Lower `--min-cluster-fraction` for finer sub-clusters.

## Install

Requires Python 3.9+.

```bash
pip install -r requirements.txt
```

## Input format

One CSV per dataset. The first column holds the genome IDs. The remaining columns are the numeric embedding. All files must have the same number of embedding dimensions, but the column headers don't need to match.

```
,dim_0,dim_1,dim_2,...
GENOME_A,-0.0619,0.1028,0.0465,...
GENOME_B,0.0293,0.1183,0.0155,...
```

## Usage

```bash
python interactive_umap_hdbscan.py \
    --input Gut      gut_embeddings.csv \
    --input Vaginal  vaginal_embeddings.csv \
    --input Urinary  urinary_embeddings.csv \
    --title "Gut vs. vaginal vs. urinary genomes" \
    --out results
```

Outputs in `results/`:

| File | Contents |
|---|---|
| `interactive_umap.html` | The interactive map (self-contained; no internet needed) |
| `coordinates_and_clusters.csv` | UMAP coordinates, cluster and dataset for every genome |
| `stability_report.txt` | Dataset sizes, cluster counts, stability summary, genomes per cluster |

Cluster `-1` means *unclustered*: the genome is in a sparse area and HDBSCAN did not assign it to any cluster. These points are shown in grey.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--min-cluster-fraction` | `0.05` | Smallest cluster, as a fraction of all genomes. Lower it (e.g. `0.017`) for more, smaller clusters. |
| `--min-samples` | `5` | How strictly HDBSCAN requires dense regions. Values of 15 or more can merge most genomes into one or two clusters. |
| `--n-neighbors` | `30` | UMAP neighbourhood size. Larger values emphasise global structure. |
| `--min-dist` | `0.1` | UMAP point spacing. |
| `--seed` | `42` | Random seed for the main result. |
| `--stability-runs` | `5` | Number of seeds used for the stability check (`0` to skip; faster). |
| `--minority-fraction` | `0.40` | Datasets below this share of genomes are drawn on top with an outline. |

Run `python interactive_umap_hdbscan.py --help` for the full list.

## How it works

1. **Combine** all embedding tables into one matrix, keeping track of which dataset each genome came from.
2. **UMAP** projects the combined embeddings to 2-D (Euclidean distance).
3. **HDBSCAN** clusters the 2-D map, so every cluster matches what you see in the plot.
4. **Stability check:** steps 2–3 are repeated with different random seeds. Agreement between runs is measured with the Adjusted Rand Index (1 = identical clusters).
   - 0.9 or higher: consistent
   - 0.75 to 0.9: mostly consistent
   - below 0.75: varied; treat cluster assignments with caution

## Limitations

- UMAP is stochastic. Cluster boundaries can shift between runs or when genomes are added or removed; the stability report tells you how much this matters for your data.
- A dataset with fewer genomes than the minimum cluster size cannot form a cluster of its own. It can only join other clusters or remain unclustered.
- UMAP distances between far-apart groups are not meaningful. Treat the map as a view of local neighbourhoods, not as exact distances.

## Citation

Manuscript in preparation. Please check back for citation details.
