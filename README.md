# Interactive UMAP + HDBSCAN maps for genome datasets

> **Status: early prototype.** Interfaces and defaults may change.

Turn genome embeddings, or UMAP/HDBSCAN results you have already computed, plus any annotations you have into a single interactive map. The output is a standalone HTML file that opens in any browser and can be emailed to collaborators or attached as supplementary material.

- **Hover** over a genome to see its ID, dataset, cluster and all of your annotations.
- **Colour by** clusters, datasets, or any annotation column. Categories such as genus get separate colours; numbers such as GC content get a colour scale.
- **Click legend entries** to show or hide clusters and datasets.
- **Smaller datasets stay visible.** Any dataset making up less than 40% of all genomes is drawn on top of the larger ones with a thin black outline, so it doesn't get hidden.
- **Bring your own analysis.** The tool doesn't need to know what your annotations mean. It displays whatever columns you give it.

## Two ways to use it

| You have… | Use | The tool… |
|---|---|---|
| **2-D coordinates** (and optionally clusters) from your own UMAP/HDBSCAN, t-SNE, PCA, etc. | `--precomputed` | builds the interactive map only |
| **Embeddings** (e.g. from DNABERT-S) | `--input` | runs UMAP + HDBSCAN, checks cluster stability, then builds the map |

Either way you can add `--annotations` and `--color-by`.

## Install

Requires Python 3.9+.

```bash
pip install -r requirements.txt
```

## Quick start (example data included)

The `examples/` folder contains a small **synthetic** dataset: 500 "gut" and 180 "vaginal" genomes with made-up embeddings, plus an annotation file with genus, GC content and genome size.

```bash
python interactive_umap_hdbscan.py --input Gut examples/example_gut_embeddings.csv --input Vaginal examples/example_vaginal_embeddings.csv --annotations examples/example_annotations.csv --out results
```

Open `results/interactive_umap.html` in a browser. To colour by genus or GC content, add `--color-by Genus` or `--color-by GC_content`.

## Input files

All inputs are CSV files with **genome IDs in the first column**.

**Precomputed coordinates** (`--precomputed`): one table for all genomes.

```
Genome_ID,UMAP_1,UMAP_2,Cluster,Dataset
GENOME_A,3.21,-1.05,0,Gut
GENOME_B,2.98,-0.87,0,Vaginal
```

- `UMAP_1`, `UMAP_2` are required. If your columns are named differently, use `--x-col` and `--y-col`.
- `Cluster` and `Dataset` are optional. If named differently, use `--cluster-col` and `--dataset-col`.
- Any other columns are treated as annotations.

**Embeddings** (`--input LABEL file.csv`): one file per dataset. Every column after the ID is a numeric embedding dimension. All files must have the same number of dimensions.

```
,dim_0,dim_1,dim_2,...
GENOME_A,-0.0619,0.1028,0.0465,...
```

**Annotations** (`--annotations`, optional): any per-genome information, matched to the main data by genome ID. Genomes with no annotation show `NA`.

```
Genome_ID,Genus,GC_content,Isolation_source
GENOME_A,Bacteroides,43.1,Stool
GENOME_B,Lactobacillus,36.4,Vaginal swab
```

## Examples

```bash
# You already have UMAP + HDBSCAN results and some metadata
python interactive_umap_hdbscan.py --precomputed my_umap.csv --annotations metadata.csv --color-by Genus

# Your coordinate columns have different names
python interactive_umap_hdbscan.py --precomputed my_tsne.csv --x-col tSNE_1 --y-col tSNE_2 --dataset-col Site

# Start from embeddings of three datasets
python interactive_umap_hdbscan.py --input Gut gut.csv --input Vaginal vaginal.csv --input Urinary urinary.csv --title "Gut vs. vaginal vs. urinary" --out results
```

On Windows, paths containing spaces must be in quotes, e.g. `--input Gut "C:\My Data\gut.csv"`.

## Outputs

Written to the `--out` folder (default `results/`):

| File | Contents |
|---|---|
| `interactive_umap.html` | The interactive map (self-contained; works offline) |
| `report.txt` | Dataset sizes, genomes per cluster and, for embedding input, a cluster stability summary |
| `coordinates_and_clusters.csv` | *Embedding input only.* Coordinates, cluster, dataset and annotations for every genome. Can be reused later with `--precomputed`, e.g. to recolour without rerunning UMAP. |

Cluster `-1` means *unclustered*: HDBSCAN did not assign the genome to any cluster because it is in a sparse area. These genomes are shown in grey.

## Options

| Option | Default | Meaning |
|---|---|---|
| `--annotations` | – | Per-genome annotation CSV. |
| `--color-by` | `Cluster` (else `Dataset`) | Column to colour by. Numeric columns with many values get a colour scale; others get separate colours. |
| `--max-categories` | `20` | When colouring by categories, the most common get their own colour and the rest are grouped as "Other". |
| `--title` | `Interactive UMAP + HDBSCAN` | Plot and page title. |
| `--minority-fraction` | `0.40` | Datasets below this share of genomes are drawn on top with an outline. |
| `--out` | `results` | Output folder. |

Embedding input only:

| Option | Default | Meaning |
|---|---|---|
| `--min-cluster-fraction` | `0.05` | Smallest cluster, as a fraction of all genomes. The default gives a coarse first-pass view (usually fewer than 10 clusters). Lower it (e.g. `0.017`) for finer sub-clusters. |
| `--min-samples` | `5` | How strictly HDBSCAN requires dense regions. Values of 15 or more can merge most genomes into one or two clusters. |
| `--n-neighbors` | `30` | UMAP neighbourhood size. Larger values emphasise global structure. |
| `--min-dist` | `0.1` | UMAP point spacing. |
| `--seed` | `42` | Random seed for the main result. |
| `--stability-runs` | `5` | Number of seeds used for the stability check (`0` to skip; faster). |

Run `python interactive_umap_hdbscan.py --help` for the full list.

## How the embedding workflow works

1. **Combine** all embedding tables into one matrix, keeping track of which dataset each genome came from.
2. **UMAP** projects the combined embeddings to 2-D (Euclidean distance).
3. **HDBSCAN** clusters the 2-D map, so every cluster matches what you see in the plot.
4. **Stability check:** steps 2–3 are repeated with different random seeds. Agreement between runs is measured with the Adjusted Rand Index (1 = identical clusters).
   - 0.9 or higher: consistent
   - 0.75 to 0.9: mostly consistent
   - below 0.75: varied; treat cluster assignments with caution

## Limitations

- UMAP is stochastic. Cluster boundaries can shift between runs or when genomes are added or removed. For embedding input, the stability report tells you how much this matters for your data.
- A dataset with fewer genomes than the minimum cluster size cannot form a cluster of its own. It can only join other clusters or remain unclustered.
- UMAP distances between far-apart groups are not meaningful. Treat the map as a view of local neighbourhoods, not as exact distances.
- For precomputed input, clusters and coordinates are shown as provided; the tool does not check them.

## Built with

This tool is a thin workflow around the following methods and open-source libraries. If you use it, please also cite the methods it relies on.

- **UMAP** (dimensionality reduction):
  - McInnes L, Healy J, Melville J. *UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.* arXiv:1802.03426 (2018).
  - Software: [umap-learn](https://github.com/lmcinnes/umap). McInnes L, Healy J, Saul N, Großberger L. *UMAP: Uniform Manifold Approximation and Projection.* Journal of Open Source Software 3(29):861 (2018).
- **HDBSCAN** (clustering):
  - Campello RJGB, Moulavi D, Sander J. *Density-Based Clustering Based on Hierarchical Density Estimates.* PAKDD 2013, LNCS 7819:160–172.
  - Software: [hdbscan](https://github.com/scikit-learn-contrib/hdbscan). McInnes L, Healy J, Astels S. *hdbscan: Hierarchical density based clustering.* Journal of Open Source Software 2(11):205 (2017).
- **scikit-learn** (Adjusted Rand Index): Pedregosa F, et al. *Scikit-learn: Machine Learning in Python.* Journal of Machine Learning Research 12:2825–2830 (2011).
- **Plotly** (interactive plots): [plotly.py](https://github.com/plotly/plotly.py), Plotly Technologies Inc.
- **pandas** and **NumPy** (data handling).

The tool works with any numeric embedding. The example use case is genome embeddings from **DNABERT-S**: Zhou Z, et al. *DNABERT-S: Pioneering Species Differentiation with Species-Aware DNA Embeddings.* arXiv:2402.08777 (2024).

## Citation

Manuscript in preparation. Please check back for citation details.
