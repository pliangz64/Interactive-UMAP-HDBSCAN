"""Interactive UMAP + HDBSCAN explorer for genome embeddings.

Combines one or more embedding tables (one row per genome), projects them with UMAP,
clusters the 2-D map with HDBSCAN, and writes a standalone interactive HTML plot.

Example:
    python interactive_umap_hdbscan.py \
        --input UHGG UHGG_embeddings.csv \
        --input VMGC vmgc_embeddings.csv \
        --out results
"""
import argparse
import itertools
import os
import sys

import hdbscan
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import umap
from sklearn.metrics import adjusted_rand_score

# Marker shapes per dataset: largest dataset first, then minority datasets
SYMBOLS = ['circle', 'diamond', 'square', 'triangle-up', 'cross', 'x', 'star', 'pentagon']
NOISE_COLOR = '#c8c8c8'


def parse_args():
    p = argparse.ArgumentParser(
        description="Interactive UMAP + HDBSCAN explorer for genome embeddings.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", nargs=2, action="append", metavar=("LABEL", "CSV"), required=True,
                   help="Dataset label and embedding CSV (first column = genome ID, remaining "
                        "columns = numeric embedding). Repeat for each dataset.")
    p.add_argument("--out", default="results", help="Output directory.")
    p.add_argument("--title", default="Interactive UMAP + HDBSCAN", help="Plot and page title.")
    p.add_argument("--min-cluster-fraction", type=float, default=0.05,
                   help="HDBSCAN min_cluster_size as a fraction of all genomes. Lower (e.g. 0.017) "
                        "for finer sub-clusters.")
    p.add_argument("--min-samples", type=int, default=5,
                   help="HDBSCAN min_samples. Values >= 15 can collapse 2-D maps into 1-2 clusters.")
    p.add_argument("--n-neighbors", type=int, default=30, help="UMAP n_neighbors.")
    p.add_argument("--min-dist", type=float, default=0.1, help="UMAP min_dist.")
    p.add_argument("--seed", type=int, default=42, help="Random seed for the main result.")
    p.add_argument("--stability-runs", type=int, default=5,
                   help="Number of UMAP seeds used to check cluster stability (0 or 1 to skip).")
    p.add_argument("--minority-fraction", type=float, default=0.40,
                   help="Datasets below this fraction of all genomes are drawn on top with an outline.")
    return p.parse_args()


def load_datasets(inputs):
    frames, labels = [], []
    for label, path in inputs:
        if not os.path.exists(path):
            sys.exit(f"Error: missing embedding file for '{label}': {path}")
        df = pd.read_csv(path, index_col=0)
        print(f"-> Loaded {label}: {df.shape[0]} genomes, {df.shape[1]} dimensions.")
        # Align by column position so differing headers ('dim_0' vs '0') don't matter
        df.columns = range(df.shape[1])
        frames.append(df)
        labels.append(label)

    widths = {f.shape[1] for f in frames}
    if len(widths) > 1:
        sys.exit(f"Error: embedding files have different numbers of dimensions: {sorted(widths)}")

    combined = pd.concat(frames, axis=0)
    datasets = pd.Series(np.concatenate([[l] * len(f) for l, f in zip(labels, frames)]),
                         index=combined.index)

    keep = combined.notna().all(axis=1).values
    if (~keep).any():
        print(f"[ALERT] Removed {(~keep).sum()} rows containing missing values.")
    return combined[keep], datasets[keep], labels


def cluster(X, args, seed):
    """2-D UMAP -> HDBSCAN. Returns (2-D coordinates, labels with -1 = noise)."""
    coords = umap.UMAP(n_neighbors=args.n_neighbors, min_dist=args.min_dist,
                       metric='euclidean', random_state=seed).fit_transform(X)
    mcs = max(2, round(args.min_cluster_fraction * len(X)))
    labels = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=args.min_samples).fit_predict(coords)
    return coords, labels


def n_clusters(labels):
    return len(set(labels)) - (1 if -1 in labels else 0)


def stability_report(X, args, primary_labels, datasets, dataset_order, minority):
    lines = [f"{len(X)} genomes"]
    fractions = datasets.value_counts(normalize=True)
    lines += [f"  {d}: {fractions[d]:.1%}" + (" (minority, drawn on top)" if d in minority else "")
              for d in dataset_order]
    lines.append(f"Main result (seed {args.seed}): {n_clusters(primary_labels)} clusters, "
                 f"{(primary_labels == -1).mean():.1%} unclustered")

    if args.stability_runs > 1:
        labelings = {args.seed: primary_labels}
        seed = 0
        while len(labelings) < args.stability_runs:
            seed += 1
            if seed not in labelings:
                print(f"   stability run with seed {seed}...")
                labelings[seed] = cluster(X, args, seed)[1]
        ari = [adjusted_rand_score(labelings[a], labelings[b])
               for a, b in itertools.combinations(labelings, 2)]
        counts = [n_clusters(l) for l in labelings.values()]
        mean_ari = np.mean(ari)
        if mean_ari >= 0.9:
            verdict = "Clusters were consistent across random runs."
        elif mean_ari >= 0.75:
            verdict = "Clusters were mostly consistent across random runs; boundaries may shift slightly."
        else:
            verdict = "Clusters varied between random runs; treat cluster assignments with caution."
        lines.append(f"Stability: {verdict}")
        lines.append(f"  {len(labelings)} runs gave {min(counts)}-{max(counts)} clusters; "
                     f"agreement (Adjusted Rand Index) mean {mean_ari:.2f}, lowest {np.min(ari):.2f}")
    return lines


def build_figure(df, minority, dataset_order, title, subtitle):
    # Hard-freeze boundaries with 5% layout margins
    x_min, x_max = df['UMAP_1'].min(), df['UMAP_1'].max()
    y_min, y_max = df['UMAP_2'].min(), df['UMAP_2'].max()
    x_pad, y_pad = (x_max - x_min) * 0.05, (y_max - y_min) * 0.05

    cluster_names = list(df['Cluster'].unique())  # already sorted numerically, noise first
    palette = px.colors.qualitative.Dark24
    color_map = {c: palette[i % len(palette)] for i, c in enumerate(c for c in cluster_names if c != '-1')}
    color_map['-1'] = NOISE_COLOR

    fig = px.scatter(
        df, x='UMAP_1', y='UMAP_2',
        color='Cluster', symbol='Dataset',
        symbol_sequence=SYMBOLS,
        category_orders={'Dataset': dataset_order},
        color_discrete_map=color_map,
        hover_data=['Genome_ID', 'Dataset'],
        title=f"{title}<br><sup>{subtitle}</sup>",
        opacity=0.65,
        template="plotly_white",
    )
    fig.update_traces(marker=dict(size=5, line=dict(width=0)), selector=dict(mode='markers'))

    def is_minority(trace):
        # px names each trace "<Cluster>, <Dataset>"
        return trace.name.split(", ")[-1] in minority

    # Minority datasets: slightly larger, more opaque, thin black outline
    for trace in fig.data:
        if is_minority(trace):
            trace.marker.size = 6
            trace.marker.opacity = 0.9
            trace.marker.line = dict(width=0.8, color='black')

    # Traces are drawn in order, so minority traces go last to sit on top
    fig.data = [t for t in fig.data if not is_minority(t)] + [t for t in fig.data if is_minority(t)]

    # Legend order is independent of drawing order: group entries by cluster, then dataset
    cluster_rank = {c: i for i, c in enumerate(cluster_names)}
    dataset_rank = {d: i for i, d in enumerate(dataset_order)}
    for trace in fig.data:
        cl, ds = trace.name.split(", ")
        trace.legendrank = 1000 + cluster_rank[cl] * len(dataset_order) + dataset_rank[ds]

    fig.update_layout(
        width=1200, height=850,
        legend_title_text='Cluster, Dataset (click to toggle)',
        legend=dict(itemsizing='constant'),
        xaxis=dict(title="UMAP 1", range=[x_min - x_pad, x_max + x_pad]),
        yaxis=dict(title="UMAP 2", range=[y_min - y_pad, y_max + y_pad]),
        uirevision='constant',
    )
    return fig


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("Step 1: loading embeddings")
    combined, datasets, labels = load_datasets(args.input)
    fractions = datasets.value_counts(normalize=True)
    minority = [d for d in labels if fractions[d] < args.minority_fraction]
    dataset_order = [d for d in labels if d not in minority] + minority

    print("Step 2: UMAP + HDBSCAN")
    X = combined.values
    coords, cluster_labels = cluster(X, args, args.seed)
    report = stability_report(X, args, cluster_labels, datasets, dataset_order, minority)

    results = pd.DataFrame(coords, columns=['UMAP_1', 'UMAP_2'], index=combined.index)
    results['Cluster'] = cluster_labels
    results['Dataset'] = datasets.values
    results.to_csv(os.path.join(args.out, "coordinates_and_clusters.csv"), index_label='Genome_ID')

    table = pd.crosstab(results['Cluster'], results['Dataset'])[labels]
    report += ["", "Genomes per cluster (-1 = unclustered):", table.to_string()]
    with open(os.path.join(args.out, "stability_report.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(report) + "\n")
    print("\n".join(report))

    print("\nStep 3: building interactive plot")
    df = results.reset_index()
    df.columns = ['Genome_ID', 'UMAP_1', 'UMAP_2', 'Cluster', 'Dataset']
    df = df.sort_values('Cluster')
    df['Cluster'] = df['Cluster'].astype(str)

    mcs = max(2, round(args.min_cluster_fraction * len(X)))
    subtitle = (f"{len(X)} genomes; UMAP n_neighbors={args.n_neighbors}, min_dist={args.min_dist}; "
                f"HDBSCAN min_cluster_size={mcs}, min_samples={args.min_samples}")
    fig = build_figure(df, minority, dataset_order, args.title, subtitle)

    html = pio.to_html(fig, full_html=True, include_plotlyjs=True)
    html = html.replace("<head>", f"<head><title>{args.title}</title>", 1)
    html_path = os.path.join(args.out, "interactive_umap.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Done. Open {html_path} in a browser.")


if __name__ == "__main__":
    main()
