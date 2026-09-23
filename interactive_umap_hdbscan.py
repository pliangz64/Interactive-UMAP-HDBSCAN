"""Interactive UMAP + HDBSCAN maps for genome datasets.

Two ways in:
  1. Embeddings  (--input LABEL CSV, repeatable): runs UMAP + HDBSCAN, then plots.
  2. Precomputed (--precomputed CSV): you already have 2-D coordinates (and optionally
     clusters); the tool only builds the interactive map.

Either way, --annotations CSV adds any per-genome columns (taxonomy, GC content, sample
site, ...) to the hover text, and --color-by colours the map by any column.

Examples:
    python interactive_umap_hdbscan.py --input Gut gut.csv --input Vaginal vaginal.csv --out results
    python interactive_umap_hdbscan.py --precomputed my_umap.csv --annotations meta.csv --color-by Genus
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

# Marker shapes per dataset: largest dataset first, then minority datasets
SYMBOLS = ['circle', 'diamond', 'square', 'triangle-up', 'cross', 'x', 'star', 'pentagon']
PALETTE = px.colors.qualitative.Dark24
NOISE_COLOR = '#c8c8c8'   # HDBSCAN cluster -1 (unclustered)
NA_COLOR = '#e0e0e0'      # missing annotation
OTHER_COLOR = '#8c8c8c'   # categories lumped together by --max-categories
HOVER_LIMIT = 15          # max annotation columns shown on hover
CONTINUOUS_MIN_UNIQUE = 13


def parse_args():
    p = argparse.ArgumentParser(
        description="Interactive UMAP + HDBSCAN maps for genome datasets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = p.add_argument_group("input (choose one)")
    src_x = src.add_mutually_exclusive_group(required=True)
    src_x.add_argument("--input", nargs=2, action="append", metavar=("LABEL", "CSV"),
                       help="Dataset label and embedding CSV (first column = genome ID, remaining "
                            "columns = numeric embedding). Repeat for each dataset.")
    src_x.add_argument("--precomputed", metavar="CSV",
                       help="Table with genome ID in the first column plus 2-D coordinates, and "
                            "optionally cluster and dataset columns. Any other columns are treated "
                            "as annotations.")

    p.add_argument("--annotations", metavar="CSV",
                   help="Optional per-genome annotations (first column = genome ID). Matched to "
                        "genomes by ID; every column is shown on hover and can be used with --color-by.")
    p.add_argument("--color-by", metavar="COLUMN",
                   help="Column to colour points by. Default: Cluster if present, otherwise Dataset.")
    p.add_argument("--max-categories", type=int, default=20,
                   help="When colouring by a categorical column, the most common categories get their "
                        "own colour; the rest are grouped as 'Other'.")
    p.add_argument("--out", default="results", help="Output directory.")
    p.add_argument("--title", default="Interactive UMAP + HDBSCAN", help="Plot and page title.")
    p.add_argument("--minority-fraction", type=float, default=0.40,
                   help="Datasets below this fraction of all genomes are drawn on top with an outline.")

    pre = p.add_argument_group("precomputed column names")
    pre.add_argument("--x-col", default="UMAP_1", help="X coordinate column.")
    pre.add_argument("--y-col", default="UMAP_2", help="Y coordinate column.")
    pre.add_argument("--cluster-col", default="Cluster", help="Cluster column (optional in the file).")
    pre.add_argument("--dataset-col", default="Dataset", help="Dataset column (optional in the file).")

    emb = p.add_argument_group("UMAP + HDBSCAN settings (embedding input only)")
    emb.add_argument("--min-cluster-fraction", type=float, default=0.05,
                     help="HDBSCAN min_cluster_size as a fraction of all genomes. Lower (e.g. 0.017) "
                          "for finer sub-clusters.")
    emb.add_argument("--min-samples", type=int, default=5,
                     help="HDBSCAN min_samples. Values >= 15 can collapse 2-D maps into 1-2 clusters.")
    emb.add_argument("--metric", default="euclidean",
                     help="UMAP distance metric for the embeddings, e.g. 'euclidean' or 'cosine' for "
                          "dense embeddings, 'jaccard' for binary presence/absence matrices. Any "
                          "metric supported by umap-learn works.")
    emb.add_argument("--n-neighbors", type=int, default=30, help="UMAP n_neighbors.")
    emb.add_argument("--min-dist", type=float, default=0.1, help="UMAP min_dist.")
    emb.add_argument("--seed", type=int, default=42, help="Random seed for the main result.")
    emb.add_argument("--stability-runs", type=int, default=5,
                     help="Number of UMAP seeds used to check cluster stability (0 or 1 to skip).")
    return p.parse_args()


# -----------------------------------------------------------------------------
# Loading
# -----------------------------------------------------------------------------

def read_table(path, what):
    if not os.path.exists(path):
        sys.exit(f"Error: {what} file not found: {path}")
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    df.index.name = 'Genome_ID'
    if df.index.duplicated().any():
        print(f"[ALERT] {what}: {df.index.duplicated().sum()} duplicate genome IDs; keeping the first of each.")
        df = df[~df.index.duplicated()]
    return df


def load_embeddings(inputs):
    frames, labels = [], []
    for label, path in inputs:
        df = read_table(path, f"embedding '{label}'")
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


def load_precomputed(args):
    df = read_table(args.precomputed, "precomputed")
    missing = [c for c in (args.x_col, args.y_col) if c not in df.columns]
    if missing:
        sys.exit(f"Error: precomputed file is missing coordinate column(s) {missing}. "
                 f"Available columns: {list(df.columns)}. Use --x-col / --y-col.")
    df = df.rename(columns={args.x_col: 'UMAP_1', args.y_col: 'UMAP_2',
                            args.cluster_col: 'Cluster', args.dataset_col: 'Dataset'})
    if 'Dataset' not in df.columns:
        df['Dataset'] = 'All genomes'
    df['Dataset'] = df['Dataset'].astype(str)
    df = df.dropna(subset=['UMAP_1', 'UMAP_2'])
    print(f"-> Loaded {len(df)} genomes with precomputed coordinates"
          + (" and clusters." if 'Cluster' in df.columns else " (no cluster column)."))
    labels = list(dict.fromkeys(df['Dataset']))  # order of first appearance
    return df, labels


def add_annotations(df, path):
    ann = read_table(path, "annotations")
    clash = [c for c in ann.columns if c in df.columns]
    if clash:
        print(f"[ALERT] Annotation columns already present in the data were skipped: {clash}")
        ann = ann.drop(columns=clash)
    matched = df.index.isin(ann.index).sum()
    print(f"-> Annotations: {ann.shape[1]} columns; matched {matched} of {len(df)} genomes by ID.")
    if matched == 0:
        print("[ALERT] No genome IDs matched. Check that the first column of the annotation file "
              "uses the same IDs as the main data.")
    return df.join(ann, how='left')


# -----------------------------------------------------------------------------
# UMAP + HDBSCAN
# -----------------------------------------------------------------------------

def cluster(X, args, seed):
    """2-D UMAP -> HDBSCAN. Returns (2-D coordinates, labels with -1 = unclustered)."""
    import hdbscan
    import umap
    coords = umap.UMAP(n_neighbors=args.n_neighbors, min_dist=args.min_dist,
                       metric=args.metric, random_state=seed).fit_transform(X)
    labels = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size(args, len(X)),
                             min_samples=args.min_samples).fit_predict(coords)
    return coords, labels


def min_cluster_size(args, n):
    return max(2, round(args.min_cluster_fraction * n))


def n_clusters(labels):
    return len(set(labels)) - (1 if -1 in labels else 0)


def stability_lines(X, args, primary_labels):
    from sklearn.metrics import adjusted_rand_score
    lines = [f"Main result (seed {args.seed}): {n_clusters(primary_labels)} clusters, "
             f"{(primary_labels == -1).mean():.1%} unclustered"]
    if args.stability_runs <= 1:
        return lines

    labelings = {args.seed: primary_labels}
    seed = 0
    while len(labelings) < args.stability_runs:
        seed += 1
        if seed not in labelings:
            print(f"   stability run with seed {seed}...")
            labelings[seed] = cluster(X, args, seed)[1]
    # Compare only genomes clustered in both runs: counting shared noise (-1) as a
    # "cluster" would inflate agreement
    ari, shared = [], []
    for a, b in itertools.combinations(labelings, 2):
        both = (labelings[a] >= 0) & (labelings[b] >= 0)
        shared.append(both.mean())
        ari.append(adjusted_rand_score(labelings[a][both], labelings[b][both]) if both.sum() > 1 else 0.0)
    counts = [n_clusters(l) for l in labelings.values()]
    mean_ari = np.mean(ari)
    # Rule-of-thumb cutoffs chosen for this tool, not an established standard
    if mean_ari >= 0.9:
        verdict = "Clusters were consistent across random runs."
    elif mean_ari >= 0.75:
        verdict = "Clusters were mostly consistent across random runs; boundaries may shift slightly."
    else:
        verdict = "Clusters varied between random runs; treat cluster assignments with caution."
    lines.append(f"Stability: {verdict}")
    lines.append(f"  {len(labelings)} runs gave {min(counts)}-{max(counts)} clusters; "
                 f"agreement (Adjusted Rand Index, genomes clustered in both runs) "
                 f"mean {mean_ari:.2f}, lowest {np.min(ari):.2f}")
    lines.append(f"  On average {np.mean(shared):.0%} of genomes were clustered in both runs of a pair.")
    lines.append("  Verdict cutoffs (0.9 / 0.75) are rule-of-thumb guides chosen for this tool.")
    return lines


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def to_label(v):
    if pd.isna(v):
        return "NA"
    if isinstance(v, (float, np.floating)) and float(v).is_integer():
        return str(int(v))
    return str(v)


def to_hover(v):
    if pd.isna(v):
        return "NA"
    if isinstance(v, (float, np.floating)):
        return f"{v:.4g}"
    return str(v)


def categorical_levels(series, max_categories):
    """String labels per genome and the category order (numeric order if numeric-like, else by size)."""
    labels = series.map(to_label)
    counts = labels[labels != "NA"].value_counts()
    try:
        order = sorted(counts.index, key=float)
    except ValueError:
        order = list(counts.index)  # most common first

    if len(order) > max_categories:
        keep = set(counts.index[:max_categories - 1])
        labels = labels.where(labels.isin(keep) | (labels == "NA"), "Other")
        order = [c for c in order if c in keep] + ["Other"]
        print(f"[NOTE] {len(counts)} categories in colour column; showing the {max_categories - 1} most "
              f"common and grouping the rest as 'Other' (see --max-categories).")
    if (labels == "NA").any():
        order.append("NA")
    return labels, order


def build_figure(df, color_col, dataset_order, minority, args, title, subtitle):
    n_ds = len(dataset_order)
    symbol = {d: SYMBOLS[i % len(SYMBOLS)] for i, d in enumerate(dataset_order)}

    # Hover: genome ID, dataset, cluster, then annotation columns
    extra = [c for c in df.columns if c not in ('UMAP_1', 'UMAP_2', 'Dataset', 'Cluster')]
    if len(extra) > HOVER_LIMIT:
        print(f"[NOTE] Showing the first {HOVER_LIMIT} of {len(extra)} annotation columns on hover.")
    hover_cols = ['Dataset'] + (['Cluster'] if 'Cluster' in df.columns else []) + extra[:HOVER_LIMIT]
    hovertemplate = ("<b>%{customdata[0]}</b><br>"
                     + "<br>".join(f"{c}: %{{customdata[{i + 1}]}}" for i, c in enumerate(hover_cols))
                     + "<extra></extra>")

    def points(sub, is_minority, **marker):
        customdata = np.column_stack([sub.index.astype(str)] + [sub[c].map(to_hover) for c in hover_cols])
        style = dict(size=6, opacity=0.9, line=dict(width=0.8, color='black')) if is_minority \
            else dict(size=5, opacity=0.65, line=dict(width=0))
        return go.Scattergl(
            x=sub['UMAP_1'], y=sub['UMAP_2'], mode='markers',
            marker=dict(**style, **marker),
            customdata=customdata, hovertemplate=hovertemplate,
        )

    values = df[color_col]
    continuous = (color_col not in ('Cluster', 'Dataset')
                  and pd.api.types.is_numeric_dtype(values)
                  and values.nunique() >= CONTINUOUS_MIN_UNIQUE)

    traces = []  # (is_minority, trace); drawing order is majority first, minority on top
    if continuous:
        for rank, ds in enumerate(dataset_order):
            sub = df[df['Dataset'] == ds]
            has = sub[color_col].notna()
            t = points(sub[has], ds in minority, symbol=symbol[ds], color=sub.loc[has, color_col],
                       coloraxis='coloraxis')
            t.update(name=ds, legendrank=1000 + rank)
            traces.append((ds in minority, t))
            if (~has).any():
                t = points(sub[~has], ds in minority, symbol=symbol[ds], color=NA_COLOR)
                t.update(name=f"NA, {ds}" if n_ds > 1 else "NA", legendrank=2000 + rank)
                traces.append((ds in minority, t))
        legend_title = "Dataset (click to toggle)"
    else:
        labels, order = categorical_levels(values, args.max_categories)
        colors, i = {}, 0
        for cat in order:
            if cat == "-1" and color_col == 'Cluster':
                colors[cat] = NOISE_COLOR
            elif cat == "NA":
                colors[cat] = NA_COLOR
            elif cat == "Other":
                colors[cat] = OTHER_COLOR
            else:
                colors[cat] = PALETTE[i % len(PALETTE)]
                i += 1
        for ci, cat in enumerate(order):
            for di, ds in enumerate(dataset_order):
                sub = df[(labels == cat) & (df['Dataset'] == ds)]
                if sub.empty:
                    continue
                name = cat if (n_ds == 1 or color_col == 'Dataset') else f"{cat}, {ds}"
                t = points(sub, ds in minority, symbol=symbol[ds], color=colors[cat])
                # Legend grouped by category, then dataset (independent of drawing order)
                t.update(name=name, legendrank=1000 + ci * n_ds + di)
                traces.append((ds in minority, t))
        legend_title = (f"{color_col}, Dataset" if n_ds > 1 and color_col != 'Dataset' else color_col) \
            + " (click to toggle)"

    fig = go.Figure([t for m, t in traces if not m] + [t for m, t in traces if m])

    # Hard-freeze boundaries with 5% layout margins
    x_min, x_max = df['UMAP_1'].min(), df['UMAP_1'].max()
    y_min, y_max = df['UMAP_2'].min(), df['UMAP_2'].max()
    x_pad, y_pad = (x_max - x_min) * 0.05, (y_max - y_min) * 0.05
    fig.update_layout(
        title=f"{title}<br><sup>{subtitle}</sup>",
        template="plotly_white",
        width=1200, height=850,
        legend=dict(title_text=legend_title, itemsizing='constant'),
        xaxis=dict(title="UMAP 1", range=[x_min - x_pad, x_max + x_pad]),
        yaxis=dict(title="UMAP 2", range=[y_min - y_pad, y_max + y_pad]),
        uirevision='constant',
    )
    if continuous:
        # Horizontal colour bar under the x-axis title, clear of the legend on the right
        fig.update_layout(
            margin=dict(b=180),
            coloraxis=dict(
                colorscale='Viridis',
                colorbar=dict(title=dict(text=color_col, side='top'), orientation='h',
                              y=-0.18, yanchor='top', len=0.5, thickness=12),
            ),
        )
    return fig


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    report = []

    print("Step 1: loading data")
    if args.input:
        combined, datasets, labels = load_embeddings(args.input)
        print("Step 2: UMAP + HDBSCAN")
        X = combined.values
        coords, cluster_labels = cluster(X, args, args.seed)
        df = pd.DataFrame(coords, columns=['UMAP_1', 'UMAP_2'], index=combined.index)
        df['Cluster'] = cluster_labels
        df['Dataset'] = datasets.values
        report += stability_lines(X, args, cluster_labels)
        subtitle = (f"{len(df)} genomes; UMAP metric={args.metric}, n_neighbors={args.n_neighbors}, "
                    f"min_dist={args.min_dist}; "
                    f"HDBSCAN min_cluster_size={min_cluster_size(args, len(df))}, "
                    f"min_samples={args.min_samples}")
    else:
        df, labels = load_precomputed(args)
        subtitle = f"{len(df)} genomes; precomputed coordinates"

    if args.annotations:
        df = add_annotations(df, args.annotations)

    color_col = args.color_by or ('Cluster' if 'Cluster' in df.columns else 'Dataset')
    if color_col not in df.columns:
        sys.exit(f"Error: --color-by column '{color_col}' not found. Available: "
                 f"{[c for c in df.columns if c not in ('UMAP_1', 'UMAP_2')]}")
    if color_col != 'Cluster':
        subtitle += f"; coloured by {color_col}"

    fractions = df['Dataset'].value_counts(normalize=True)
    minority = [d for d in labels if fractions[d] < args.minority_fraction] if len(labels) > 1 else []
    dataset_order = [d for d in labels if d not in minority] + minority

    # Summary report
    summary = [f"{len(df)} genomes"]
    summary += [f"  {d}: {fractions[d]:.1%}" + (" (minority, drawn on top)" if d in minority else "")
                for d in dataset_order]
    report = summary + report
    if 'Cluster' in df.columns:
        table = pd.crosstab(df['Cluster'], df['Dataset'])[labels]
        report += ["", "Genomes per cluster (-1 = unclustered):", table.to_string()]
    with open(os.path.join(args.out, "report.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(report) + "\n")
    print("\n".join(report))

    if args.input:
        # Includes any annotations, so it can be reused later with --precomputed
        df.to_csv(os.path.join(args.out, "coordinates_and_clusters.csv"), index_label='Genome_ID')

    print("\nStep 3: building interactive plot")
    fig = build_figure(df, color_col, dataset_order, minority, args, args.title, subtitle)
    html = pio.to_html(fig, full_html=True, include_plotlyjs=True)
    html = html.replace("<head>", f"<head><title>{args.title}</title>", 1)
    html_path = os.path.join(args.out, "interactive_umap.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Done. Open {html_path} in a browser.")


if __name__ == "__main__":
    main()
