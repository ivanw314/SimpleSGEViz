"""SGE Visualization Pipeline

Usage:
    python pipeline.py <input_dir> <output_dir> [--format html|png|svg]

The input directory must contain:
    *allscores.tsv      SNV + deletion fitness scores (combined file)
    *modelparams.tsv    SGE model thresholds
    *snvcounts.tsv      Per-replicate SNV counts
    *delcounts.tsv      Per-replicate deletion counts

Optional (allele frequency figures generated only if detected):
    *{gene}*gnomAD*     gnomAD allele frequencies (CSV or Excel)
    *{gene}*Regeneron*  Regeneron allele frequencies (CSV or Excel)

Outputs (saved to output_dir):
    {gene}_histogram_stripplot    Score distribution histogram + strip plot
    {gene}_correlation_heatmap    Replicate Pearson r heatmap
    {gene}_scores_across_gene     Per-exon scatter plot of scores vs. genomic position
    {gene}_maf_vs_score           Allele frequency vs. score heatmap (if AF files present)

PNG and SVG output require vl-convert-python (pip install vl-convert-python).
"""

import argparse
import sys
from pathlib import Path

from sgeviz import io, process
from sgeviz.figures import correlation, histogram_strip, maf_score, scores_gene


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate standard SGE visualization figures from raw score files."
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing input TSV files",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory to write output figures",
    )
    parser.add_argument(
        "--format",
        choices=["html", "png", "svg"],
        default="png",
        help="Output format for figures (default: html). "
             "PNG/SVG require vl-convert-python.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.input_dir.is_dir():
        sys.exit(f"Error: input directory not found: {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fmt = args.format

    # --- Discover genes ---
    print(f"Scanning for gene datasets in: {args.input_dir}")
    genes = io.find_genes(args.input_dir)
    print(f"  Found {len(genes)} gene(s): {', '.join(genes)}")

    # --- Process each gene ---
    for gene, files in genes.items():
        print(f"\n[{gene}] Loading data...")
        scores_df, thresholds = process.load_scores(files)
        counts_df = io.load_counts(files)
        print(f"  {len(scores_df)} variants loaded")

        print(f"[{gene}] Generating figures (format: {fmt})")

        hist, strip = histogram_strip.make_figures(scores_df, thresholds, gene=gene)
        io.save_figure(
            histogram_strip.combine(hist, strip),
            args.output_dir / f"{gene}_histogram_stripplot.{fmt}",
        )

        r_df = correlation.compute_correlations(counts_df)
        io.save_figure(
            correlation.make_heatmap(r_df, gene=gene),
            args.output_dir / f"{gene}_correlation_heatmap.{fmt}",
        )

        io.save_figure(
            scores_gene.make_plot(scores_df, thresholds, gene=gene),
            args.output_dir / f"{gene}_scores_across_gene.{fmt}",
        )

        maf_df = process.load_allele_freqs(files, scores_df)
        if maf_df is not None:
            io.save_figure(
                maf_score.make_plot(maf_df, gene=gene),
                args.output_dir / f"{gene}_maf_vs_score.{fmt}",
            )
        else:
            print(f"[{gene}] No allele frequency files found, skipping MAF figure.")

    print(f"\nDone. Figures saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
