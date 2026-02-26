"""SGE Visualization Pipeline

Usage:
    python pipeline.py <input_dir> <output_dir> [--format html|png|svg] [--excel]

The input directory must contain:
    *allscores.tsv      SNV + deletion fitness scores (combined file)
    *modelparams.tsv    SGE model thresholds
    *snvcounts.tsv      Per-replicate SNV counts
    *delcounts.tsv      Per-replicate deletion counts

Optional (figures generated only if detected):
    *{gene}*gnomAD*     gnomAD allele frequencies (CSV or Excel)
    *{gene}*Regeneron*  Regeneron allele frequencies (CSV or Excel)
    *{gene}*editrates*  Library edit rates (TSV with target_rep + edit_rate columns)
    *{gene}*cartoon*    Gene cartoon Excel file (sheets: exon_coords, metadata, optionally lib_coords)
    *{gene}*vep*        VEP Excel output (.xlsx) with AlphaMissense, REVEL, CADD, SpliceAI scores

Outputs (saved to output_dir):
    {gene}_histogram_stripplot    Score distribution histogram + strip plot
    {gene}_correlation_heatmap    Replicate Pearson r heatmap
    {gene}_scores_across_gene     Per-exon scatter plot of scores vs. genomic position
    {gene}_aa_heatmap             Amino acid substitution heatmap (if amino_acid_change column present)
    {gene}_predictor_scatter      Predictor vs. fitness score scatter panels (if predictor columns present)
    {gene}_clinvar_strip          ClinVar classification strip plot (if *{gene}*ClinVar*SNV* file present)
    {gene}_clinvar_roc            ROC curve for SGE score B/LB vs P/LP classification (if ClinVar file present)
    {gene}_maf_vs_score           Allele frequency vs. score heatmap (if AF files present)
    {gene}_edit_rate_barplot      Library edit rate bar plot by target (if *{gene}*editrates* file present)
    {gene}_exon_cartoon           Exon structure cartoon (if *{gene}*cartoon* file present, no lib_coords sheet)
    {gene}_library_cartoon        Exon + library design cartoon (if *{gene}*cartoon* file with lib_coords sheet)
    {gene}_data.xlsx              Multi-sheet Excel workbook (if --excel flag is set)

PNG and SVG output require vl-convert-python (pip install vl-convert-python).
Excel output requires openpyxl (pip install openpyxl).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from sgeviz import io, process
from sgeviz.figures import aa_heatmap, clinvar_strip, correlation, edit_rate_barplot, gene_cartoon, histogram_strip, maf_score, predictor_scatter, scores_gene


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
    parser.add_argument(
        "--excel",
        action="store_true",
        default=False,
        help="Also write a multi-sheet Excel workbook ({gene}_annotated_data.xlsx) "
             "containing scores, thresholds, counts, and any optional datasets. "
             "Requires openpyxl.",
    )
    parser.add_argument(
        "--protein-length",
        type=int,
        default=None,
        metavar="N",
        help="Known full protein length (aa). If the data covers fewer residues, "
             "the x-axis will be extended to this length. If omitted, you will be "
             "prompted interactively for each gene.",
    )
    parser.add_argument(
        "--px-per-aa",
        type=int,
        default=4,
        metavar="N",
        help="Pixels per amino acid column in the AA heatmap (default: 4). "
             "Reduce to produce a narrower figure.",
    )
    parser.add_argument(
        "--gene-name",
        type=str,
        default=None,
        metavar="NAME",
        help="Override the gene name used in figure titles and output filenames. "
             "Cannot be used when multiple gene datasets are detected.",
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

    if args.gene_name is not None:
        if len(genes) > 1:
            sys.exit(
                f"Error: --gene-name cannot be used when multiple gene datasets are detected "
                f"({', '.join(genes)})."
            )
        original = next(iter(genes))
        genes = {args.gene_name: genes[original]}
        print(f"  Gene name overridden: '{original}' -> '{args.gene_name}'")

    # --- Resolve protein lengths ---
    # If --protein-length was not supplied, prompt interactively for each gene.
    protein_lengths: dict[str, int | None] = {}
    for gene in genes:
        if args.protein_length is not None:
            protein_lengths[gene] = args.protein_length
        else:
            raw = input(f"Protein length for {gene} (press Enter to estimate from data): ").strip()
            protein_lengths[gene] = int(raw) if raw else None

    # --- Process each gene ---
    for gene, files in genes.items():
        print(f"\n[{gene}] Loading data...")
        scores_df, thresholds = process.load_scores(files)
        scores_df = process.load_vep(files, scores_df)
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

        if "amino_acid_change" in scores_df.columns:
            domains_path = files.get("domains")
            if domains_path is not None:
                print(f"[{gene}] Domain file detected: {domains_path.name}")
            else:
                print(f"[{gene}] No domain file detected (place a *{gene}*domain* file in input dir to enable).")
            io.save_figure(
                aa_heatmap.make_plot(
                    scores_df, gene=gene, thresholds=thresholds,
                    domains_path=domains_path,
                    protein_length=protein_lengths[gene],
                    px_per_aa=args.px_per_aa,
                ),
                args.output_dir / f"{gene}_aa_heatmap.{fmt}",
            )
        else:
            print(f"[{gene}] No amino_acid_change column, skipping AA heatmap.")

        pred_plot = predictor_scatter.make_plot(scores_df, thresholds, gene=gene)
        if pred_plot is not None:
            io.save_figure(
                pred_plot,
                args.output_dir / f"{gene}_predictor_scatter.{fmt}",
            )
        else:
            print(f"[{gene}] No predictor score columns found, skipping predictor scatter.")

        clinvar_df = process.load_clinvar(files, scores_df)
        if clinvar_df is not None:
            io.save_figure(
                clinvar_strip.make_strip(clinvar_df, thresholds, gene=gene),
                args.output_dir / f"{gene}_clinvar_strip.{fmt}",
            )
            roc = clinvar_strip.make_roc(clinvar_df, gene=gene)
            if roc is not None:
                io.save_figure(
                    roc,
                    args.output_dir / f"{gene}_clinvar_roc.{fmt}",
                )
            else:
                print(f"[{gene}] Insufficient B/LB and P/LP variants for ROC curve.")
        else:
            print(f"[{gene}] No ClinVar file found, skipping ClinVar figures.")

        maf_df = process.load_allele_freqs(files, scores_df)
        if maf_df is not None:
            io.save_figure(
                maf_score.make_plot(maf_df, gene=gene),
                args.output_dir / f"{gene}_maf_vs_score.{fmt}",
            )
        else:
            print(f"[{gene}] No allele frequency files found, skipping MAF figure.")

        cartoon_data = io.load_cartoon(files)
        if cartoon_data is not None:
            exon_df, lib_df, meta_df = cartoon_data
            if lib_df is not None and not lib_df.empty:
                cartoon_chart = gene_cartoon.make_library_cartoon(exon_df, lib_df, meta_df)
                cartoon_name = f"{gene}_library_cartoon"
            else:
                cartoon_chart = gene_cartoon.make_exon_cartoon(exon_df, meta_df)
                cartoon_name = f"{gene}_exon_cartoon"
            io.save_figure(
                cartoon_chart,
                args.output_dir / f"{cartoon_name}.{fmt}",
            )
        else:
            print(f"[{gene}] No cartoon file found, skipping gene cartoon.")

        edit_rates_df = io.load_edit_rates(files)
        if edit_rates_df is not None:
            io.save_figure(
                edit_rate_barplot.make_plot(edit_rates_df, gene=gene),
                args.output_dir / f"{gene}_edit_rate_barplot.{fmt}",
            )
        else:
            print(f"[{gene}] No edit rates file found, skipping edit rate bar plot.")

        if args.excel:
            print(f"[{gene}] Writing Excel workbook...")
            thresh_df = pd.DataFrame({
                "non_functional_threshold": [thresholds[0]],
                "functional_threshold": [thresholds[1]],
            })

            # Build merged scores sheet: left-join optional annotations by pos_id
            merged_scores = scores_df.copy()
            if clinvar_df is not None:
                merged_scores = pd.merge(
                    merged_scores,
                    clinvar_df[["pos_id", "Germline classification"]],
                    on="pos_id",
                    how="left",
                )
            if maf_df is not None:
                af_wide = (
                    maf_df[["pos_id", "dataset", "Allele Frequency"]]
                    .pivot_table(
                        index="pos_id",
                        columns="dataset",
                        values="Allele Frequency",
                        aggfunc="first",
                    )
                    .reset_index()
                )
                af_wide.columns.name = None
                af_wide = af_wide.rename(columns={
                    "gnomAD": "gnomad_af",
                    "Regeneron": "regeneron_af",
                })
                merged_scores = pd.merge(merged_scores, af_wide, on="pos_id", how="left")

            sheets = {"scores": merged_scores, "thresholds": thresh_df, "counts": counts_df}
            if edit_rates_df is not None:
                sheets["edit_rates"] = edit_rates_df
            io.save_excel(sheets, args.output_dir / f"{gene}_data.xlsx")

    print(f"\nDone. Figures saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
