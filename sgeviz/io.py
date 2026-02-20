from pathlib import Path

import altair as alt
import pandas as pd


def find_genes(input_dir: Path) -> dict:
    """Discover all gene datasets in the input directory.

    Detects genes by finding all *allscores.tsv files and extracting the gene
    name from each filename (e.g. 20260129_RAD51Dallscores.tsv -> RAD51D).
    Companion files (*modelparams.tsv, *snvcounts.tsv, *delcounts.tsv) are
    then located using the gene name as a search key.

    Returns a dict mapping gene name -> files dict, e.g.:
        {"RAD51D": {"all_scores": Path(...), "snv_counts": Path(...), ...}}
    """
    allscores_files = sorted(input_dir.glob("*allscores.tsv"))
    if not allscores_files:
        raise FileNotFoundError(f"No '*allscores.tsv' files found in {input_dir}")

    def find_one(pattern):
        matches = list(input_dir.glob(pattern))
        if not matches:
            raise FileNotFoundError(f"Could not find '{pattern}' in {input_dir}")
        if len(matches) > 1:
            raise ValueError(
                f"Multiple files match '{pattern}': "
                + ", ".join(str(m) for m in matches)
            )
        return matches[0]

    def find_optional(pattern):
        matches = list(input_dir.glob(pattern))
        return matches[0] if len(matches) == 1 else None

    genes = {}
    for allscores_path in allscores_files:
        # e.g. '20260129_RAD51Dallscores' -> 'RAD51D'
        gene = allscores_path.stem.split("_")[-1].replace("allscores", "")
        genes[gene] = {
            "all_scores": allscores_path,
            "model_params": find_one(f"*{gene}modelparams.tsv"),
            "snv_counts": find_one(f"*{gene}snvcounts.tsv"),
            "del_counts": find_one(f"*{gene}delcounts.tsv"),
            # Optional allele frequency files (CSV or Excel)
            "gnomad": find_optional(f"*{gene}*gnomAD*"),
            "regeneron": find_optional(f"*{gene}*Regeneron*"),
        }

    return genes


def load_counts(files: dict) -> pd.DataFrame:
    """Load and merge SNV and deletion counts files, standardizing column names."""
    snv_df = pd.read_csv(files["snv_counts"], sep="\t")
    del_df = pd.read_csv(files["del_counts"], sep="\t")

    del_df["var_type"] = "3bp_del"
    snv_df["start"] = snv_df["pos"]
    snv_df["end"] = snv_df["pos"]
    snv_df["var_type"] = "snv"
    snv_df = snv_df.drop("pos", axis=1)

    df = pd.concat([snv_df, del_df], ignore_index=True)

    # Strip gene prefix from target names (e.g. BARD1_X1B -> 1B, RAD51D_X1 -> 1)
    df["target"] = df["target"].str.replace(r"^[A-Z0-9]+_X", "", regex=True)

    df = df.rename(
        columns={
            "D05_R1": "D05 R1",
            "D05_R2": "D05 R2",
            "D05_R3": "D05 R3",
            "D13_R1": "D13 R1",
            "D13_R2": "D13 R2",
            "D13_R3": "D13 R3",
            "D17_R1": "D17 R1",
            "D17_R2": "D17 R2",
            "D17_R3": "D17 R3",
        }
    )
    return df


def save_figure(chart: alt.Chart, path: Path):
    """Save an Altair chart. Format is inferred from the file extension.

    HTML is fully self-contained and interactive.
    PNG and SVG require vl-convert-python to be installed.
    """
    chart.save(str(path))
    print(f"  Saved: {path.name}")
