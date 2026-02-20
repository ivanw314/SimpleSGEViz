from pathlib import Path

import altair as alt
import pandas as pd


def find_files(input_dir: Path) -> dict:
    """Discover required input TSV files in the given directory by suffix pattern."""
    patterns = {
        "all_scores": "*allscores.tsv",
        "model_params": "*modelparams.tsv",
        "snv_counts": "*snvcounts.tsv",
        "del_counts": "*delcounts.tsv",
    }

    found = {}
    for key, pattern in patterns.items():
        matches = list(input_dir.glob(pattern))
        if not matches:
            raise FileNotFoundError(f"Could not find '{pattern}' in {input_dir}")
        if len(matches) > 1:
            raise ValueError(
                f"Multiple files match '{pattern}' in {input_dir}: "
                + ", ".join(str(m) for m in matches)
            )
        found[key] = matches[0]

    return found


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
