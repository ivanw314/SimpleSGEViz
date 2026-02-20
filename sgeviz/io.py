from pathlib import Path

import altair as alt
import pandas as pd


def find_files(input_dir: Path) -> dict:
    """Discover required input TSV files in the given directory by suffix pattern."""
    patterns = {
        "snv_scores": "*snvscores.tsv",
        "del_scores": "*delscores.tsv",
        "model_params": "*modelparams.tsv",
        "counts": "*counts.tsv",
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


def load_counts(counts_path: Path) -> pd.DataFrame:
    """Load SNV counts file and standardize column names."""
    df = pd.read_csv(counts_path, sep="\t")
    df["pos_id"] = df["pos"].astype(str) + ":" + df["alt"]
    df["target"] = df["target"].transform(lambda x: x[7:])  # strip 'BARD1_X' prefix
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
