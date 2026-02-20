from pathlib import Path

import pandas as pd


def load_scores(files: dict):
    """Load and process SNV + deletion scores into a single merged dataframe.

    Returns:
        df: combined scores dataframe with standardized consequence labels
        thresholds: [non-functional threshold, functional threshold]
    """
    thresholds = _read_thresholds(files["model_params"])
    snv_df = _read_snv_scores(files["snv_scores"])
    del_df = _read_del_scores(files["del_scores"])

    df = pd.concat([snv_df, del_df], ignore_index=True)
    df = _rename_consequences(df)

    return df, thresholds


def _read_thresholds(path: Path) -> list:
    df = pd.read_csv(path, sep="\t")
    return [df["thresh_abnormal"][0], df["thresh_normal"][0]]


def _read_snv_scores(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.loc[df["score"] >= 0, "functional_consequence"] = "functionally_normal"
    df["var_type"] = "snv"
    df["pos"] = df["pos"].astype(int)
    df["start"] = df["pos"]
    df["end"] = df["pos"]
    df["pos_id"] = df["pos"].astype(str) + ":" + df["alt"]
    df = df[df["variant_qc_flag"] != "WARN"]
    return df


def _read_del_scores(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["var_type"] = "3bp_del"
    df["start"] = df["pos"] + 1
    df["end"] = df["pos"] + 3
    df["pos_id"] = df["start"].astype(str) + "-" + df["end"].astype(str)
    df = df.astype({"pos": int, "start": int, "end": int})
    return df


def _rename_consequences(df: pd.DataFrame) -> pd.DataFrame:
    """Rename raw VEP consequence terms to display labels."""
    df = df.rename(columns={"consequence": "Consequence"})

    df.loc[df["Consequence"].str.contains("missense", na=False), "Consequence"] = "Missense"
    df.loc[df["Consequence"] == "synonymous_variant", "Consequence"] = "Synonymous"
    df.loc[df["Consequence"] == "intron_variant", "Consequence"] = "Intron"
    df.loc[df["Consequence"] == "stop_gained", "Consequence"] = "Stop Gained"
    df.loc[df["Consequence"] == "stop_lost", "Consequence"] = "Stop Lost"
    df.loc[df["Consequence"].str.contains("site", na=False), "Consequence"] = "Canonical Splice"
    df.loc[df["Consequence"].str.contains("ing_var", na=False), "Consequence"] = "Splice Region"
    df.loc[df["Consequence"].str.contains("UTR", na=False), "Consequence"] = "UTR Variant"
    df.loc[df["Consequence"] == "start_lost", "Consequence"] = "Start Lost"
    df.loc[df["var_type"] == "3bp_del", "Consequence"] = "3bp Deletion"

    return df
