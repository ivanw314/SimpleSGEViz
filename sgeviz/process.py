from pathlib import Path

import numpy as np
import pandas as pd


def load_scores(files: dict):
    """Load and process scores into a single merged dataframe.

    SNVs and deletions are read from the same *allscores.tsv file and
    distinguished by ref allele length (1 = SNV, 4 = 3bp deletion).
    Thresholds are derived from the functional_consequence classifications
    already present in the scores file (GMM approach).

    Returns:
        df: combined scores dataframe with standardized consequence labels
        thresholds: [non-functional threshold, functional threshold]
    """
    snv_df = _read_snv_scores(files["all_scores"])
    del_df = _read_del_scores(files["all_scores"])

    thresholds = _get_gmm_thresholds(snv_df)

    df = pd.concat([snv_df, del_df], ignore_index=True)
    df = _rename_consequences(df)

    return df, thresholds


def _get_gmm_thresholds(snv_df: pd.DataFrame) -> list:
    """Derive thresholds from functional_consequence classifications in the scores data.

    Uses the highest-scoring functionally_abnormal variant as the lower threshold
    and the lowest-scoring functionally_normal variant as the upper threshold.
    """
    ab_df = snv_df.loc[snv_df["functional_consequence"] == "functionally_abnormal"]
    norm_df = snv_df.loc[snv_df["functional_consequence"] == "functionally_normal"]
    return [ab_df["score"].max(), norm_df["score"].min()]


def _read_snv_scores(path: Path) -> pd.DataFrame:
    """Read SNV rows from allscores file (ref allele length == 1)."""
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["ref"].str.len() == 1]
    df.loc[df["score"] >= 0, "functional_consequence"] = "functionally_normal"
    df["var_type"] = "snv"
    df["pos"] = df["pos"].astype(int)
    df["start"] = df["pos"]
    df["end"] = df["pos"]
    df["pos_id"] = df["pos"].astype(str) + ":" + df["alt"]
    df = df[df["variant_qc_flag"] != "WARN"]
    return df


def _read_del_scores(path: Path) -> pd.DataFrame:
    """Read 3bp deletion rows from allscores file (ref allele length == 4)."""
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["ref"].str.len() == 4]
    df["var_type"] = "3bp_del"
    df["start"] = df["pos"] + 1
    df["end"] = df["pos"] + 3
    df["pos_id"] = df["start"].astype(str) + "-" + df["end"].astype(str)
    df = df.astype({"pos": int, "start": int, "end": int})
    return df


def load_allele_freqs(files: dict, scores_df: pd.DataFrame):
    """Load gnomAD and/or Regeneron allele frequency files and merge with scores.

    Both CSV and Excel formats are supported. Files are detected by the optional
    'gnomad' and 'regeneron' keys in the files dict (set by find_genes).

    Returns a combined long-form dataframe with columns including 'Allele Frequency',
    'log_AF', and 'dataset', or None if neither file is present.
    """
    gnomad_path = files.get("gnomad")
    regeneron_path = files.get("regeneron")

    if gnomad_path is None and regeneron_path is None:
        return None

    # SNV scores only — deletions have no gnomAD/Regeneron entries
    snv_scores = scores_df.loc[scores_df["var_type"] == "snv"].copy()

    dfs = []

    if gnomad_path is not None:
        gnomad_df = _read_gnomad(gnomad_path)
        merged = pd.merge(snv_scores, gnomad_df, on="pos_id", how="inner")
        merged["dataset"] = "gnomAD"
        merged = merged.rename(columns={"gnomad_af": "Allele Frequency"})
        dfs.append(merged)
        print(f"  gnomAD: {len(merged)} variants with allele frequency data")

    if regeneron_path is not None:
        reg_df = _read_regeneron(regeneron_path)
        merged = pd.merge(snv_scores, reg_df, on="pos_id", how="inner")
        merged["dataset"] = "Regeneron"
        merged = merged.rename(columns={"regeneron_maf": "Allele Frequency"})
        dfs.append(merged)
        print(f"  Regeneron: {len(merged)} variants with allele frequency data")

    combined = pd.concat(dfs, ignore_index=True)
    combined["log_AF"] = np.log10(combined["Allele Frequency"])
    return combined


def _read_file(path: Path) -> pd.DataFrame:
    """Read CSV or Excel file based on extension."""
    if path.suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_excel(path)


def _read_gnomad(path: Path) -> pd.DataFrame:
    """Read gnomAD allele frequency file and return pos_id + gnomad_af columns.

    Supports gnomAD ID format: '2-214728604-T-A' -> pos_id '214728604:A'
    """
    df = _read_file(path)
    df = df[["gnomAD ID", "Allele Frequency"]].copy()
    df["pos_id"] = df["gnomAD ID"].transform(
        lambda x: x.split("-")[1] + ":" + x.split("-")[3]
    )
    df = df.rename(columns={"Allele Frequency": "gnomad_af"})
    return df[["pos_id", "gnomad_af"]]


def _read_regeneron(path: Path) -> pd.DataFrame:
    """Read Regeneron allele frequency file and return pos_id + regeneron_maf columns.

    Supports Variant format: '2:214728582:C:G' -> pos_id '214728582:G'
    """
    df = _read_file(path)
    df = df[["Variant", "AAF"]].copy()
    df = df.rename(columns={"AAF": "regeneron_maf", "Variant": "pos_id"})
    df["pos_id"] = df["pos_id"].transform(
        lambda x: x.split(":")[1] + ":" + x.split(":")[3]
    )
    return df[["pos_id", "regeneron_maf"]]


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
