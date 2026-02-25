import fnmatch
from pathlib import Path

import altair as alt
import pandas as pd


def find_genes(input_dir: Path) -> dict:
    """Discover all gene datasets in the input directory.

    Detects genes by finding all *allscores.tsv files and extracting the gene
    name from each filename. Both dot-separated (e.g. CTCF.allscores.tsv) and
    run-together (e.g. 20260129_RAD51Dallscores.tsv) naming conventions are
    supported. Companion files (*modelparams.tsv, *snvcounts.tsv, *delcounts.tsv)
    are then located using the gene name as a search key.

    Returns a dict mapping gene name -> files dict, e.g.:
        {"RAD51D": {"all_scores": Path(...), "snv_counts": Path(...), ...}}
    """
    allscores_files = sorted(input_dir.glob("*allscores.tsv"))
    if not allscores_files:
        raise FileNotFoundError(f"No '*allscores.tsv' files found in {input_dir}")

    def find_one(*patterns):
        for pattern in patterns:
            matches = list(input_dir.glob(pattern))
            if len(matches) == 1:
                return matches[0]
            if len(matches) > 1:
                raise ValueError(
                    f"Multiple files match '{pattern}': "
                    + ", ".join(str(m) for m in matches)
                )
        raise FileNotFoundError(
            f"Could not find any of {patterns} in {input_dir}"
        )

    def find_optional(pattern):
        matches = list(input_dir.glob(pattern))
        return matches[0] if len(matches) == 1 else None

    def find_optional_icase(pattern, exclude=None):
        matches = [
            p for p in input_dir.iterdir()
            if not p.name.startswith("~$")
            and fnmatch.fnmatch(p.name.lower(), pattern.lower())
            and (exclude is None or not fnmatch.fnmatch(p.name.lower(), exclude.lower()))
        ]
        return matches[0] if len(matches) == 1 else None

    genes = {}
    for allscores_path in allscores_files:
        # Handles both GENE.allscores.tsv and GENEallscores.tsv (with optional prefix)
        stem_part = allscores_path.stem.split("_")[-1]
        gene = stem_part.removesuffix(".allscores").removesuffix("allscores")
        genes[gene] = {
            "all_scores": allscores_path,
            "model_params": find_one(f"*{gene}modelparams.tsv", f"*{gene}.modelparams.tsv"),
            "snv_counts": find_one(f"*{gene}snvcounts.tsv", f"*{gene}.snvcounts.tsv"),
            "del_counts": find_one(f"*{gene}delcounts.tsv", f"*{gene}.delcounts.tsv"),
            # Optional allele frequency files (CSV or Excel)
            "gnomad": find_optional(f"*{gene}*gnomAD*"),
            "regeneron": find_optional(f"*{gene}*Regeneron*"),
            # Optional ClinVar SNV file (tab-delimited .txt from ClinVar download)
            "clinvar": find_optional_icase(f"*{gene}*clinvar*snv*"),
            # Optional domain annotation file (Excel or CSV with region_name + aa_residues cols)
            # Exclude files that also contain "cartoon" — those are gene-cartoon files, not domain files.
            "domains": find_optional_icase(f"*{gene}*domain*", exclude=f"*cartoon*"),
            # Optional library edit rates file (*editrates*.tsv)
            "edit_rates": find_optional_icase(f"*{gene}*editrates*"),
            # Optional gene cartoon file (Excel with exon_coords, metadata, and optionally lib_coords)
            "cartoon": find_optional_icase(f"*{gene}*cartoon*"),
            # Optional VEP output file (Excel .xlsx or tab-delimited .txt) for predictor scores
            "vep": find_optional_icase(f"*{gene}*vep*"),
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


def load_cartoon(files: dict):
    """Load a gene cartoon Excel file if present.

    Expected sheets: ``exon_coords`` (required), ``metadata`` (required),
    ``lib_coords`` (optional).

    Returns ``(exon_df, lib_df, metadata_df)`` where ``lib_df`` is ``None``
    when the ``lib_coords`` sheet is absent.  Returns ``None`` if no cartoon
    file was detected.
    """
    path = files.get("cartoon")
    if path is None:
        return None
    xl = pd.ExcelFile(path)
    exon_df = xl.parse("exon_coords")
    lib_df = xl.parse("lib_coords") if "lib_coords" in xl.sheet_names else None
    meta_df = xl.parse("metadata")
    return exon_df, lib_df, meta_df


def load_edit_rates(files: dict):
    """Load an edit rates TSV file if present.

    Expected columns: target_rep, edit_rate.
    Returns the raw DataFrame or None if no edit rates file was detected.
    """
    path = files.get("edit_rates")
    if path is None:
        return None
    return pd.read_csv(path, sep="\t")


def save_figure(chart: alt.Chart, path: Path):
    """Save an Altair chart. Format is inferred from the file extension.

    HTML is fully self-contained and interactive.
    PNG and SVG require vl-convert-python to be installed.
    """
    chart.save(str(path))
    print(f"  Saved: {path.name}")


def save_excel(sheets: dict, path: Path):
    """Write a multi-sheet Excel workbook. Requires openpyxl.

    Args:
        sheets: Dict mapping sheet name -> DataFrame (insertion order preserved).
        path: Output path (.xlsx).
    """
    with pd.ExcelWriter(str(path), engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"  Saved: {path.name}")
