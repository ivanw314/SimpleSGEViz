from pathlib import Path

import altair as alt
import pandas as pd

_DOMAIN_PALETTE = [
    "#B9DBF4",  # light blue
    "#C8DBC8",  # light green
    "#FF9A00",  # orange
    "#ffc976",  # amber
    "#D5A6BD",  # pink
    "#A4C2F4",  # periwinkle
    "#FFD966",  # gold
    "#B6D7A8",  # sage green
]

_DEL_TYPES = ["Inframe Indel", "Stop Gained"]
_DEL_PALETTE = ["black", "#ffc007"]
_DEL_CONSEQUENCE_MAP = {
    # process.py already maps inframe_indel → "Inframe Indel" and stop_gained → "Stop Gained";
    # merge the rarer cases into "Inframe Indel" for the del panel display.
    "Start Lost": "Inframe Indel",
    "Stop Lost": "Inframe Indel",
}

_AA_ORDER = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
    "*", "Mis. Min.", "Mis. Mean",
]

_VEP_COLS = {
    "am_score": "AlphaMissense",
    "revel_score": "REVEL",
    "MutPred2": "MutPred2",
}


def _load_domains(path, prot_length: int):
    """Load a domain annotation file and return (tier0_df, tier1_df).

    Expected columns: region_name, aa_residues (format "start-end").
    Optional columns:
      color (hex)  – auto-assigned from _DOMAIN_PALETTE if absent.
      tier (int)   – 0 = main strip (default), 1 = sub-feature strip.

    Each returned DataFrame has start, end, label, color columns.
    tier1_df is None when no tier-1 rows exist.
    """
    suffix = Path(path).suffix.lower()
    raw = pd.read_excel(path) if suffix in (".xlsx", ".xls") else pd.read_csv(path)

    raw[["start", "end"]] = raw["aa_residues"].str.split("-", expand=True).astype(int)
    raw = raw.sort_values("start").reset_index(drop=True)

    if "tier" not in raw.columns:
        raw["tier"] = 0

    if "color" not in raw.columns:
        palette = _DOMAIN_PALETTE * ((len(raw) // len(_DOMAIN_PALETTE)) + 1)
        raw["color"] = palette[: len(raw)]

    def _fill_gaps(rows, gap_color):
        segments = []
        pos = 1
        for _, row in rows.iterrows():
            if pos < row["start"]:
                segments.append({"start": pos, "end": row["start"], "label": "", "color": gap_color})
            segments.append({
                "start": row["start"],
                "end": row["end"],
                "label": row["region_name"],
                "color": row["color"],
            })
            pos = row["end"]
        if pos <= prot_length:
            segments.append({"start": pos, "end": prot_length + 1, "label": "", "color": gap_color})
        return pd.DataFrame(segments)

    tier0 = raw.loc[raw["tier"] == 0].sort_values("start").reset_index(drop=True)
    tier1 = raw.loc[raw["tier"] == 1].sort_values("start").reset_index(drop=True)

    tier0_df = _fill_gaps(tier0, "#CCCCCC")
    tier1_df = _fill_gaps(tier1, "#FFFFFF") if len(tier1) > 0 else None

    return tier0_df, tier1_df


def _make_domain_cartoon(tier0_df: pd.DataFrame, tier1_df, prot_length: int, width: int) -> alt.Chart:
    """Colored rect strip(s) with domain name labels centered in each segment.

    tier0_df is the main strip; tier1_df (optional) is a thinner sub-feature strip
    rendered below it (e.g. Walker motifs inside a RecA domain).
    """
    x_scale = alt.Scale(domain=[1, prot_length + 1])

    def _strip(df, height, font_size):
        df = df.copy()
        df["center"] = (df["start"] + df["end"]) / 2
        label_df = df.loc[df["label"] != ""].copy()
        label_df["y_mid"] = 0.5  # maps to vertical center in a [0,1] scale

        rects = (
            alt.Chart(df)
            .mark_rect(stroke="black", strokeWidth=0.5)
            .encode(
                x=alt.X("start:Q", axis=None, scale=x_scale),
                x2="end:Q",
                y=alt.value(0),
                y2=alt.value(height),
                color=alt.Color("color:N", scale=None, legend=None),
                tooltip=[
                    alt.Tooltip("label:N", title="Domain"),
                    alt.Tooltip("start:Q", title="Start"),
                    alt.Tooltip("end:Q", title="End"),
                ],
            )
            .properties(width=width, height=height)
        )

        labels = (
            alt.Chart(label_df)
            .mark_text(fontSize=font_size, fontWeight="bold", baseline="middle", align="center")
            .encode(
                x=alt.X("center:Q", axis=None, scale=x_scale),
                y=alt.Y("y_mid:Q", scale=alt.Scale(domain=[0, 1]), axis=None),
                text="label:N",
            )
            .properties(width=width, height=height)
        )

        return alt.layer(rects, labels).resolve_scale(y="independent")

    strip0 = _strip(tier0_df, 25, 13)
    if tier1_df is None:
        return strip0

    strip1 = _strip(tier1_df, 15, 10)
    return alt.vconcat(strip0, strip1, spacing=1)


def _make_del_panel(
    del_df: pd.DataFrame, thresholds, prot_length: int, width: int
) -> alt.Chart:
    """Scatter plot of 3bp deletion fitness scores by protein position."""
    del_df = del_df.copy()
    if "amino_acid_change" in del_df.columns:
        del_df = del_df.loc[~del_df["amino_acid_change"].isin(["---"])]

    del_df["_cds_start"] = del_df["CDS_position"].str.split("-").str[0]
    del_df["_cds_end"] = del_df["CDS_position"].str.split("-").str[1]
    del_df = del_df.loc[
        ~del_df["_cds_start"].isin(["?"]) & ~del_df["_cds_end"].isin(["?"])
    ]
    del_df["_cds_start"] = del_df["_cds_start"].astype(int)
    del_df["ps_aa_start"] = ((del_df["_cds_start"] + 2) / 3).round(2)
    del_df["Consequence"] = del_df["Consequence"].replace(_DEL_CONSEQUENCE_MAP)

    scatter = (
        alt.Chart(del_df)
        .mark_point(strokeWidth=3, size=75, opacity=1)
        .encode(
            x=alt.X(
                "ps_aa_start:Q",
                title="",
                axis=alt.Axis(ticks=False, labels=False),
                scale=alt.Scale(domain=[1, prot_length + 1], nice=False),
            ),
            y=alt.Y(
                "score:Q",
                title="Fitness Score",
                axis=alt.Axis(titleFontSize=20, labelFontSize=18),
                scale=alt.Scale(domain=[-0.5, 0.1]),
            ),
            color=alt.Color(
                "Consequence:N",
                scale=alt.Scale(domain=_DEL_TYPES, range=_DEL_PALETTE),
                legend=alt.Legend(
                    title="Consequence", labelFontSize=18, titleFontSize=20, orient="right"
                ),
            ),
            shape=alt.Shape(
                "Consequence:N",
                scale=alt.Scale(domain=_DEL_TYPES, range=["square", "triangle"]),
                legend=None,
            ),
            order=alt.Order("Consequence:N", sort="ascending"),
        )
        .properties(width=width, height=150)
    )

    layers = [scatter]
    if thresholds is not None:
        layers.append(
            alt.Chart(del_df)
            .mark_rule(color="red", strokeDash=[8, 8], strokeWidth=1)
            .encode(y=alt.datum(thresholds[0]))
        )
        layers.append(
            alt.Chart(del_df)
            .mark_rule(color="blue", strokeDash=[8, 8], strokeWidth=1)
            .encode(y=alt.datum(thresholds[1]))
        )
    return alt.layer(*layers)


def make_plot(df: pd.DataFrame, gene: str = "", thresholds=None, domains_path=None) -> alt.Chart:
    """Generate amino acid substitution heatmap of SGE fitness scores.

    X-axis: amino acid position (per-residue bins).
    Y-axis: amino acid substitution, sorted by 1-letter code then *, Mis. Min., Mis. Mean.
    Color: fitness score (bluepurple reversed, clamped to [-0.2, 0]).

    If am_score, revel_score, or MutPred2 columns are present in df, a VEP predictor
    sub-panel is appended below the main heatmap with an independent color scale.

    Args:
        df: scores dataframe from process.load_scores (must contain amino_acid_change).
        gene: Gene name used in the plot title.
    """
    alt.data_transformers.disable_max_rows()

    snv_df = df.loc[df["var_type"] == "snv"].copy()
    snv_df = snv_df.loc[~snv_df["amino_acid_change"].isin(["-", "---"])]
    snv_df = snv_df.dropna(subset=["amino_acid_change"])

    # Parse amino acid change string: e.g. 'A123G' -> og_AA='A', AA_change='G', AApos=123
    snv_df["og_AA"] = snv_df["amino_acid_change"].str[0]
    snv_df["AA_change"] = snv_df["amino_acid_change"].str[-1]
    snv_df["AApos"] = pd.to_numeric(
        snv_df["amino_acid_change"].str[1:-1], errors="coerce"
    )
    snv_df = snv_df.dropna(subset=["AApos"])
    snv_df["AApos"] = snv_df["AApos"].astype(int)

    if "max_SpliceAI" in snv_df.columns:
        snv_df = snv_df.loc[snv_df["max_SpliceAI"] <= 0.2]

    prot_length = snv_df["AApos"].max()
    n = len(snv_df)
    n_del = int((df["var_type"] == "3bp_del").sum()) if "var_type" in df.columns else 0
    width = 4 * prot_length
    height_per_row = 25

    # Missense min/mean rows (exclude stop gained)
    mis_df = snv_df.loc[snv_df["Consequence"] != "Stop Gained"]
    min_df = mis_df.groupby("AApos")["score"].min().reset_index()
    min_df["AA_change"] = "Mis. Min."
    mean_df = mis_df.groupby("AApos")["score"].mean().reset_index()
    mean_df["AA_change"] = "Mis. Mean"

    plot_df = pd.concat(
        [snv_df[["AApos", "AA_change", "score"]], min_df, mean_df],
        ignore_index=True,
    )

    gene_label = f" ({gene})" if gene else ""
    count_label = f"n = {n_del} deletions, {n} SNVs" if n_del > 0 else f"n = {n}"
    title = f"Deletion and Heatmap{gene_label} ({count_label})"

    vep_cols_present = {k: v for k, v in _VEP_COLS.items() if k in snv_df.columns}
    has_vep = bool(vep_cols_present)

    heatmap = alt.Chart(plot_df).mark_rect().encode(
        x=alt.X(
            "AApos:Q",
            title="" if has_vep else "Amino Acid Position",
            axis=alt.Axis(
                labels=not has_vep,
                ticks=not has_vep,
                values=[] if has_vep else list(range(0, prot_length, 25)),
                titleFontSize=22,
                labelFontSize=18,
            ),
            scale=alt.Scale(domain=[1, prot_length + 1]),
            bin=alt.Bin(maxbins=prot_length + 1, minstep=1),
        ),
        y=alt.Y(
            "AA_change:N",
            title="Amino Acid Substitution",
            sort=_AA_ORDER,
            axis=alt.Axis(labelFontSize=18, titleFontSize=20),
        ),
        color=alt.Color(
            "score:Q",
            title="Fitness Score",
            scale=alt.Scale(domain=[-0.2, 0], clamp=True, scheme="bluepurple", reverse=True),
            legend=alt.Legend(titleFontSize=18, labelFontSize=16),
        ),
    ).properties(width=width, height=height_per_row * len(_AA_ORDER))

    if has_vep:
        vep_summary = (
            snv_df.groupby("AApos")[list(vep_cols_present.keys())]
            .mean()
            .reset_index()
            .rename(columns=vep_cols_present)
        )
        vep_df = vep_summary.melt(id_vars=["AApos"], var_name="Predictor", value_name="score")
        vep_order = list(vep_cols_present.values())

        vep_map = alt.Chart(vep_df).mark_rect().encode(
            x=alt.X(
                "AApos:Q",
                title="Amino Acid Position",
                axis=alt.Axis(
                    values=list(range(0, prot_length, 25)),
                    titleFontSize=22,
                    labelFontSize=18,
                ),
                scale=alt.Scale(domain=[1, prot_length + 1]),
                bin=alt.Bin(maxbins=prot_length + 1, minstep=1),
            ),
            y=alt.Y(
                "Predictor:N",
                title="",
                sort=vep_order,
                axis=alt.Axis(labelFontSize=16),
            ),
            color=alt.Color(
                "score:Q",
                title="Predictor Score",
                scale=alt.Scale(domain=[0, 1], clamp=True, scheme="bluepurple"),
                legend=alt.Legend(
                    titleFontSize=18,
                    labelFontSize=16,
                    gradientLength=height_per_row * len(vep_cols_present),
                    values=[0, 1.0],
                ),
            ),
        ).properties(width=width, height=height_per_row * len(vep_cols_present))

        lower = alt.vconcat(heatmap, vep_map, spacing=9).resolve_scale(color="independent")
    else:
        lower = heatmap

    # Build panel stack top-to-bottom: domains → deletions → heatmap(+vep)
    panels = []

    if domains_path is not None:
        tier0_df, tier1_df = _load_domains(domains_path, prot_length)
        panels.append(_make_domain_cartoon(tier0_df, tier1_df, prot_length, width))

    has_del = (
        "var_type" in df.columns
        and "CDS_position" in df.columns
        and (df["var_type"] == "3bp_del").any()
    )
    if has_del:
        del_df = df.loc[df["var_type"] == "3bp_del"].copy()
        panels.append(_make_del_panel(del_df, thresholds, prot_length, width))

    panels.append(lower)

    if len(panels) == 1:
        chart = panels[0]
    else:
        chart = alt.vconcat(*panels, spacing=9).resolve_scale(
            x="shared", color="independent", shape="independent"
        )

    return (
        chart
        .properties(title=alt.TitleParams(text=title, fontSize=22, anchor="middle"))
        .configure_axis(grid=False)
        .configure_view(stroke=None)
    )
