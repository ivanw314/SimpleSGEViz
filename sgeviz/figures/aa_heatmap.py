import altair as alt
import pandas as pd


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


def make_plot(df: pd.DataFrame, gene: str = "") -> alt.Chart:
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

    title = f"Amino Acid Heatmap{' (' + gene + ')' if gene else ''} (n = {n})"

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

        chart = alt.vconcat(heatmap, vep_map, spacing=9).resolve_scale(color="independent")
    else:
        chart = heatmap

    return (
        chart
        .properties(title=alt.TitleParams(text=title, fontSize=22))
        .configure_axis(grid=False)
        .configure_view(stroke=None)
    )
