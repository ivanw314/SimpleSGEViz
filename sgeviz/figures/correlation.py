import altair as alt
import pandas as pd
from natsort import natsorted


_COMBOS = [
    ("D05 R1", "D05 R2"),
    ("D05 R1", "D05 R3"),
    ("D05 R2", "D05 R3"),
    ("D13 R1", "D13 R2"),
    ("D13 R1", "D13 R3"),
    ("D13 R2", "D13 R3"),
]


def compute_correlations(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Compute pairwise Pearson r between all replicate combinations per SGE target.

    Returns a long-form dataframe with columns: Targets, Tests, r_correlation.
    """
    records = []
    rep_cols = ["D05 R1", "D05 R2", "D05 R3", "D13 R1", "D13 R2", "D13 R3"]

    for target, group in counts_df.groupby("target"):
        group_reps = group[rep_cols]
        for col1, col2 in _COMBOS:
            r = round(group_reps[col1].corr(group_reps[col2]), 3)
            records.append(
                {"Targets": target, "Tests": f"{col1} vs {col2}", "r_correlation": r}
            )

    return pd.DataFrame(records)


def make_heatmap(df: pd.DataFrame) -> alt.Chart:
    """Generate Pearson r heatmap across SGE targets and replicate comparisons."""
    targets = natsorted(df["Targets"].unique().tolist())

    base = alt.Chart(
        df,
        title=alt.TitleParams(text="Correlation of Replicates", fontSize=32),
    ).encode(
        x=alt.X(
            "Tests:N",
            axis=alt.Axis(
                title="",
                titleFontSize=28,
                labelFontSize=24,
                labelLimit=300,
                labelAngle=45,
            ),
        ),
        y=alt.Y(
            "Targets:N",
            sort=targets,
            axis=alt.Axis(
                title="SGE Target Region", titleFontSize=28, labelFontSize=24
            ),
        ),
    )

    rect = base.mark_rect().encode(
        color=alt.Color(
            "r_correlation:Q",
            scale=alt.Scale(domain=[0.2, 1]),
            legend=alt.Legend(
                title="Pearson's r", titleFontSize=24, labelFontSize=22
            ),
        ),
        tooltip=[alt.Tooltip("r_correlation", title="Pearson's r: ")],
    ).properties(width=600, height=900)

    text_color = (
        alt.when(alt.datum.r_correlation > 0.5)
        .then(alt.value("white"))
        .otherwise(alt.value("black"))
    )

    text = base.mark_text(baseline="middle", fontSize=20).encode(
        text=alt.Text("r_correlation:Q", format="0.2f"),
        color=text_color,
    ).transform_filter("isValid(datum.r_correlation)")

    return rect + text
