# app.py
from __future__ import annotations

import numpy as np
import pandas as pd

from dash import Dash, html, dcc, Input, Output
from dash import dash_table
import plotly.express as px

from scipy import stats


CELL_TYPES = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def fdr_bh(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini Hochberg FDR adjustment.
    Returns q values aligned to pvals, NaN where pvals are NaN.
    """
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)

    mask = ~np.isnan(p)
    pv = p[mask]
    if pv.size == 0:
        return q

    order = np.argsort(pv)
    ranked = pv[order]
    m = ranked.size

    q_ranked = ranked * m / (np.arange(1, m + 1))
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_ranked = np.clip(q_ranked, 0, 1)

    q_out = np.empty_like(ranked)
    q_out[order] = q_ranked
    q[mask] = q_out
    return q


def find_column_case_insensitive(df: pd.DataFrame, name: str) -> str | None:
    target = name.strip().lower()
    for c in df.columns:
        if str(c).strip().lower() == target:
            return c
    return None


def ensure_canonical_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Creates canonical columns used by the app if they exist under alternate names.
    It does not delete original columns.

    Canonical names used in this app:
      sample, subject, project, disease, sample_type, tissue, treatment, response, sex, time_from_treatment_start
    """
    out = df.copy()

    # sample
    if "sample" not in out.columns:
        c = find_column_case_insensitive(out, "sample")
        if c is not None:
            out.rename(columns={c: "sample"}, inplace=True)

    # subject: accept subject or subject_id
    if "subject" not in out.columns:
        c = find_column_case_insensitive(out, "subject")
        if c is not None:
            out.rename(columns={c: "subject"}, inplace=True)
        else:
            c2 = find_column_case_insensitive(out, "subject_id")
            if c2 is not None:
                out.rename(columns={c2: "subject"}, inplace=True)

    # other metadata
    for canon in ["project", "disease", "treatment", "response", "sex", "time_from_treatment_start", "sample_type", "tissue"]:
        if canon not in out.columns:
            c = find_column_case_insensitive(out, canon)
            if c is not None:
                out.rename(columns={c: canon}, inplace=True)

    # strip common string cols
    for col in ["project", "disease", "treatment", "response", "sex", "sample_type", "tissue"]:
        if col in out.columns:
            out[col] = out[col].astype(str).str.strip()

    return out


def build_summary_table(cell_count_path: str) -> pd.DataFrame:
    """
    Reads a wide cell count table and returns the long summary table with:
    sample, total_count, population, count, percentage
    """
    cc = pd.read_csv(cell_count_path)
    cc = ensure_canonical_columns(cc)

    missing = [c for c in (["sample"] + CELL_TYPES) if c not in cc.columns]
    if missing:
        raise ValueError(f"cell-count.csv missing columns: {missing}")

    cc["total_count"] = cc[CELL_TYPES].sum(axis=1)

    long_df = cc.melt(
        id_vars=["sample", "total_count"],
        value_vars=CELL_TYPES,
        var_name="population",
        value_name="count",
    )

    long_df["percentage"] = np.where(
        long_df["total_count"] > 0,
        100.0 * long_df["count"] / long_df["total_count"],
        np.nan,
    )

    long_df["population"] = pd.Categorical(long_df["population"], categories=CELL_TYPES, ordered=True)
    long_df = long_df.sort_values(["sample", "population"]).reset_index(drop=True)
    return long_df


def load_and_join_data(cell_count_path: str, metadata_path: str) -> pd.DataFrame:
    summary = build_summary_table(cell_count_path)

    meta = pd.read_csv(metadata_path)
    meta = ensure_canonical_columns(meta)

    if "sample" not in meta.columns:
        raise ValueError("sample-metadata.csv must contain a sample column (any case variation is ok).")

    df = summary.merge(meta, on="sample", how="left")
    df = ensure_canonical_columns(df)
    return df


def detect_pbmc_column(df: pd.DataFrame) -> str | None:
    if "tissue" in df.columns:
        return "tissue"
    if "sample_type" in df.columns:
        return "sample_type"
    return None


def filter_core(df: pd.DataFrame, disease: str, treatment: str, pbmc_only: bool) -> pd.DataFrame:
    out = df.copy()

    if disease != "All" and "disease" in out.columns:
        out = out[out["disease"].str.contains(disease, case=False, na=False)]

    if treatment != "All" and "treatment" in out.columns:
        out = out[out["treatment"].str.lower() == treatment.lower()]

    if pbmc_only:
        pbmc_col = detect_pbmc_column(out)
        if pbmc_col is not None:
            out = out[out[pbmc_col].str.contains("pbmc", case=False, na=False)]

    return out


def subject_level_means(df: pd.DataFrame) -> pd.DataFrame:
    """
    Produces one row per subject per population.
    Mean counts and mean percentages are averages across that subject's samples.
    Response and sex are taken as the first non null value per subject.
    """
    if "subject" not in df.columns:
        raise ValueError("Missing metadata column: subject")

    needed = ["subject", "population", "count", "percentage"]
    for c in needed:
        if c not in df.columns:
            raise ValueError(f"Missing column: {c}")

    tmp = df.copy()
    tmp["count"] = pd.to_numeric(tmp["count"], errors="coerce")
    tmp["percentage"] = pd.to_numeric(tmp["percentage"], errors="coerce")

    # subject level fields if present
    subject_fields = []
    for c in ["response", "sex"]:
        if c in tmp.columns:
            subject_fields.append(c)

    # aggregate mean per subject per population
    agg = (
        tmp.groupby(["subject", "population"], as_index=False)
        .agg(mean_count=("count", "mean"), mean_percentage=("percentage", "mean"))
    )

    # attach subject level response and sex
    if subject_fields:
        subj = (
            tmp[["subject"] + subject_fields]
            .replace("nan", np.nan)
            .dropna(subset=["subject"])
            .drop_duplicates(subset=["subject"])
        )
        agg = agg.merge(subj, on="subject", how="left")

    agg["population"] = pd.Categorical(agg["population"], categories=CELL_TYPES, ordered=True)
    agg = agg.sort_values(["population", "subject"]).reset_index(drop=True)
    return agg


def welch_stats_by_population_subject_mean_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Welch t test per population comparing subject level mean_count between response yes and no.
    """
    subdf = subject_level_means(df)

    if "response" not in subdf.columns:
        raise ValueError("Missing metadata column: response")

    resp = subdf["response"].astype(str).str.strip().str.lower()

    rows = []
    for pop in CELL_TYPES:
        d = subdf[subdf["population"] == pop].copy()
        r = resp.loc[d.index]

        a = d[r == "yes"]["mean_count"].to_numpy(dtype=float)
        b = d[r == "no"]["mean_count"].to_numpy(dtype=float)

        a = a[~np.isnan(a)]
        b = b[~np.isnan(b)]

        if len(a) < 2 or len(b) < 2:
            t_stat = np.nan
            pval = np.nan
        else:
            t_stat, pval = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")

        rows.append(
            {
                "population": pop,
                "n_subjects_responders": int(len(a)),
                "n_subjects_nonresponders": int(len(b)),
                "mean_count_responders": float(np.mean(a)) if len(a) else np.nan,
                "mean_count_nonresponders": float(np.mean(b)) if len(b) else np.nan,
                "t_stat": float(t_stat) if t_stat == t_stat else np.nan,
                "p_value": float(pval) if pval == pval else np.nan,
            }
        )

    res = pd.DataFrame(rows)
    res["q_value"] = fdr_bh(res["p_value"].to_numpy(dtype=float))
    res["significant_fdr_0_05"] = res["q_value"] < 0.05

    return res.sort_values(["significant_fdr_0_05", "q_value"], ascending=[False, True]).reset_index(drop=True)


def make_datatable(df: pd.DataFrame, page_size: int = 15) -> dash_table.DataTable:
    return dash_table.DataTable(
        data=df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in df.columns],
        page_size=page_size,
        sort_action="native",
        filter_action="native",
        style_table={"overflowX": "auto"},
        style_cell={"fontFamily": "Arial", "fontSize": 13, "padding": "6px"},
        style_header={"fontWeight": "bold"},
    )


def pick_default(options: list[str], desired_lower: str) -> str:
    for x in options:
        if str(x).strip().lower() == desired_lower:
            return x
    return "All"


# Load once at startup
DF = load_and_join_data("data/cell-count.csv", "data/sample-metadata.csv")

disease_options = ["All"]
if "disease" in DF.columns:
    disease_options += sorted({d for d in DF["disease"].dropna().unique() if str(d).strip().lower() != "nan"})

treatment_options = ["All"]
if "treatment" in DF.columns:
    treatment_options += sorted({t for t in DF["treatment"].dropna().unique() if str(t).strip().lower() != "nan"})


app = Dash(__name__)
app.title = "Cell Frequency Dashboard"

app.layout = html.Div(
    style={"maxWidth": "1200px", "margin": "0 auto", "padding": "16px"},
    children=[
        html.H2("Miraclib PBMC Cell Population Dashboard"),
        html.Div(
            style={"display": "flex", "gap": "16px", "flexWrap": "wrap"},
            children=[
                html.Div(
                    children=[
                        html.Label("Disease"),
                        dcc.Dropdown(
                            id="disease_dd",
                            options=[{"label": x, "value": x} for x in disease_options],
                            value=pick_default(disease_options, "melanoma"),
                            clearable=False,
                        ),
                    ],
                    style={"minWidth": "260px"},
                ),
                html.Div(
                    children=[
                        html.Label("Treatment"),
                        dcc.Dropdown(
                            id="treatment_dd",
                            options=[{"label": x, "value": x} for x in treatment_options],
                            value=pick_default(treatment_options, "miraclib"),
                            clearable=False,
                        ),
                    ],
                    style={"minWidth": "260px"},
                ),
                html.Div(
                    children=[
                        html.Label("PBMC only"),
                        dcc.RadioItems(
                            id="pbmc_radio",
                            options=[{"label": "Yes", "value": "yes"}, {"label": "No", "value": "no"}],
                            value="yes",
                            inline=True,
                        ),
                    ],
                    style={"minWidth": "220px"},
                ),
            ],
        ),

        html.Hr(),

        html.H3("Part 2  Data Overview"),
        html.Div(id="overview_table_container"),

        html.Hr(),

        html.H3("Part 3  Responders vs Nonresponders"),
        dcc.Graph(id="boxplot_graph"),
        html.H4("Welch t test on subject level mean counts"),
        html.Div(id="stats_table_container"),

        html.Hr(),

        html.H3("Part 4  Baseline subset analysis"),
        html.Div(
            style={"display": "flex", "gap": "16px", "flexWrap": "wrap"},
            children=[
                html.Div(
                    children=[html.H4("Samples per project"), html.Div(id="baseline_project_table")],
                    style={"flex": "1", "minWidth": "320px"},
                ),
                html.Div(
                    children=[html.H4("Subjects by response"), html.Div(id="baseline_response_table")],
                    style={"flex": "1", "minWidth": "320px"},
                ),
                html.Div(
                    children=[html.H4("Subjects by sex"), html.Div(id="baseline_sex_table")],
                    style={"flex": "1", "minWidth": "320px"},
                ),
            ],
        ),
    ],
)


@app.callback(
    Output("overview_table_container", "children"),
    Output("boxplot_graph", "figure"),
    Output("stats_table_container", "children"),
    Output("baseline_project_table", "children"),
    Output("baseline_response_table", "children"),
    Output("baseline_sex_table", "children"),
    Input("disease_dd", "value"),
    Input("treatment_dd", "value"),
    Input("pbmc_radio", "value"),
)
def refresh(disease: str, treatment: str, pbmc_yesno: str):
    pbmc_only = pbmc_yesno == "yes"
    filtered = filter_core(DF, disease=disease, treatment=treatment, pbmc_only=pbmc_only)

    # Part 2 overview table
    overview_cols = ["sample", "total_count", "population", "count", "percentage"]
    overview = filtered[overview_cols].copy()
    overview_table = make_datatable(overview, page_size=12)

    # Part 3 subject level boxplot and Welch stats
    try:
        subdf = subject_level_means(filtered)

        if "response" in subdf.columns:
            fig = px.box(
                subdf,
                x="population",
                y="mean_count",
                color="response",
                points="all",
                title="Subject level mean cell counts by population, responders vs nonresponders",
            )
            fig.update_layout(xaxis_title="Population", yaxis_title="Mean cell count per subject")
        else:
            fig = px.box(title="Missing metadata column: response")

        stats_df = welch_stats_by_population_subject_mean_counts(filtered)
        show_cols = [
            "population",
            "n_subjects_responders",
            "n_subjects_nonresponders",
            "mean_count_responders",
            "mean_count_nonresponders",
            "t_stat",
            "p_value",
            "q_value",
            "significant_fdr_0_05",
        ]
        stats_table = make_datatable(stats_df[show_cols], page_size=8)

    except Exception as e:
        fig = px.box(title="Unable to compute Part 3")
        stats_table = html.Div(str(e))

    # Part 4 baseline subset analysis
    base = filtered.copy()

    if "time_from_treatment_start" in base.columns:
        t = pd.to_numeric(base["time_from_treatment_start"], errors="coerce")
        base = base[t == 0]

    # Samples per project
    if "project" in base.columns:
        sp = (
            base[["project", "sample"]]
            .drop_duplicates()
            .groupby("project")["sample"]
            .nunique()
            .reset_index()
            .rename(columns={"sample": "n_samples"})
            .sort_values("n_samples", ascending=False)
        )
        project_table = make_datatable(sp, page_size=10)
    else:
        project_table = html.Div("Missing metadata column: project")

    # Subjects by response
    if "subject" in base.columns and "response" in base.columns:
        sr = (
            base[["subject", "response"]]
            .drop_duplicates()
            .groupby("response")["subject"]
            .nunique()
            .reset_index()
            .rename(columns={"subject": "n_subjects"})
            .sort_values("n_subjects", ascending=False)
        )
        response_table = make_datatable(sr, page_size=10)
    else:
        missing_bits = []
        if "subject" not in base.columns:
            missing_bits.append("subject")
        if "response" not in base.columns:
            missing_bits.append("response")
        response_table = html.Div("Missing metadata column: " + " or ".join(missing_bits))

    # Subjects by sex
    if "subject" in base.columns and "sex" in base.columns:
        ss = (
            base[["subject", "sex"]]
            .drop_duplicates()
            .groupby("sex")["subject"]
            .nunique()
            .reset_index()
            .rename(columns={"subject": "n_subjects"})
            .sort_values("n_subjects", ascending=False)
        )
        sex_table = make_datatable(ss, page_size=10)
    else:
        missing_bits = []
        if "subject" not in base.columns:
            missing_bits.append("subject")
        if "sex" not in base.columns:
            missing_bits.append("sex")
        sex_table = html.Div("Missing metadata column: " + " or ".join(missing_bits))

    return overview_table, fig, stats_table, project_table, response_table, sex_table


if __name__ == "__main__":
    app.run(debug=True)

