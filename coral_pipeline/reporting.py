"""Concise quality tables, preprocessing figures, and PCA-readiness validation."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .data import GROUP_LABELS, INTERPRETATION_COLUMNS, RAW_COLUMNS, SUBSTRATE_COLUMNS, TARGET


def configure_display():
    """Use restrained report styling without hiding warnings or changing data."""
    pd.set_option("display.max_rows", 70, "display.max_columns", 12,
                  "display.max_colwidth", 75, "display.float_format", "{:,.4f}".format)
    plt.rcParams.update({
        "figure.dpi": 110, "axes.spines.top": False, "axes.spines.right": False,
        "axes.axisbelow": True, "axes.grid": True, "grid.color": "#e8e8e8",
        "font.size": 10, "axes.titlesize": 11, "figure.constrained_layout.use": True,
    })


def missingness_table(frame, columns):
    """Count missing inputs before imputation, in descending order."""
    counts = frame[columns].isna().sum()
    result = pd.DataFrame({"feature": columns, "missing": counts.to_numpy(),
                           "missing_percent": 100 * counts.to_numpy() / len(frame)})
    return result.sort_values("missing", ascending=False, kind="stable").reset_index(drop=True)


def plot_missingness(table):
    """Show the twelve most affected inputs; the full audit retains every feature."""
    values = table.loc[table.missing.gt(0)].head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, max(2.5, 0.32 * len(values) + 1)))
    ax.barh(values.feature, values.missing_percent, color="#287b8e")
    ax.set(xlabel="Sampling events with missing values (%)",
           title="Inputs with greatest missingness before sparse-event exclusion")
    ax.grid(axis="y", visible=False)
    return fig


def group_summary(frame):
    """Describe cohort sizes and imputed depth by the original bleaching cut-points."""
    grouped = frame.groupby("Bleaching_Group", observed=False)
    result = grouped.size().reindex(GROUP_LABELS).rename("samples").to_frame()
    result["percent"] = 100 * result.samples / len(frame)
    if "Depth_m_was_missing" in frame:
        result["depth_imputed"] = grouped["Depth_m_was_missing"].sum().reindex(GROUP_LABELS)
        result["depth_imputed_percent"] = 100 * result.depth_imputed / result.samples
    return result.rename_axis("bleaching_group").reset_index()


def depth_missing_by_source(frame):
    """Retain the source-related missingness caveat without comparing environmental effects."""
    return frame.groupby("Data_Source", dropna=False).agg(
        samples=("Sample_ID", "size"), depth_missing=("Depth_m_was_missing", "sum")
    ).assign(depth_missing_percent=lambda t: 100 * t.depth_missing / t.samples).reset_index()


def profile_features(frame, features):
    """Summarise original-unit ranges, skew, zeros, and IQR extremes; do not trim data."""
    rows = []
    for name in features:
        values = frame[name].dropna()
        q1, q3 = values.quantile([0.25, 0.75])
        iqr = q3 - q1
        extreme = (values.lt(q1 - 1.5 * iqr) | values.gt(q3 + 1.5 * iqr)).mean()
        rows.append({"feature": name, "min": values.min(), "median": values.median(),
                     "max": values.max(), "skew": values.skew(),
                     "zero_percent": 100 * values.eq(0).mean(),
                     "iqr_extreme_percent": 100 * extreme if iqr > 0 else np.nan})
    return pd.DataFrame(rows)


def plot_transform_diagnostics(audit, imputed, scaled):
    """Display skew reduction and two illustrative before/after distributions."""
    top = audit.sort_values("skew_before", key=abs, ascending=False).head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    positions = np.arange(len(top))
    ax.hlines(positions, top.skew_after.abs(), top.skew_before.abs(), color="#c9d0d2")
    ax.scatter(top.skew_before.abs(), positions, label="Before", color="#287b8e", zorder=3)
    ax.scatter(top.skew_after.abs(), positions, label="After", color="#c56b3f", zorder=3)
    ax.set_yticks(positions, top.feature)
    ax.set(xlabel="Absolute sample skewness", title="Twelve inputs with the greatest initial skew")
    ax.legend()
    ax.grid(axis="y", visible=False)

    examples = ["Distance_to_Shore", "TSA_DHW"]
    examples = [name for name in examples if name in imputed]
    distributions, axes = plt.subplots(2, len(examples), figsize=(10, 5.5), squeeze=False)
    for column, name in enumerate(examples):
        axes[0, column].hist(imputed[name], bins=40, color="#287b8e")
        axes[0, column].set(title=name, xlabel="Imputed measurement in original units", ylabel="Events")
        axes[1, column].hist(scaled[name], bins=40, color="#c56b3f")
        axes[1, column].set(xlabel="Transformed, standardised value", ylabel="Events")
    return fig, distributions


def retention_table(raw, samples, labelled, period_samples, retained, min_sample_year=1998):
    """Reconcile each cohort change without double-counting substrate rows as samples."""
    return pd.DataFrame([
        {"stage": "Raw substrate records", "rows": len(raw), "removed_from_previous": 0},
        {"stage": "One row per Sample_ID", "rows": len(samples),
         "removed_from_previous": len(raw) - len(samples)},
        {"stage": "Observed bleaching outcome", "rows": len(labelled),
         "removed_from_previous": len(samples) - len(labelled)},
        {"stage": f"Sampling year >= {min_sample_year}", "rows": len(period_samples),
         "removed_from_previous": len(labelled) - len(period_samples)},
        {"stage": "At most 50% of environmental inputs missing", "rows": len(retained),
         "removed_from_previous": len(period_samples) - len(retained)},
    ])


def validate_ready(cleaned, ready, features, min_sample_year=None):
    """Require a finite standardised feature matrix and consistent interpretation data.

    These checks establish numerical suitability; they do not fit PCA, guarantee
    normality, or establish independence among samples from the same site.
    """
    values = ready[features]
    array = values.to_numpy(dtype=float)
    target = ready[TARGET]
    expected_groups = pd.cut(target, [0, 10, 50, 100], labels=GROUP_LABELS, include_lowest=True)
    dates = pd.to_datetime(ready.Date, format="%Y-%m-%d", errors="raise")
    numeric = all(pd.api.types.is_numeric_dtype(values[name]) for name in features)
    metadata_columns = [name for name in ready if name not in features]
    checks = {
        "Feature list is nonempty and unique": bool(features) and len(set(features)) == len(features),
        "No metadata or outcome in the feature matrix": not set(features).intersection(INTERPRETATION_COLUMNS),
        "Every feature has a numeric dtype": numeric,
        "No missing or infinite feature values": bool(np.isfinite(array).all()),
        "No constant feature": bool(values.nunique().gt(1).all()),
        "All feature means within 1e-6 of zero": bool(np.allclose(array.mean(axis=0), 0, atol=1e-6)),
        "All population standard deviations within 1e-6 of one": bool(np.allclose(array.std(axis=0), 1, atol=1e-6, rtol=0)),
        "One complete Sample_ID per row": ready.Sample_ID.notna().all() and ready.Sample_ID.is_unique,
        "Complete Site_ID values": ready.Site_ID.notna().all(),
        "No missing or out-of-range bleaching outcomes": target.notna().all() and target.between(0, 100).all(),
        "Bleaching groups match the declared cut-points": ready.Bleaching_Group.notna().all()
        and ready.Bleaching_Group.astype(str).eq(expected_groups.astype(str)).all(),
        "Depth imputation flag is binary": ready.Depth_m_was_missing.isin([0, 1]).all(),
        "Cleaned and standardised rows align": cleaned.Sample_ID.equals(ready.Sample_ID),
        "All non-PCA fields are unchanged by transformation": cleaned[metadata_columns].equals(ready[metadata_columns]),
        "Dates agree with year, month, and day": dates.notna().all()
        and dates.dt.year.eq(ready.Date_Year).all() and dates.dt.month.eq(ready.Date_Month).all()
        and dates.dt.day.eq(ready.Date_Day).all(),
        "Latitude and longitude remain in valid ranges": ready.Latitude_Degrees.between(-90, 90).all()
        and ready.Longitude_Degrees.between(-180, 180).all(),
        "More observations than inputs": len(ready) > len(features),
        "Cleaned and standardised schemas match": list(cleaned) == list(ready),
        "All sample-level source columns are preserved": (RAW_COLUMNS - set(SUBSTRATE_COLUMNS)).issubset(ready.columns),
    }
    if min_sample_year is not None:
        checks[f"All events sampled from {min_sample_year} onward"] = (
            ready.Date_Year.ge(min_sample_year).all() and dates.dt.year.ge(min_sample_year).all()
        )
    result = pd.DataFrame({"check": list(checks), "passed": [bool(value) for value in checks.values()]})
    if not result.passed.all():
        raise ValueError(f"Readiness checks failed: {result.loc[~result.passed, 'check'].tolist()}")
    return result
