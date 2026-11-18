"""Descriptive tables, captions, and computed prose for the final study report.

These helpers present existing data and PCA results. They do not perform new
cleaning, feature selection, statistical testing, or PCA fitting.
"""

from collections.abc import Mapping
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from IPython.display import Markdown, display

from .data import IDENTIFIERS


def raw_data_profile(
    raw: pd.DataFrame, typed: pd.DataFrame, dictionary: pd.DataFrame
) -> pd.DataFrame:
    """Profile every raw field after recognising the existing missing markers.

    ``raw`` and ``typed`` must have the same rows and columns. Missingness is
    measured on ``typed`` before aggregation, filtering, or imputation. Numeric
    summaries are omitted for identifiers, whose useful summary is distinctness.
    The original dtype and the number of literal ``nd`` markers remain visible.
    """
    if list(raw.columns) != list(typed.columns) or not raw.index.equals(typed.index):
        raise ValueError("Raw and typed frames must have matching rows and columns.")
    if raw.empty:
        raise ValueError("The raw dataset must contain observations.")
    rows = []
    for name in raw:
        values = typed[name]
        numeric = pd.api.types.is_numeric_dtype(values)
        valid = values.dropna()
        frequencies = valid.value_counts()
        numeric_measurement = numeric and name not in IDENTIFIERS
        rows.append(
            {
                "variable": name,
                "raw_dtype": str(raw[name].dtype),
                "typed_dtype": str(values.dtype),
                "nonmissing_records": int(values.notna().sum()),
                "missing_records": int(values.isna().sum()),
                "missing_percent": float(100 * values.isna().mean()),
                "nd_markers": int(
                    raw[name].astype("string").str.strip().str.casefold().eq("nd").sum()
                ),
                "distinct_nonmissing": int(values.nunique(dropna=True)),
                "minimum": valid.min() if numeric_measurement and len(valid) else np.nan,
                "median": valid.median() if numeric_measurement and len(valid) else np.nan,
                "maximum": valid.max() if numeric_measurement and len(valid) else np.nan,
                "most_frequent": frequencies.index[0] if len(frequencies) else pd.NA,
                "most_frequent_count": int(frequencies.iloc[0]) if len(frequencies) else 0,
            }
        )
    result = pd.DataFrame(rows).set_index("variable")
    available = [name for name in ["description", "units"] if name in dictionary]
    return result.join(dictionary[available], how="left").reset_index()


def source_collection_summary(typed: pd.DataFrame) -> pd.DataFrame:
    """Describe recorded source coverage without inferring collection procedures.

    Counts and missingness percentages refer to raw substrate records. Distinct
    event and site counts are reported separately because record multiplicity
    varies between sources. Collection methods require external documentation.
    """
    required = {"Data_Source", "Sample_ID", "Site_ID", "Date_Year", "Percent_Bleaching", "Depth_m"}
    missing = sorted(required - set(typed))
    if missing:
        raise ValueError(f"Source summary requires columns: {missing}")
    return (
        typed.groupby("Data_Source", dropna=False)
        .agg(
            substrate_records=("Sample_ID", "size"),
            distinct_events=("Sample_ID", "nunique"),
            distinct_sites=("Site_ID", "nunique"),
            first_recorded_year=("Date_Year", "min"),
            last_recorded_year=("Date_Year", "max"),
            missing_bleaching_record_percent=("Percent_Bleaching", lambda s: 100 * s.isna().mean()),
            missing_depth_record_percent=("Depth_m", lambda s: 100 * s.isna().mean()),
        )
        .sort_values("substrate_records", ascending=False)
        .reset_index()
    )


def cohort_audit(
    raw: pd.DataFrame,
    samples: pd.DataFrame,
    labelled: pd.DataFrame,
    period_samples: pd.DataFrame,
    retained: pd.DataFrame,
) -> pd.DataFrame:
    """Reconcile counts and stage percentages while distinguishing the unit change.

    The raw-to-event transition reduces substrate records through consolidation;
    its percentage is a record reduction, not sampling-event exclusion. Later
    removal/retention percentages use events entering each filtering stage.
    """
    counts = [len(raw), len(samples), len(labelled), len(period_samples), len(retained)]
    if any(count == 0 for count in counts) or any(
        later > earlier for earlier, later in zip(counts, counts[1:])
    ):
        raise ValueError("Cohort stages must be nonempty and cannot increase in size.")
    descriptions = [
        ("Raw observations", "substrate records", "Source dataset"),
        ("Consolidated sampling events", "sampling events", "One event per Sample_ID"),
        (
            "Observed bleaching outcome",
            "sampling events",
            "Exclude missing outcomes; do not impute",
        ),
        ("Sampling year 1998 onward", "sampling events", "Exclude 1997 and earlier"),
        ("Complete enough environmental profile", "sampling events", "Exclude >50% missing inputs"),
    ]
    rows = []
    for i, ((stage, unit, rule), count) in enumerate(zip(descriptions, counts)):
        change = counts[i - 1] - count if i else 0
        event_filter = i >= 2
        rows.append(
            {
                "stage": stage,
                "unit": unit,
                "retained_count": count,
                "rule": rule,
                "removed_or_consolidated_count": change,
                "removed_percent_of_entering_events": (
                    100 * change / counts[i - 1] if event_filter else np.nan
                ),
                "retained_percent_of_entering_events": (
                    100 * count / counts[i - 1] if event_filter else np.nan
                ),
                "record_reduction_percent": 100 * change / counts[0] if i == 1 else np.nan,
                "percent_of_distinct_events": (100 * count / counts[1] if i >= 1 else np.nan),
            }
        )
    return pd.DataFrame(rows)


def workflow_figure() -> Figure:
    """Return the study-specific workflow diagram without hard-coded PCA findings."""
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    width, height = 0.28, 0.29
    positions = [(0.03, 0.60), (0.36, 0.60), (0.69, 0.60), (0.69, 0.12), (0.36, 0.12), (0.03, 0.12)]
    labels = [
        "1. Raw coral observations\nSubstrate records, survey sources\nand environmental measurements",
        "2. Sampling-event consolidation\nRecognise missing markers;\nverify agreement within Sample_ID",
        "3. Study cohort\nObserved bleaching outcome;\nsampling year 1998 onward",
        "4. Environmental preparation\nSparse-profile exclusion; medians;\nselective power transform + scaling",
        "5. Environmental PCA\n35 explicit inputs; parallel analysis;\nscores, loadings and contributions",
        "6. Interpretation and checks\nBleaching/source comparisons;\nsite, robustness and time audits",
    ]
    for (x, y), label in zip(positions, labels):
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.01,rounding_size=0.015",
            facecolor="#e7f0f8",
            edgecolor="#2a78d6",
            linewidth=1.4,
        )
        ax.add_patch(patch)
        ax.text(x + width / 2, y + height / 2, label, ha="center", va="center", fontsize=10)
    connectors = [
        ((0.32, 0.745), (0.35, 0.745)),
        ((0.65, 0.745), (0.68, 0.745)),
        ((0.83, 0.59), (0.83, 0.42)),
        ((0.68, 0.265), (0.65, 0.265)),
        ((0.35, 0.265), (0.32, 0.265)),
    ]
    for start, end in connectors:
        ax.add_patch(
            FancyArrowPatch(
                start, end, arrowstyle="-|>", mutation_scale=15, linewidth=1.5, color="#40586e"
            )
        )
    fig.tight_layout()
    return fig


def display_table(frame: pd.DataFrame, number: str | int, title: str, note: str = "") -> None:
    """Display a numbered table and optional explanatory note without changing data."""
    display(Markdown(f"**Table {number}. {title}**"))
    # Wrap complete descriptions instead of truncating important report cells.
    styled = frame.style.format(
        formatter=lambda value: f"{value:,.4g}" if isinstance(value, (float, np.floating)) else str(value),
        na_rep="—", escape="html",
    ).set_properties(**{"white-space": "normal", "max-width": "360px", "text-align": "left"})
    styled = styled.set_table_styles([
        {"selector": "th", "props": [("white-space", "normal"), ("text-align", "left"),
                                        ("background-color", "#edf3f8"), ("padding", "6px")]},
        {"selector": "td", "props": [("padding", "6px"), ("vertical-align", "top")]},
    ])
    if isinstance(frame.index, pd.RangeIndex) and frame.index.name is None:
        styled = styled.hide(axis="index")
    display(styled)
    if note:
        display(Markdown(f"*Note.* {note}"))


def display_figure(fig: Figure, number: str | int, title: str, note: str = "") -> None:
    """Display a numbered figure with a caption, then close its interactive handle."""
    try:
        display(fig)
        display(Markdown(f"**Figure {number}. {title}**" + (f"\n\n{note}" if note else "")))
    finally:
        plt.close(fig)


def _report_values(samples: pd.DataFrame, result: Mapping[str, Any], tables: Mapping[str, Any]):
    """Read the existing result objects used by the abstract and conclusion."""
    count = int(result["n_retained"])
    variance = float(result["variance"].iloc[count - 1]["cumulative_%"])
    group = tables["pca_group_separation"]
    names = group.cliffs_delta.abs().nlargest(2).index
    effects = "; ".join(f"{name}: δ = {group.loc[name, 'cliffs_delta']:+.2f}" for name in names)
    return count, variance, names, effects


def abstract_markdown(
    samples: pd.DataFrame, result: Mapping[str, Any], tables: Mapping[str, Any]
) -> str:
    """Return a short study abstract with current cohort, methods, findings and recommendation."""
    count, variance, names, effects = _report_values(samples, result, tables)
    features = len(result["loadings"])
    return (
        "Coral survey datasets contain correlated environmental measurements that can obscure "
        "the patterns associated with bleaching severity. This study asks whether principal "
        "component analysis (PCA) can provide an interpretable summary of those patterns. "
        f"The analysis used {len(samples):,} sampling events from "
        f"{int(samples.Date_Year.min())}–{int(samples.Date_Year.max())} and {features} "
        "environmental inputs. Sampling events were consolidated, remaining input gaps were "
        "median-imputed, and selected power transformations and standardisation preceded PCA. "
        "Parallel analysis guided component retention; bleaching-group contrasts, source "
        "stratification and sensitivity checks supported interpretation. "
        f"The {count} retained components explained {variance:.1f}% of environmental variance. "
        f"The strongest pooled High-versus-Low score differences were {effects}. "
        "The components summarise environmental variation but do not fully preserve observed "
        "bleaching-related separation. Original indicators should therefore be retained alongside "
        "PCA summaries, and associations should be validated across independent sources and "
        "periods before operational use."
    )


def conclusion_markdown(
    samples: pd.DataFrame, result: Mapping[str, Any], tables: Mapping[str, Any]
) -> str:
    """Answer the study question using current numerical results and stated inferential limits."""
    count, variance, names, effects = _report_values(samples, result, tables)
    label = " and ".join(names)
    text = (
        "**Answer to the research question: partly.** For this cohort, PCA reveals an "
        "interpretable environmental structure associated with bleaching severity, but it "
        "provides an incomplete summary of the observed bleaching-group differences. "
        f"Across {len(samples):,} events, {count} retained components explain {variance:.1f}% "
        f"of environmental variance. {label} show the strongest pooled High-versus-Low "
        f"contrasts ({effects}). These associations should be read with the source-stratified "
        "comparisons and uncertainty audits, because component scores also vary with source "
        "and geography."
    )
    if "hotelling_high_vs_low" in tables:
        ratio = float(tables["hotelling_high_vs_low"].share_of_all_variables_D2.iloc[1])
        text += (
            f" The retained components preserve {100 * ratio:.1f}% of the sample-level squared "
            "Mahalanobis separation between High and Low relative to all original inputs. "
            "This is a descriptive separation ratio, not a measure of predictive accuracy."
        )
    return text + (
        " The findings support retaining original environmental indicators alongside the "
        "component summaries. They do not establish causal effects, operational warning "
        "thresholds, or independent future-year performance."
    )
