"""Deterministic cleaning, sampling-event consolidation, and feature selection."""

import numpy as np
import pandas as pd

from .data import (
    DATE_COLUMNS, ENVIRONMENTAL_COLUMNS,
    GROUP_LABELS, IDENTIFIERS, NUMERIC_COLUMNS, SUBSTRATE_COLUMNS, TARGET,
)


def normalise_types(raw):
    """Strip text, replace case-insensitive nd tokens, and strictly parse known numbers."""
    frame = raw.copy(deep=True)
    rows = []
    for name in frame:
        series = frame[name]
        count = 0
        if not pd.api.types.is_numeric_dtype(series):
            series = series.astype("string").str.strip()
            marker = series.str.casefold().eq("nd").fillna(False)
            count = int(marker.sum())
            series = series.mask(marker)
        converted = name in NUMERIC_COLUMNS and not pd.api.types.is_numeric_dtype(raw[name])
        if name in NUMERIC_COLUMNS:
            # errors='raise' avoids silently turning unexpected numeric tokens into missing data.
            series = pd.to_numeric(series, errors="raise").astype(float)
            if np.isinf(series.to_numpy()).any():
                raise ValueError(f"Infinite values in {name}.")
        frame[name] = series
        rows.append({"column": name, "nd_replaced": count, "converted_to_numeric": converted,
                     "missing_after": int(series.isna().sum())})
    for name in IDENTIFIERS + DATE_COLUMNS[1:]:
        if frame[name].isna().any() or not frame[name].eq(np.floor(frame[name])).all():
            raise ValueError(f"{name} must contain complete integer values.")
        frame[name] = frame[name].astype("int64")
    return frame, pd.DataFrame(rows)


def validate_measurements(frame):
    """Fail on impossible recorded values without deleting plausible extremes."""
    checks = {
        "bleaching between 0 and 100": frame[TARGET].dropna().between(0, 100).all(),
        "cover between 0 and 100": frame["Percent_Cover"].dropna().between(0, 100).all(),
        "latitude in range": frame["Latitude_Degrees"].notna().all()
        and frame["Latitude_Degrees"].between(-90, 90).all(),
        "longitude in range": frame["Longitude_Degrees"].notna().all()
        and frame["Longitude_Degrees"].between(-180, 180).all(),
    }
    for name in ENVIRONMENTAL_COLUMNS:
        # Anomalies can be negative; distances, depths, variability, and counts cannot.
        nonnegative = name in ["Distance_to_Shore", "Turbidity", "Depth_m", "Cyclone_Frequency", "Windspeed"]
        nonnegative |= any(token in name for token in ["Standard_Deviation", "Frequency", "DHW"])
        if nonnegative:
            checks[f"{name} nonnegative"] = frame[name].dropna().ge(0).all()
        if name in ["ClimSST", "Temperature_Kelvin", "Temperature_Mean", "Temperature_Minimum", "Temperature_Maximum"]:
            checks[f"{name} positive Kelvin"] = frame[name].dropna().gt(0).all()
    result = pd.DataFrame({"check": list(checks), "passed": [bool(v) for v in checks.values()]})
    if not result.passed.all():
        raise ValueError(f"Invalid recorded measurements: {result.loc[~result.passed, 'check'].tolist()}")
    return result


def collapse_samples(frame):
    """Keep one event per Sample_ID only after checking all retained fields agree.

    Substrate-specific rows are not independent samples. Repeated outcome values
    are averaged, as in the original workflow; distinct event IDs are never merged.
    """
    exact_duplicates = int(frame.duplicated().sum())
    unique = frame.drop_duplicates()
    grouped = unique.groupby("Sample_ID", sort=True, dropna=False)
    variation = grouped.nunique(dropna=False).gt(1).sum()
    allowed = SUBSTRATE_COLUMNS + [TARGET]
    conflicts = variation.drop(labels=allowed, errors="ignore")
    if conflicts.gt(0).any():
        raise ValueError(f"Conflicting fields within Sample_ID: {conflicts[conflicts.gt(0)].to_dict()}")
    target_counts = grouped[TARGET].nunique(dropna=True)
    target_conflicts = grouped[TARGET].agg(["count", "min", "max", "mean"])
    target_conflicts = target_conflicts.loc[target_counts.gt(1)].reset_index()
    aggregate = {name: "first" for name in unique if name not in ["Sample_ID", *allowed]}
    aggregate[TARGET] = "mean"
    samples = grouped.agg(aggregate).reset_index()
    non_ids = samples.drop(columns=IDENTIFIERS)
    audit = pd.DataFrame([
        {"measure": "Raw rows", "value": len(frame)},
        {"measure": "Exact duplicate rows removed", "value": exact_duplicates},
        {"measure": "Distinct sampling events", "value": len(samples)},
        {"measure": "Events with differing recorded outcomes (averaged)", "value": len(target_conflicts)},
        {"measure": "Events identical on non-ID fields (retained)", "value": int(non_ids.duplicated().sum())},
    ])
    return samples, audit, target_conflicts


def retain_labelled_samples(samples):
    """Preserve the original labelled cohort; never impute bleaching outcomes."""
    frame = samples.loc[samples[TARGET].notna()].copy()
    if frame.empty or not frame[TARGET].between(0, 100).all():
        raise ValueError("No valid labelled sampling events, or bleaching outside 0–100%.")
    frame["Bleaching_Group"] = pd.cut(
        frame[TARGET], bins=[0, 10, 50, 100], labels=GROUP_LABELS,
        include_lowest=True, right=True, ordered=True,
    )
    return frame


def validate_dates(samples):
    """Parse ISO dates and require exact agreement with the supplied date parts."""
    frame = samples.copy()
    dates = pd.to_datetime(frame["Date"], format="%Y-%m-%d", errors="raise")
    if dates.isna().any():
        raise ValueError("Missing sampling dates.")
    for field, values in [("Date_Year", dates.dt.year), ("Date_Month", dates.dt.month),
                          ("Date_Day", dates.dt.day)]:
        if not frame[field].eq(values).all():
            raise ValueError(f"Date disagrees with {field}.")
    frame["Date"] = dates
    return frame


def filter_samples_from_year(samples, min_year=1998):
    """Retain events from min_year onward, before fitting any preprocessing.

    Validate agreement between the recorded date and year, and return an audit
    of excluded events. The original input frame is not modified.
    """
    if isinstance(min_year, bool) or not isinstance(min_year, (int, np.integer)):
        raise ValueError("The minimum sampling year must be an integer.")
    dated = validate_dates(samples)
    eligible = dated["Date_Year"].ge(min_year)
    excluded = dated.loc[~eligible, IDENTIFIERS + ["Date", "Date_Year"]].copy()
    retained = dated.loc[eligible].copy().reset_index(drop=True)
    if retained.empty:
        raise ValueError(f"No sampling events remain from {min_year} onward.")
    return retained, excluded


def select_features(samples, correlation_threshold=0.999):
    """Remove constants and later near-duplicate inputs using pairwise Spearman rho.

    Original source column order is retained. Other correlated environmental
    variables remain, preserving the original intent and weighting.
    """
    if not 0 < correlation_threshold <= 1:
        raise ValueError("Correlation threshold must be in (0, 1].")
    removed = {}
    for name in ENVIRONMENTAL_COLUMNS:
        if samples[name].nunique(dropna=True) <= 1:
            removed[name] = "Constant or entirely missing in the labelled cohort."
    candidates = [name for name in ENVIRONMENTAL_COLUMNS if name not in removed]
    correlation = samples[candidates].corr(method="spearman")
    pairs = []
    retained = []
    for name in candidates:
        duplicates = [previous for previous in retained
                      if abs(correlation.loc[previous, name]) >= correlation_threshold]
        if duplicates:
            previous = duplicates[0]
            removed[name] = f"Near-duplicate of {previous}: absolute Spearman rho >= {correlation_threshold}."
        else:
            retained.append(name)
    for i, left in enumerate(candidates):
        for right in candidates[i + 1:]:
            rho = correlation.loc[left, right]
            if abs(rho) >= 0.90:
                pairs.append({"feature_a": left, "feature_b": right, "spearman_rho": rho,
                              "near_duplicate": bool(abs(rho) >= correlation_threshold)})
    return retained, removed, pd.DataFrame(pairs, columns=[
        "feature_a", "feature_b", "spearman_rho", "near_duplicate"
    ]).sort_values("spearman_rho", key=abs, ascending=False)


def filter_sparse_samples(samples, features, max_missing_fraction=0.5):
    """Drop events missing more than half the selected inputs and flag missing depth."""
    if not 0 <= max_missing_fraction < 1:
        raise ValueError("Missingness threshold must be in [0, 1).")
    sparse = samples[features].isna().mean(axis=1).gt(max_missing_fraction)
    excluded = samples.loc[sparse, IDENTIFIERS].copy()
    excluded["missing_feature_count"] = samples.loc[sparse, features].isna().sum(axis=1)
    frame = samples.loc[~sparse].copy().reset_index(drop=True)
    frame["Depth_m_was_missing"] = frame["Depth_m"].isna().astype("int8")
    if frame.empty:
        raise ValueError("All sampling events failed the missingness rule.")
    return frame, excluded


def sample_metadata(samples, features):
    """Copy every non-PCA field, including comments and excluded environmental fields."""
    return samples[[name for name in samples if name not in features]].copy()


def substrate_records(typed_records, retained_samples):
    """Preserve original substrate rows for the final cohort without arbitrary averaging."""
    columns = IDENTIFIERS + SUBSTRATE_COLUMNS + [TARGET]
    return typed_records.loc[
        typed_records.Sample_ID.isin(retained_samples.Sample_ID), columns
    ].copy().reset_index(drop=True)
