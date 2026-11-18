"""Dataset schema, read-only input loading, and reproducible output writing."""

from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import re

import joblib
import pandas as pd


IDENTIFIERS = ["Sample_ID", "Site_ID"]
DATE_COLUMNS = ["Date", "Date_Day", "Date_Month", "Date_Year"]
GEO_COLUMNS = ["Latitude_Degrees", "Longitude_Degrees"]
TARGET = "Percent_Bleaching"
GROUP_LABELS = ["Low", "Moderate", "High"]
CATEGORICAL_COLUMNS = [
    "Data_Source", "Ocean_Name", "Reef_ID", "Realm_Name", "Ecoregion_Name",
    "Country_Name", "State_Island_Province_Name", "City_Town_Name", "Site_Name",
    "Exposure", "Substrate_Name", "Bleaching_Level",
]
TEXT_COLUMNS = ["Site_Comments", "Sample_Comments", "Bleaching_Comments"]
SUBSTRATE_COLUMNS = ["Substrate_Name", "Percent_Cover"]
ENVIRONMENTAL_COLUMNS = [
    "Distance_to_Shore", "Turbidity", "Cyclone_Frequency", "Depth_m", "ClimSST",
    "Temperature_Kelvin", "Temperature_Mean", "Temperature_Minimum",
    "Temperature_Maximum", "Temperature_Kelvin_Standard_Deviation", "Windspeed",
    "SSTA", "SSTA_Standard_Deviation", "SSTA_Mean", "SSTA_Minimum", "SSTA_Maximum",
    "SSTA_Frequency", "SSTA_Frequency_Standard_Deviation", "SSTA_FrequencyMax",
    "SSTA_FrequencyMean", "SSTA_DHW", "SSTA_DHW_Standard_Deviation", "SSTA_DHWMax",
    "SSTA_DHWMean", "TSA", "TSA_Standard_Deviation", "TSA_Minimum", "TSA_Maximum",
    "TSA_Mean", "TSA_Frequency", "TSA_Frequency_Standard_Deviation",
    "TSA_FrequencyMax", "TSA_FrequencyMean", "TSA_DHW", "TSA_DHW_Standard_Deviation",
    "TSA_DHWMax", "TSA_DHWMean",
]
NUMERIC_COLUMNS = (
    IDENTIFIERS + GEO_COLUMNS + DATE_COLUMNS[1:] + [TARGET, "Percent_Cover"]
    + ENVIRONMENTAL_COLUMNS
)
RAW_COLUMNS = set(NUMERIC_COLUMNS + CATEGORICAL_COLUMNS + TEXT_COLUMNS + ["Date"])
INTERPRETATION_COLUMNS = (
    IDENTIFIERS + [TARGET, "Bleaching_Group", "Depth_m_was_missing"]
    + DATE_COLUMNS + GEO_COLUMNS
)


def file_sha256(path):
    """Return a streaming SHA-256 fingerprint without modifying the file."""
    digest = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_project_root(start=None):
    """Find this project from its root or a child working directory."""
    start = Path.cwd() if start is None else Path(start).resolve()
    for candidate in (start, *start.parents):
        if (candidate / "coral.csv").is_file() and (candidate / "coral_pipeline").is_dir():
            return candidate
    raise FileNotFoundError("Run from the project folder or a child folder containing coral.csv.")


def load_raw(path):
    """Read coral.csv and reject schema changes instead of guessing column roles."""
    frame = pd.read_csv(path, low_memory=False)
    missing, unexpected = RAW_COLUMNS - set(frame), set(frame) - RAW_COLUMNS
    if missing or unexpected or not frame.columns.is_unique:
        raise ValueError(f"Unexpected raw schema: missing={sorted(missing)}, extra={sorted(unexpected)}")
    if frame.empty:
        raise ValueError("The raw dataset is empty.")
    return frame


def load_pca_inputs(root):
    """Load the combined export and its exact ordered 35-feature metadata list.

    No numeric inference, imputation, or automatic rescaling is performed here.
    Categories and comments remain in the returned combined table.
    """
    import numpy as np

    output = Path(root) / "processed"
    metadata = json.loads((output / "coral_preprocessing_metadata.json").read_text(encoding="utf-8"))
    features = metadata.get("pca_input_columns", [])
    if len(features) != 35 or len(set(features)) != 35:
        raise ValueError("Preprocessing metadata must specify 35 unique ordered PCA inputs.")
    frame = pd.read_csv(output / "coral_pca_ready.csv", low_memory=False)
    required = features + IDENTIFIERS + DATE_COLUMNS + GEO_COLUMNS + [
        TARGET, "Bleaching_Group", "Depth_m_was_missing", "Data_Source", "Exposure",
        "Ocean_Name", "Realm_Name", "Country_Name", "Ecoregion_Name",
    ]
    missing = sorted(set(required) - set(frame))
    if missing:
        raise ValueError(f"PCA input table is missing required columns: {missing}")
    excluded = set(metadata.get("interpretation_columns_not_for_pca", []))
    if excluded.intersection(features) or not set(features).issubset(ENVIRONMENTAL_COLUMNS):
        raise ValueError("Invalid feature roles in preprocessing metadata.")
    values = frame[features]
    if not all(pd.api.types.is_numeric_dtype(values[name]) for name in features):
        raise ValueError("All PCA inputs must be numeric.")
    if not np.isfinite(values.to_numpy()).all() or values.nunique().le(1).any():
        raise ValueError("PCA inputs contain missing, infinite, or constant values.")
    if not np.allclose(values.mean(), 0, atol=1e-6) or not np.allclose(values.std(ddof=0), 1, atol=1e-6):
        raise ValueError("Inputs are not standardised; rerun preprocessing instead of scaling silently.")
    from .cleaning import validate_dates
    frame = validate_dates(frame)
    if (not frame.Sample_ID.is_unique or frame.Sample_ID.isna().any()
            or frame.Site_ID.isna().any() or not frame.Date_Year.ge(1998).all()
            or not frame.Bleaching_Group.isin(GROUP_LABELS).all()):
        raise ValueError("Invalid sampling events, groups, or 1998-onward cohort.")
    expected = pd.cut(frame[TARGET], [0, 10, 50, 100], labels=GROUP_LABELS, include_lowest=True)
    if not expected.astype(str).eq(frame.Bleaching_Group).all():
        raise ValueError("Bleaching outcomes and groups disagree.")
    return frame, values.copy(), metadata


def load_cv_measurements(root, samples, features):
    """Recover missing original-unit inputs for fold-specific preprocessing in CV.

    Cohort eligibility and the study's 35 named variables remain fixed. Medians,
    skew decisions, power parameters, scaling and PCA must be fitted inside folds.
    """
    from .cleaning import normalise_types, collapse_samples
    typed, _ = normalise_types(load_raw(Path(root) / "coral.csv"))
    raw_samples, _, _ = collapse_samples(typed)
    aligned = raw_samples.set_index("Sample_ID").loc[samples.Sample_ID, features].copy()
    aligned.index = samples.index
    return aligned


def load_dictionary(path):
    """Read the supplied fixed-width Markdown dictionary into descriptions and units."""
    entries, spans, current = {}, None, None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if re.fullmatch(r"-{3,}( +-{3,}){2,}", text):
            spans = [(match.start(), match.end()) for match in re.finditer(r"-+", line)]
            current = None
            continue
        if not text or set(text) == {"-"}:
            current = None
            continue
        if spans is None or not line.startswith("  "):
            continue
        cells = [line[a:b].strip() for a, b in spans[:-1]] + [line[spans[-1][0]:].strip()]
        if cells[0].startswith("`"):
            current = cells[0].strip("`")
            entries[current] = {"description": cells[1], "units": cells[2]}
        elif current and not cells[0]:
            entries[current]["description"] += " " + cells[1]
            entries[current]["units"] += " " + cells[2]
    result = pd.DataFrame.from_dict(entries, orient="index")
    return result.apply(lambda col: col.str.replace(r"\s+", " ", regex=True).str.strip())


def column_roles(columns, dictionary, removed=None):
    """Distinguish column preservation from eligibility for the numeric PCA matrix."""
    removed = {} if removed is None else removed
    rows = []
    for name in columns:
        if name in IDENTIFIERS:
            role, reason = "Identifier", "Trace records; exclude from environmental inputs."
        elif name in DATE_COLUMNS + GEO_COLUMNS:
            role, reason = "Interpretation", "Retain date or location; exclude from environmental inputs."
        elif name == TARGET:
            role, reason = "Outcome", "Retain for later interpretation; never impute or use as an input."
        elif name in SUBSTRATE_COLUMNS:
            role, reason = "Substrate detail (separate table)", "Multiple substrate records per event; preserve in coral_substrate_records.csv."
        elif name in CATEGORICAL_COLUMNS:
            role, reason = "Category (retained)", "Keep as metadata; no numeric encoding is assumed."
            if name == "Exposure":
                reason = "Sometimes is a distinct exposure mechanism, not a midpoint on a numeric scale."
        elif name in TEXT_COLUMNS:
            role, reason = "Text (retained)", "Keep as metadata only; bleaching comments may reveal outcome derivation."
        else:
            role, reason = "Environmental input", "Numeric environmental measurement."
        if name in removed:
            role, reason = "Environmental metadata (retained)", removed[name] + " Retain in original units, outside PCA inputs."
        rows.append({"column": name, "role": role, "reason": reason,
                     "in_sample_table": name not in SUBSTRATE_COLUMNS,
                     "in_pca_input": name in ENVIRONMENTAL_COLUMNS and name not in removed})
    return pd.DataFrame(rows).set_index("column").join(dictionary)


def save_outputs(root, raw_hash, cleaned, ready, sample_metadata, bundle, metadata, tables):
    """Write derived files only, with feature order, fitted state, audits, and fingerprints."""
    root = Path(root)
    if file_sha256(root / "coral.csv") != raw_hash:
        raise ValueError("Raw input changed during processing; outputs have not been written.")
    output = root / "processed"
    output.mkdir(exist_ok=True)
    artifacts = {
        "coral_cleaned.csv": cleaned,
        "coral_pca_ready.csv": ready,
        "coral_sample_metadata.csv": sample_metadata,
        "coral_pca_features.csv": ready[bundle["pca_input_columns"]],
        **{f"{name}.csv": table for name, table in tables.items()},
    }
    # All audit tables have ordinary columns, including their row identifiers.
    for name, frame in artifacts.items():
        frame.to_csv(output / name, index=False, date_format="%Y-%m-%d")
    transformer_path = output / "coral_transformers.joblib"
    joblib.dump(bundle, transformer_path)
    metadata = {
        **metadata,
        "source_file": "coral.csv", "source_sha256": raw_hash,
        "rows": len(ready), "unit_of_analysis": "one row per Sample_ID",
        "pca_input_columns": bundle["pca_input_columns"],
        "interpretation_columns_not_for_pca": [name for name in ready if name not in bundle["pca_input_columns"]],
        "combined_table_columns": list(ready.columns),
        "metadata_missing_values": "Preserved; only selected PCA input gaps are imputed.",
        "package_versions": {name: version(name) for name in
                             ["numpy", "pandas", "scipy", "scikit-learn", "matplotlib", "joblib"]},
        "output_sha256": {name: file_sha256(output / name) for name in artifacts},
        "transformers_sha256": file_sha256(transformer_path),
    }
    metadata_path = output / "coral_preprocessing_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return pd.DataFrame([
        {"file": name, "rows": len(frame), "columns": len(frame.columns)}
        for name, frame in artifacts.items()
    ])


def preprocessing_metadata(roles, transformations, retention, max_missing_fraction=0.5,
                           correlation_threshold=0.999, min_sample_year=1998):
    """Collect methodological decisions and executed audit values for downstream use."""
    excluded = roles.loc[~roles.in_pca_input, "reason"].to_dict()
    removed = roles.loc[~roles.in_sample_table, "reason"].to_dict()
    return {
        "scope": "Cleaning and preprocessing only; no PCA or subsequent analysis fitted.",
        "reference_population": f"Retained sampling events from {min_sample_year} onward with an observed bleaching outcome.",
        "minimum_sampling_year_inclusive": min_sample_year,
        "sampling_year_filter": "Applied after outcome availability and date validation, before feature selection or fitted preprocessing.",
        "fit_scope": "Full retained cohort for subsequent descriptive PCA.",
        "predictive_validation_requirement": (
            "Split raw sample-level data by site, and by time if needed, before fitting feature selection, "
            "medians, power transforms, scaling, and subsequent dimensionality reduction. "
            "Do not cross-validate predictions using this globally fitted table."
        ),
        "source_sha256_in_original_notebook": "78cc014b1e7887f26f4c44a40d4d9213212f3c6217541c6a1b69e723b61bfacb",
        "missing_tokens": ["nd"], "imputation": "Column medians fitted after sparse-event exclusion.",
        "max_missing_fraction_per_event": max_missing_fraction,
        "near_duplicate_absolute_spearman_threshold": correlation_threshold,
        "other_high_correlations_pruned": False,
        "outliers_removed_or_clipped": False,
        "outcome_aggregation": "Arithmetic mean across distinct substrate records within Sample_ID.",
        "exact_duplicate_policy": "Remove exact duplicate records before sample-level aggregation.",
        "bleaching_groups": {"Low": "0 <= value <= 10", "Moderate": "10 < value <= 50", "High": "50 < value <= 100"},
        "removed_columns": removed,
        "excluded_from_pca_columns": excluded,
        "column_preservation_policy": "Retain all sample-level source columns. Preserve substrate-specific fields in a separate detail table. Transform and impute only the selected PCA inputs.",
        "power_transformed": transformations.loc[transformations.method.eq("yeo-johnson"), "feature"].tolist(),
        "standardised_only": transformations.loc[~transformations.method.eq("yeo-johnson"), "feature"].tolist(),
        "transformations": json.loads(transformations.to_json(orient="records", double_precision=15)),
        "cohort_retention": json.loads(retention.to_json(orient="records", double_precision=15)),
        "limitations": [
            "Outcome-based exclusion preserves the original research cohort but can introduce selection bias.",
            "Several events may share a site; one row per event does not establish statistical independence.",
            "Median imputation can alter covariance; source-related depth missingness is retained in metadata.",
            "Source documentation substitutes day 15 when sampling day is unknown; prefer year and month for interpretation.",
            "Repeated thermal-summary families retain unequal weighting; zero-heavy inputs need not become normal.",
            "Historical summary windows may include later observations; assess temporal availability before prediction.",
            "Windspeed units in the supplied dictionary require verification before physical interpretation.",
        ],
    }
