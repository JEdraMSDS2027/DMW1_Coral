"""Fit median imputation and per-feature transforms; reuse fitted state on new data."""

import warnings

import numpy as np
import pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PowerTransformer, RobustScaler, StandardScaler


def build_processed_tables(samples, imputed, scaled):
    """Preserve sample columns, replacing only the selected feature measurements.

    The cleaned copy uses imputed original-unit inputs; the ready copy uses the
    transformed inputs. All non-PCA values and missingness remain unchanged.
    """
    if (list(imputed.columns) != list(scaled.columns)
            or not samples.index.equals(imputed.index)
            or not samples.index.equals(scaled.index)):
        raise ValueError("Sample rows and feature tables must align before constructing outputs.")
    if not set(imputed.columns).issubset(samples.columns):
        raise ValueError("Feature tables contain columns absent from the sample table.")
    cleaned = samples.copy(deep=True)
    cleaned[imputed.columns] = imputed
    ready = cleaned.copy(deep=True)
    ready[scaled.columns] = scaled
    return cleaned, ready


def _fit_power(values):
    """Try the original Yeo–Johnson fit, then median/IQR rescaling for stability."""
    failures = []
    for rescale in (False, True):
        power = PowerTransformer(method="yeo-johnson", standardize=False)
        model = make_pipeline(RobustScaler(), power, StandardScaler()) if rescale else make_pipeline(power, StandardScaler())
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", RuntimeWarning)
                transformed = model.fit_transform(values).ravel()
            lam = float(power.lambdas_[0])
            acceptable = (np.isfinite(transformed).all() and abs(transformed.mean()) < 1e-6
                          and abs(transformed.std(ddof=0) - 1) < 1e-6 and abs(lam) <= 10)
            if acceptable:
                return model, transformed, lam, rescale, "; ".join(failures)
            failures.append(f"{'rescaled' if rescale else 'direct'} fit failed stability criteria")
            failures.extend(str(warning.message) for warning in caught)
        except (ValueError, FloatingPointError, OverflowError, RuntimeError) as exc:
            failures.append(f"{'rescaled' if rescale else 'direct'} fit: {exc}")
    return None, None, None, False, "; ".join(failures)


def fit_preprocessing(samples, features, skew_threshold=0.75):
    """Fit on the supplied reference cohort and return cleaned/scaled inputs and state.

    For descriptive work the full retained cohort is the reference. For predictive
    validation, call this only on training data after splitting by site or time;
    feature selection must also use that training data alone.
    """
    if skew_threshold < 0 or not features or len(set(features)) != len(features):
        raise ValueError("Use a nonnegative skew threshold and unique, nonempty features.")
    original = samples[features].astype(float)
    if np.isinf(original.to_numpy()).any():
        raise ValueError("Infinite environmental measurements.")
    medians = original.median()
    if medians.isna().any():
        raise ValueError(f"Cannot impute all-missing features: {medians[medians.isna()].index.tolist()}")
    imputed = original.fillna(medians)
    if imputed.nunique().le(1).any():
        raise ValueError("A retained feature is constant after cohort filtering; review feature selection.")
    scaled = pd.DataFrame(index=samples.index)
    models, rows = {}, []
    for name in features:
        values = imputed[[name]].to_numpy()
        skew = float(imputed[name].skew())
        model, transformed, lam, rescaled, detail = None, None, None, False, ""
        method = "standardise only"
        if abs(skew) > skew_threshold:
            model, transformed, lam, rescaled, detail = _fit_power(values)
            method = "yeo-johnson" if model is not None else "standardise only (power fit failed)"
        if model is None:
            model = StandardScaler().fit(values)
            transformed = model.transform(values).ravel()
        models[name], scaled[name] = model, transformed
        rows.append({"feature": name, "median": float(medians[name]),
                     "imputed_cells": int(original[name].isna().sum()),
                     "imputed_percent": float(100 * original[name].isna().mean()),
                     "method": method, "lambda": lam, "robust_rescale_first": rescaled,
                     "skew_before": skew, "skew_after": float(pd.Series(transformed).skew()),
                     "fit_note": detail})
    bundle = {"pca_input_columns": list(features), "medians": medians.to_dict(),
              "transformers": models, "skew_threshold": skew_threshold,
              "reference_samples": len(samples)}
    return imputed, scaled, bundle, pd.DataFrame(rows)


def transform_features(samples, bundle):
    """Apply stored training medians and transforms without refitting or using outcomes.

    Pass cleaned, sample-level numeric measurements in original units. New rows
    must satisfy the same quality rules, including the sparse-profile exclusion.
    """
    features = bundle["pca_input_columns"]
    missing = set(features) - set(samples)
    if missing:
        raise ValueError(f"Missing required features: {sorted(missing)}")
    values = samples[features].astype(float).fillna(bundle["medians"])
    if not np.isfinite(values.to_numpy()).all():
        raise ValueError("Unresolved missing or infinite inputs.")
    output = pd.DataFrame(index=samples.index)
    for name in features:
        output[name] = bundle["transformers"][name].transform(values[[name]].to_numpy()).ravel()
    if not np.isfinite(output.to_numpy()).all():
        raise ValueError("The fitted transform produced non-finite values on new data.")
    return output
